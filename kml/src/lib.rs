use indicatif::{ProgressBar, ProgressStyle};
use kmer::{kmer::KmerGenerator, numeric_to_kmer, Kmer};
use ktio::seq::{get_reader, SeqFormat, Sequences};
use std::{
    cmp::{max, min},
    collections::{HashMap, HashSet},
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::Path,
    sync::{Arc, Mutex},
};

type SeqArc = Arc<Mutex<Sequences<BufReader<Box<dyn Read + Sync + Send>>>>>;

pub struct CountComputer {
    in_path: String,
    out_dir: String,
    ksize: usize,
    threads: usize,
    records: SeqArc,
    chunks: u64,
    n_parts: u64,
    memory_ceil_gb: f64,
    seq_count: u64,
    debug: bool,
    acgt: bool,
}

impl CountComputer {
    pub fn new(in_path: String, out_dir: String, ksize: usize) -> Self {
        let format = SeqFormat::get(&in_path).unwrap();
        let reader = ktio::seq::get_reader(&in_path).unwrap();
        let records = Sequences::new(format, reader).unwrap();

        Self {
            in_path,
            out_dir,
            ksize,
            threads: rayon::current_num_threads(),
            records: Arc::new(Mutex::new(records)),
            chunks: 0,
            n_parts: 0,
            seq_count: 0,
            memory_ceil_gb: 6_f64,
            debug: false,
            acgt: false,
        }
    }

    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads;
    }

    pub fn set_max_memory(&mut self, memory_ceil_gb: f64) {
        self.memory_ceil_gb = memory_ceil_gb;
    }

    pub fn set_acgt_output(&mut self, acgt: bool) {
        self.acgt = acgt;
    }

    // Setter para asignar un chunk único desde fuera
    pub fn set_chunk(&mut self, chunk: u64) {
        self.chunks = chunk;
    }

    pub fn count(&mut self) {
        self.init();
        let pbar = ProgressBar::new(self.seq_count);
        pbar.set_style(
            ProgressStyle::with_template(
                "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} ({percent}%) {msg}"
            )
            .unwrap()
            .progress_chars("#>-"),
        );

        let mut combined_seq = String::new();
        let mut label = String::new();
        {
            let mut records = self.records.lock().unwrap();
            // Procesamos el primer registro para extraer el label
            if let Some(record) = records.next() {
                let seq_str = std::str::from_utf8(&record.seq).unwrap();
                combined_seq.push_str(seq_str);
                // Se asume que record.description y record.id existen
                label = record.desc.replace(&record.id, "").trim().to_string();
                pbar.inc(1);
            }
            // Los registros restantes se concatenan sin modificar el label
            while let Some(record) = records.next() {
                let seq_str = std::str::from_utf8(&record.seq).unwrap();
                combined_seq.push_str(seq_str);
                pbar.inc(1);
            }
        }
        pbar.finish();

        let chunk = self.chunks;
        self.chunks += 1;

        // Se obtiene el genome_id (por ejemplo, a partir del nombre del archivo)
        let genome_id = Path::new(&self.in_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        // Guardar el ID del genoma y el label en genome_ids.txt (separados por tabulador)
        let genome_ids_path = format!("{}/genome_ids.txt", self.out_dir);
        let mut file = File::options()
            .create(true)
            .append(true)
            .open(genome_ids_path)
            .unwrap();
        writeln!(file, "{}\t{}", genome_id, label).unwrap();

        let mut part_counts = vec![HashMap::new(); self.n_parts as usize];
        for (fmer, rmer) in KmerGenerator::new(combined_seq.as_bytes(), self.ksize) {
            let min_mer = min(fmer, rmer);
            let part = (min_mer % self.n_parts) as usize;
            *part_counts[part].entry(min_mer).or_insert(0) += 1;
        }

        for (part, counts) in part_counts.into_iter().enumerate() {
            let path = format!("{}/temp_kmers.part_{}_chunk_{}", self.out_dir, part, chunk);
            let mut buff = BufWriter::new(File::create(path).unwrap());
            for (kmer, count) in counts {
                writeln!(buff, "{}\t{}", count, kmer).unwrap();
            }
        }
    }

    fn init(&mut self) {
        let reader = get_reader(&self.in_path).unwrap();
        let format = SeqFormat::get(&self.in_path).unwrap();
        let stats = Sequences::seq_stats(format, reader);
        let data_size_gb = stats.total_length as f64 / (1 << 30) as f64;
        let n_parts = max(
            if self.debug { 1 } else { self.threads as u64 },
            (8_f64 * data_size_gb / (2_f64 * self.memory_ceil_gb)).ceil() as u64,
        );
        self.n_parts = n_parts;
        self.seq_count = stats.seq_count as u64;
    }
}

// Función independiente para fusionar todos los archivos temporales y generar el archivo final
pub fn merge_all(out_dir: &str, ksize: usize, acgt: bool, delete: bool) {
    // Leer el archivo genome_ids.txt, ahora con dos columnas: genome_id y label.
    let genome_ids_path = format!("{}/genome_ids.txt", out_dir);
    let genome_records = match fs::read_to_string(&genome_ids_path) {
        Ok(content) => {
            content
                .lines()
                .map(|line| {
                    let mut parts = line.split('\t');
                    let genome_id = parts.next().unwrap_or("").to_string();
                    let label = parts.next().unwrap_or("").to_string();
                    (genome_id, label)
                })
                .collect::<Vec<_>>()
        }
        Err(_) => Vec::new(),
    };

    let mut all_kmers = HashSet::new();
    // Recorremos cada chunk (fila) usando el número de líneas en genome_ids.txt
    for chunk in 0..genome_records.len() {
        let mut part = 0;
        loop {
            let path = format!("{}/temp_kmers.part_{}_chunk_{}", out_dir, part, chunk);
            if !Path::new(&path).exists() {
                break;
            }
            let file = File::open(&path).unwrap();
            let reader = BufReader::new(file);
            for line in reader.lines().filter_map(Result::ok) {
                let kmer: Kmer = line.split('\t').nth(1).unwrap().parse().unwrap();
                all_kmers.insert(kmer);
            }
            part += 1;
        }
    }

    let mut sorted_kmers: Vec<Kmer> = all_kmers.into_iter().collect();
    sorted_kmers.sort_unstable();

    let outf = File::create(format!("{}/kmers.counts", out_dir)).unwrap();
    let mut buff = BufWriter::new(outf);

    // Escribir el header con las columnas: file, label, y luego los k-mers
    write!(buff, "file\tlabel").unwrap();
    for kmer in &sorted_kmers {
        let kmer_str = if acgt {
            numeric_to_kmer(*kmer, ksize)
        } else {
            kmer.to_string()
        };
        write!(buff, "\t{}", kmer_str).unwrap();
    }
    writeln!(buff).unwrap();

    // Escribir los conteos por genoma, incluyendo genome_id y label
    for (chunk, (genome_id, label)) in genome_records.iter().enumerate() {
        let mut genome_counts = HashMap::new();
        let mut part = 0;
        loop {
            let path = format!("{}/temp_kmers.part_{}_chunk_{}", out_dir, part, chunk);
            if !Path::new(&path).exists() {
                break;
            }
            let file = File::open(&path).unwrap();
            let reader = BufReader::new(file);
            for line in reader.lines().filter_map(Result::ok) {
                let mut parts = line.split('\t');
                let count: u32 = parts.next().unwrap().parse().unwrap();
                let kmer: Kmer = parts.next().unwrap().parse().unwrap();
                *genome_counts.entry(kmer).or_insert(0) += count;
            }
            part += 1;
        }

        write!(buff, "{}\t{}", genome_id, label).unwrap();
        for kmer in &sorted_kmers {
            write!(buff, "\t{}", genome_counts.get(kmer).unwrap_or(&0)).unwrap();
        }
        writeln!(buff).unwrap();

        if delete {
            let mut part = 0;
            loop {
                let path = format!("{}/temp_kmers.part_{}_chunk_{}", out_dir, part, chunk);
                if fs::remove_file(&path).is_err() {
                    break;
                }
                part += 1;
            }
        }
    }

    if delete {
        fs::remove_file(genome_ids_path).ok();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ktio::fops::{create_directory, load_lines_sorted};

    const PATH_FQ: &str = "../test_data/reads.fq";

    #[test]
    fn count_test() {
        create_directory("../test_data/computed_counts").expect("Directory must be creatable");
        let mut ctr = CountComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_counts".to_owned(),
            15,
        );
        ctr.debug = true;
        ctr.count();
        assert_eq!(ctr.n_parts, 1);
        assert_eq!(ctr.chunks, 1);
        let exp = load_lines_sorted("../test_data/expected_counts.part_0_chunk_0");
        let res = load_lines_sorted("../test_data/computed_counts/temp_kmers.part_0_chunk_0");
        println!("Result  : {:?}", res);
        println!("Expected: {:?}", exp);
        assert_eq!(exp, res);
    }

    #[test]
    fn merge_test() {
        let out_dir = "../test_data/computed_counts_test";
        let ksize = 15;
        let mut ctr = CountComputer::new(PATH_FQ.to_owned(), out_dir.to_owned(), ksize);
        ctr.debug = true;
        ctr.count();
        merge_all(out_dir, ksize, false, false);
        // Verificar el archivo kmers.counts
    }
}
