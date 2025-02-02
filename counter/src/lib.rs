use indicatif::{ProgressBar, ProgressStyle};
use kmer::{kmer::KmerGenerator, numeric_to_kmer, Kmer};
use ktio::seq::{get_reader, SeqFormat, Sequences};
use std::{
    cmp::{max, min},
    collections::{HashMap, HashSet},
    fs,
    io::{BufRead, BufReader, BufWriter, Read, Write},
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
    genome_ids: Arc<Mutex<Vec<String>>>,
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
            genome_ids: Arc::new(Mutex::new(Vec::new())),
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

    pub fn count(&mut self) {
        self.init();
        let pbar = ProgressBar::new(self.seq_count);
        pbar.set_style(
            ProgressStyle::with_template(
                "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} ({percent}%) {msg}",
            )
            .unwrap()
            .progress_chars("#>-"),
        );

        // Wrap pbar en Arc<Mutex<>> para compartir entre hilos
        let pbar_arc = Arc::new(Mutex::new(pbar));

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();

        let records_clone = Arc::clone(&self.records);
        let genome_ids_clone = Arc::clone(&self.genome_ids);
        let pbar_clone = Arc::clone(&pbar_arc); // Clone para el hilo

        pool.scope(|scope| {
            scope.spawn(move |_| {
                while let Some(record) = records_clone.lock().unwrap().next() {
                    pbar_clone.lock().unwrap().inc(1); // Acceso seguro al progress bar

                    // Usar record.id en lugar de record.header
                    let genome_id = record.id.split_whitespace().next().unwrap().to_string();
                    genome_ids_clone.lock().unwrap().push(genome_id);

                    let mut part_counts = vec![HashMap::new(); self.n_parts as usize];

                    for (fmer, rmer) in KmerGenerator::new(&record.seq, self.ksize) {
                        let min_mer = min(fmer, rmer);
                        let part = (min_mer % self.n_parts) as usize;
                        *part_counts[part].entry(min_mer).or_insert(0) += 1;
                    }

                    let chunk = self.chunks;
                    self.chunks += 1;

                    for (part, counts) in part_counts.into_iter().enumerate() {
                        let path = format!(
                            "{}/temp_kmers.part_{}_chunk_{}",
                            self.out_dir, part, chunk
                        );
                        let mut buff = BufWriter::new(fs::File::create(path).unwrap());
                        for (kmer, count) in counts {
                            writeln!(buff, "{}\t{}", count, kmer).unwrap();
                        }
                    }
                }
            });
        });

        pbar_arc.lock().unwrap().finish(); // Finalizar desde el hilo principal
    }

    pub fn merge(&self, delete: bool) {
        let genome_ids = self.genome_ids.lock().unwrap().clone();
        let outf = fs::File::create(format!("{}/kmers.counts", self.out_dir)).unwrap();
        let mut buff = BufWriter::new(outf);

        // Colectar todos los k-mers únicos
        let mut all_kmers = HashSet::new();
        for chunk in 0..genome_ids.len() {
            for part in 0..self.n_parts {
                let path = format!(
                    "{}/temp_kmers.part_{}_chunk_{}",
                    self.out_dir, part, chunk
                );
                if let Ok(file) = fs::File::open(&path) {
                    let buff = BufReader::new(file);
                    for line in buff.lines().filter_map(Result::ok) {
                        let kmer: Kmer = line.split('\t').nth(1).unwrap().parse().unwrap();
                        all_kmers.insert(kmer);
                    }
                }
            }
        }

        // Ordenar k-mers y escribir encabezado
        let mut sorted_kmers: Vec<Kmer> = all_kmers.into_iter().collect();
        sorted_kmers.sort_unstable();
        write!(buff, "ID_Genome").unwrap();
        for kmer in &sorted_kmers {
            let kmer_str = if self.acgt {
                numeric_to_kmer(*kmer, self.ksize)
            } else {
                kmer.to_string()
            };
            write!(buff, "\t{}", kmer_str).unwrap();
        }
        writeln!(buff).unwrap();

        // Escribir conteos por genoma
        for (chunk, genome_id) in genome_ids.iter().enumerate() {
            let mut genome_counts = HashMap::new();
            for part in 0..self.n_parts {
                let path = format!(
                    "{}/temp_kmers.part_{}_chunk_{}",
                    self.out_dir, part, chunk
                );
                if let Ok(file) = fs::File::open(&path) {
                    let buff = BufReader::new(file);
                    for line in buff.lines().filter_map(Result::ok) {
                        let mut parts = line.split('\t');
                        let count: u32 = parts.next().unwrap().parse().unwrap();
                        let kmer: Kmer = parts.next().unwrap().parse().unwrap();
                        genome_counts.insert(kmer, count);
                    }
                }
            }

            write!(buff, "{}", genome_id).unwrap();
            for kmer in &sorted_kmers {
                write!(buff, "\t{}", genome_counts.get(kmer).unwrap_or(&0)).unwrap();
            }
            writeln!(buff).unwrap();

            if delete {
                for part in 0..self.n_parts {
                    let path = format!(
                        "{}/temp_kmers.part_{}_chunk_{}",
                        self.out_dir, part, chunk
                    );
                    fs::remove_file(path).ok();
                }
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
        let mut ctr = CountComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_counts_test".to_owned(),
            15,
        );
        ctr.chunks = 2;
        ctr.n_parts = 2;
        ctr.merge(false);
        let exp = load_lines_sorted("../test_data/expected_counts_test.counts");
        let res = load_lines_sorted("../test_data/computed_counts_test/kmers.counts");
        println!("Result  : {:?}", res);
        println!("Expected: {:?}", exp);
        assert_eq!(exp, res);
    }

    #[test]
    fn merge_acgt_test() {
        let mut ctr = CountComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_counts_acgt_test".to_owned(),
            15,
        );
        ctr.chunks = 2;
        ctr.n_parts = 2;
        ctr.set_acgt_output(true);
        ctr.merge(false);
        let exp = load_lines_sorted("../test_data/expected_counts_acgt_test.counts");
        let res = load_lines_sorted("../test_data/computed_counts_acgt_test/kmers.counts");
        println!("Result  : {:?}", res);
        println!("Expected: {:?}", exp);
        assert_eq!(exp, res);
    }
}