use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use kmer::{kmer::KmerGenerator, numeric_to_kmer, Kmer};
use ktio::seq::{get_reader, SeqFormat, Sequences};

pub type GenomeResult = (String, String, Vec<HashMap<Kmer, u32>>);

pub struct CountComputer {
    in_path: String,
    ksize: usize,
    threads: usize,
    chunks: u64,
    n_parts: u64,
    memory_ceil_gb: f64,
    seq_count: u64,
    debug: bool,
    acgt: bool,
    temp_dir: Option<String>, // New field for temporary directory
}

impl CountComputer {
    // Creates a new instance.
    pub fn new(in_path: String, _out_dir: String, ksize: usize) -> Self {
        let format = SeqFormat::get(&in_path).expect("Unable to determine sequence format");
        let reader = get_reader(&in_path).expect("Unable to open file");
        // Create an iterator to calculate statistics only
        let _ = Sequences::new(format, reader).expect("Error creating sequences");
        Self {
            in_path,
            ksize,
            threads: rayon::current_num_threads(),
            chunks: 0,
            n_parts: 0,
            seq_count: 0,
            memory_ceil_gb: 6_f64,
            debug: false,
            acgt: false,
            temp_dir: None, // Initialize as None
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

    /// Allows setting a chunk number (useful for processing multiple executions)
    pub fn set_chunk(&mut self, chunk: u64) {
        self.chunks = chunk;
    }

    /// Set a custom temporary directory for intermediate files
    pub fn set_temp_dir(&mut self, temp_dir: String) {
        self.temp_dir = Some(temp_dir);
    }

    /// Get the directory to use for temporary files
    pub fn get_temp_dir(&self) -> &str {
        self.temp_dir.as_deref().unwrap_or("")
    }

    /// --- New method: count_with_chunks ---
    /// Processes sequences sequentially, accumulating counts in RAM until the limit is reached.
    /// When the limit is reached, the current chunk is flushed (saved) and the counters reset.
    /// Returns a vector of GenomeResult (one per chunk).
    pub fn count_with_chunks(&mut self) -> Vec<GenomeResult> {
        self.init(); // Initializes n_parts and seq_count (although seq_count is no longer used for progress)
        let mut chunks_results = Vec::new();

        // Obtain genome_id from the file name.
        let genome_id = Path::new(&self.in_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        // Initialize label (set from the first sequence)
        let mut label: Option<String> = None;

        // Initialize counters for this chunk: a vector with one HashMap per partition.
        let mut part_counts: Vec<HashMap<Kmer, u32>> =
            (0..self.n_parts).map(|_| HashMap::new()).collect();

        // Counter for bases accumulated in the current chunk.
        let mut total_bases: u64 = 0;
        // Convert the memory limit to number of bases (assuming 8 bytes per k-mer)
        let mem_limit_bases = (1_000_000_000_f64 * self.memory_ceil_gb / 8.0) as u64;

        // Create the sequences iterator (sequentially)
        let reader = get_reader(&self.in_path).expect("Unable to get reader");
        let format = SeqFormat::get(&self.in_path).expect("Unable to determine format");
        let mut sequences = Sequences::new(format, reader).expect("Error creating sequences");

        // Iterate over the sequences.
        while let Some(sequence) = sequences.next() {
            // Assign label if not already assigned.
            if label.is_none() {
                label = Some(if sequence.desc.trim().is_empty() {
                    "".to_string()
                } else {
                    sequence.desc.replace(&sequence.id, "").trim().to_string()
                });
            }
            let seq_len = sequence.seq.len() as u64;
            total_bases += seq_len;
            // Process the k-mers of the sequence and accumulate in the corresponding partition.
            for (fmer, rmer) in KmerGenerator::new(sequence.seq.as_slice(), self.ksize) {
                let min_mer = std::cmp::min(fmer, rmer);
                let part = (min_mer % self.n_parts as Kmer) as usize;
                *part_counts[part].entry(min_mer).or_insert(0) += 1;
            }
            // If the memory limit for this chunk is reached or exceeded, flush the chunk.
            if total_bases >= mem_limit_bases {
                chunks_results.push((
                    genome_id.clone(),
                    label.clone().unwrap_or_default(),
                    part_counts,
                ));
                // Reset counters for the next chunk.
                part_counts = (0..self.n_parts).map(|_| HashMap::new()).collect();
                total_bases = 0;
            }
        }
        // If there are remaining data that have not been flushed, add them as the last chunk.
        if total_bases > 0 {
            chunks_results.push((
                genome_id,
                label.unwrap_or_default(),
                part_counts,
            ));
        }
        chunks_results
    }

    /// Initializes statistics and configures `n_parts` and `seq_count`.
    fn init(&mut self) {
        let reader = get_reader(&self.in_path).expect("Unable to get reader in init");
        let format = SeqFormat::get(&self.in_path).expect("Unable to determine format in init");
        let stats = Sequences::seq_stats(format, reader);
        let data_size_gb = stats.total_length as f64 / ((1 << 30) as f64);
        self.n_parts = std::cmp::max(
            if self.debug { 1 } else { self.threads as u64 },
            ((8_f64 * data_size_gb) / (2_f64 * self.memory_ceil_gb)).ceil() as u64,
        );
        self.seq_count = stats.seq_count as u64;
    }
}

/// Merges all in-memory accumulated results and generates the final `kmers.counts` file.
/// Accepts a vector of GenomeResult and writes the final output.
/// This function processes the already flushed chunks.
pub fn merge_all(
    genome_results: Vec<GenomeResult>,
    out_dir: &str,
    ksize: usize,
    acgt: bool,
) {
    // 1. Globally merge all k-mers and counts per genome.
    let mut global_kmers = HashSet::new();
    let mut merged_genome_counts = Vec::new();

    for (_genome_id, _label, part_counts) in &genome_results {
        let mut genome_map: HashMap<Kmer, u32> = HashMap::new();
        for part in part_counts {
            for (&kmer, &count) in part.iter() {
                *genome_map.entry(kmer).or_insert(0) += count;
                global_kmers.insert(kmer);
            }
        }
        merged_genome_counts.push(genome_map);
    }

    // 2. Globally sort the k-mers.
    let mut sorted_kmers: Vec<Kmer> = global_kmers.into_iter().collect();
    sorted_kmers.sort_unstable();

    // 3. Write the final counts file.
    let out_file_path = format!("{}/kmers.counts", out_dir);
    let outf = File::create(&out_file_path)
        .expect("Error creating kmers.counts file");
    let mut buff = BufWriter::new(outf);

    // Write the header: "file\tlabel" followed by each k-mer.
    write!(buff, "file\tlabel").expect("Error writing header");
    for kmer in &sorted_kmers {
        let kmer_str = if acgt {
            numeric_to_kmer(*kmer, ksize)
        } else {
            kmer.to_string()
        };
        write!(buff, "\t{}", kmer_str).expect("Error writing k-mer in header");
    }
    writeln!(buff).expect("Error finishing header line");

    // 4. For each chunk (genome), write a line with its id, label, and counts for each k-mer.
    for (i, (genome_id, label, _)) in genome_results.iter().enumerate() {
        write!(buff, "{}\t{}", genome_id, label)
            .expect("Error writing genome_id and label");
        let genome_counts = &merged_genome_counts[i];
        for kmer in &sorted_kmers {
            let count = genome_counts.get(kmer).unwrap_or(&0);
            write!(buff, "\t{}", count)
                .expect("Error writing count for k-mer");
        }
        writeln!(buff).expect("Error finishing genome count line");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ktio::fops::create_directory;
    use std::path::Path;

    const PATH_FQ: &str = "../test_data/reads.fq";

    #[test]
    fn count_test() {
        // Create the output directory for the test.
        create_directory("../test_data/computed_counts")
            .expect("Directory must be creatable");
        let mut ctr = CountComputer::new(
            PATH_FQ.to_owned(),
            "../test_data/computed_counts".to_owned(),
            15,
        );
        ctr.debug = true;
        let chunks = ctr.count_with_chunks();
        // Verify that at least one chunk is generated and each chunk has counts.
        assert!(!chunks.is_empty(), "There should be at least one chunk");
        for (genome_id, _label, part_counts) in chunks {
            assert!(!genome_id.is_empty(), "Genome id should not be empty");
            assert!(!part_counts.is_empty(), "Part counts vector should not be empty");
        }
    }

    #[test]
    fn merge_test() {
        let out_dir = "../test_data/computed_counts_test";
        create_directory(out_dir).expect("Directory must be creatable");
        let ksize = 15;
        let mut ctr = CountComputer::new(PATH_FQ.to_owned(), out_dir.to_owned(), ksize);
        ctr.debug = true;
        let genome_results = ctr.count_with_chunks();
        merge_all(genome_results, out_dir, ksize, false);
        // Verify that the final kmers.counts file is generated.
        let kmers_counts_path = format!("{}/kmers.counts", out_dir);
        assert!(
            Path::new(&kmers_counts_path).exists(),
            "kmers.counts file should exist"
        );
    }
}
