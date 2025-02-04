```markdown
# kmertools: DNA Vectorisation Tool (Fork with New KML Function)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Cargo tests](https://github.com/anuradhawick/kmertools/actions/workflows/rust_test.yml/badge.svg)](https://github.com/anuradhawick/kmertools/actions/workflows/rust_test.yml)
[![Clippy check](https://github.com/anuradhawick/kmertools/actions/workflows/clippy_check.yml/badge.svg)](https://github.com/anuradhawick/kmertools/actions/workflows/clippy_check.yml)
[![codecov](https://codecov.io/gh/anuradhawick/kmertools/graph/badge.svg?token=IDGRE54SSQ)](https://codecov.io/gh/anuradhawick/kmertools)

<div align="center">
<pre>
$$\   $$\                                   $$$$$$$$\                     $$\           
$$ | $$  |                                  \__$$  __|                    $$ |          
$$ |$$  / $$$$$$\$$$$\   $$$$$$\   $$$$$$\     $$ |    $$$$$$\   $$$$$$\  $$ | $$$$$$$\ 
$$$$$  /  $$  _$$  _$$\ $$  __$$\ $$  __$$\    $$ |   $$  __$$\ $$  __$$\ $$ |$$  _____|
$$  $$<   $$ / $$ / $$ |$$$$$$$$ |$$ |  \__|   $$ |   $$ /  $$ |$$ /  $$ |$$ |\$$$$$$\  
$$ |\$$\  $$ | $$ | $$ |$$   ____|$$ |         $$ |   $$ |  $$ |$$ |  $$ |$$ | \____$$\ 
$$ | \$$\ $$ | $$ | $$ |\$$$$$$$\ $$ |         $$ |   \$$$$$$  |\$$$$$$  |$$ |$$$$$$$  |
\__|  \__|\__| \__| \__| \_______|\__|         \__|    \______/  \______/ \__|\_______/ 
</pre>
</div>

## Overview

This repository is a **fork** of the original [Kmertools](https://github.com/anuradhawick/kmertools) project that adds a new function: **KML**.  
The **KML** function generates a table from DNA vectorisation data, streamlining the creation of datasets ready for machine learning applications. All other original features of Kmertools remain intact, including:

- **Oligonucleotide Frequency Vectors:** Generate frequency vectors for oligonucleotides.
- **Minimiser Binning:** Efficiently bin sequences using minimisers to reduce data complexity.
- **Chaos Game Representation (CGR):** Compute CGR vectors for DNA sequences using k-mers or full sequence transformation.
- **Coverage Histograms:** Create coverage histograms to analyze sequencing read depth.
- **Python Binding:** Import kmertools functionality using `import pykmertools as kt`.

## Installation from Sources

To install `kmertools` (including the new KML functionality) from source, follow these steps. Make sure you have [Rust](https://www.rust-lang.org/tools/install) and its package manager `cargo` installed.

```bash
# Clone the repository
git clone https://github.com/your-repository/kmertools.git
cd kmertools

# Build the project in release mode
cargo build --release
```

Add the generated binary to your PATH (you can modify your `~/.bashrc` or `~/.zshrc`):

```sh
# To add it to the current terminal session:
export PATH=$PATH:$(pwd)/target/release/

# To add it permanently (for Bash):
echo "export PATH=\$PATH:$(pwd)/target/release/" >> ~/.bashrc
source ~/.bashrc

# For Zsh (e.g., on macOS):
echo "export PATH=\$PATH:$(pwd)/target/release/" >> ~/.zshrc
source ~/.zshrc
```

## Testing the Installation

After installation, verify everything is set up correctly by running:

```bash
kmertools --help
```

This command should display the help message, confirming that the program (including the new KML functionality) is working as expected.
