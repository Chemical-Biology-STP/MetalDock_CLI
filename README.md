# MetalDock CLI

MetalDock CLI is a command-line interface for setting up and running metal coordination docking workflows, supporting both single and batch ligand screening.

## 📥 Installation Instructions

Follow these steps carefully, even if you're new to command-line tools:

### 1. Clone the Repository

First, make sure you have Git installed. Then open a terminal and run:

```bash
git clone git@github.com:Chemical-Biology-STP/MetalDock_CLI.git
cd MetalDock_CLI
```

### 2. Install Conda (if not already installed)

Download and install Miniconda (recommended for beginners) from:  
👉 [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

Follow the installation instructions for your OS.

### 3. Create and Activate the Conda Environment

```bash
conda env create -f environment.yml
conda activate metaldock_cli
```

### 4. Install the CLI Tool

From within the `MetalDock_CLI` directory:

```bash
pip install .
```

You can now run the CLI with:

```bash
metaldock_cli
```

## 🚀 Usage

The tool will guide you step-by-step through either a **single-ligand** or **batch-ligand** docking setup. You will need:

- A protein structure in `.pdb` format.
- One or more ligands in `.xyz` format.
- ORCA and MetalDock binaries in your system `$PATH`.

> The script will prompt you to enter all relevant docking parameters interactively, and generate input files and SBATCH scripts for HPC submission.

## ❗ System Requirements

- ORCA (quantum chemistry software)
- MetalDock installed and accessible in `$PATH`
- SLURM (for SBATCH job submission)

If either ORCA or MetalDock is missing, the script will notify you with instructions on how to get help.

## 📖 Citation

If you use this code in any context, please cite:

> Y. M. Yip (2025). *MetalDock CLI* (v0.1.0) [Computer software].  
> DOI: [<zenodo_doi>](https://doi.org/<zenodo_doi>)
