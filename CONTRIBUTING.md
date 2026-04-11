# Contributing to BDB Genomics Pipelines

First off, thank you for considering contributing to BDB Genomics! It's people like you that make open-source bioinformatics powerful, robust, and accessible.

Our goal is to build the most rigorous, completely modular standard for epigenomic and transcriptomic data processing. We welcome contributions from the community--whether you are fixing a bug, adding a new experimental utility, or improving our documentation.

## How Can I Contribute?

### 1\. Reporting Bugs & Suggesting Enhancements 
If you encounter a crashed rule, a missing conda dependency, or just have a great idea for a new feature (like integrating a new peak caller), please check if an issue already exists.
- If not, open a new **Issue**!
- - Provide as much detail as possible, including snakemake logs, environment specs, and a minimum reproducible fastq example if relevant.
 
  - ### 2\. Submitting Pull Requests (PRs)
  - We love pull requests! If you want to contribute code:
  - 1\. **Fork** the repository and create your branch from `main`.
  - 2\. **Develop** your isolated Snakemake rule (`.smk`) or script.
  - 3\. **Use Conda**: Ensure your rule relies strictly on a deterministic conda environment descriptor.
  - 4\. **Test**: Run your changes locally against a small dataset. Ensure it passes the automatic configuration validation (`validate_config.py`).
  - 5\. **Open a PR**: Describe what your code does, why it is needed, and any parameters you added.
 
  - ### 3\. Improving Documentation
  - Bioinformatics is complex. If you wrote a great tutorial on how to run our pipelines on a specific HPC setting (like Slurm or LSF), please submit a documentation PR.
 
  - *Thank you for helping us build better science.*
