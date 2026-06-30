# Characterization of the Functional Landscape of Somatic Evolution using Protein Language Models

Master Thesis project to earn the Master's Degree in Computational Biology from the Universidad Politécnica de Madrid.

> **⚠️ Important Note:**
> Due to a major incident on June 2026 affecting HPC cluster used for this project, several scripts and supplementary outputs referenced in the thesis are temporarily unavailable. These will be uploaded to this repository as soon as cluster access is restored.

---

## Contents

### 1. Systematic pLM Comparison
Scripts implementing the high-throughput benchmark of six protein language models (ProtT5, ProstT5, ESMv1, VESM 650M, ESM2, ESMC) on an _in-silico_ saturation mutagenesis experiment across oncoKB cancer genes (Bandlamudi et al., 2026), including embedding extraction, deviation profile computation, and statistical association with functional annotations and AlphaMissense pathogenicity scores.

### 2. SAE-Based Functional Interpretation and Human-SAE training
Analysis pipeline for transforming pLM embeddings into interpretable sparse representations using the [InterPLM](https://github.com/ElanaPearl/InterPLM/) toolkit, including validation against oncoKB/MSK-50K pathogenicity annotations and human proteome concept-feature mapping.

### 3. plm_dissect
A command-line tool for the automated extraction and analysis of SAE-derived interpretable features from wild-type and mutant protein sequences.

**Installation:**
```bash
git clone https://github.com/lcanoalmarza/Characterization-of-the-Functional-Landscape-of-Somatic-Evolution-using-Protein-Language-Models
cd plm-dissect_project
pip install -e .
```

**Usage:** Given a multifasta file containing wild-type and mutant sequences, plm_dissect provides four analytical modules: `score` (disruption quantification), `visualize` (functional interpretation of deviating features), `cluster` (variant clustering by disruption profile similarity), and `profile` (comparison against oncoKB-annotated variants).

Example analyses demonstrating the tool's functionality are currently stored on the affected HPC cluster and could not be included here.

---

## Contact

For questions regarding this repository, please contact laura.cano@alumnos.upm.es
