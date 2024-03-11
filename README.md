# GGOT: Detecting Disease Critical Transitions Using Gaussian Graphical Optimal Transport Embedded With Protein Interaction Networks.

This repository contains our proposed methods *Gaussian Graphical Optimal Transport* (GGOT) for detecting critical transitions and identifying trigger molecules of diseases, the corresponding pre-processed diseases datasets from TCGA and GEO databases, and pre-processed Protein-Protein Interaction (PPI) Networks of the Human and Mouse Genomes.

Code for the paper:

> Wenbo Hua, Ruixia Cui, Heran Yang, Jingyao Zhang, Chang Liu, Jian Sun. "Detecting Disease Critical Transitions Using Gaussian Graphical Optimal Transport Embedded With Protein Interaction Networks"

<!-- [[arxiv]](https://arxiv.org/abs/1907.03907) -->

## Overview

GGOT uses Gaussian graphical model, incorporating the gene interaction network, to model the data distributions at different disease stages. Then we use population-level optimal transport to calculate the Wasserstein distance and transport maps between stages, enabling us to detect critical transitions. By analyzing the per-molecule transport map, we quantify the importance of each molecule and identify the trigger molecules. Moreover, GGOT predicts the occurrence of critical transitions in unseen samples and visualizes the disease progression process.

<img src="assets/Overview.png" alt="Overview" style="zoom: 25%;" />

## Datasets

Clicking on the name of the corresponding dataset will redirect you to the website to download the corresponding dataset. The dataset of **XJTUSepsis** needs to contact hwb0856@stu.xjtu.edu.cn to get access.

| Dataset                                                                | Description               | Species      | stages | Type of disease                       |
| ---------------------------------------------------------------------- | ------------------------- | ------------ | ------ | ------------------------------------- |
| [GSE48452](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48452)   | non-alcoholic fatty liver | Home sapiens | 4      | Chronic progressive benign disease    |
| [GSE2565](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2565)     | lung injury               | Mus musculus | 10     | Acute progressive noncritical disease |
| [LUAD](https://portal.gdc.cancer.gov/projects/TCGA-LUAD)                  | lung adenocarcinoma       | Home sapiens | 8      | Chronic progressive malignant disease |
| [COAD](https://portal.gdc.cancer.gov/projects/TCGA-COAD)                  | colon adenocarcinoma      | Home sapiens | 8      | Chronic progressive malignant disease |
| XJTUSepsis                                                             | sepsis                    | Home sapiens | 8      | Acute progressive critical disease    |
| [GSE154918](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154918) | sepsis                    | Home sapiens | 4      | Acute progressive critical disease    |

## Prerequisites

1 Install required packages

```bash
pip install -r requirements.txt
```

2 Download the pre-processed Protein-Protein Interaction (PPI) Networks of the Human and Mouse Genomes into "`PPI/`". We construct the PPI networks using [string](https://string-db.org) database.

PPI Networks of two species:[[GGOT_PPI](https://drive.google.com/file/d/1Gxp5MpbQQ3l4wRtxBOYSaGVxdYS6e7Xy/view?usp=sharing)]

3 Download the RNA-seq dataset of different diseases into "`data/source/`".

Raw datasets: [[GGOT_datasets](https://drive.google.com/file/d/1c0SqU3dq22lE5qNlW7wbKG5UqoHKS83_/view?usp=sharing)]

**Make sure your data file tree is similar to the following.**

```
GGOT
├──run_model.py
├──data/
│  ├──source
│  │  ├──GSE48452
│  │  │  ├──GSE48452.csv (RNA-seq expression)
│  │  │  ├──group.csv (stage information)
│  │  ├──GSE2565
│  │  │  ├──GSE48452.csv
│  │  │  ├──group.csv
│  │  ├──LUAD
│  │  │  ├──LUAD.csv
│  │  │  ├──group.csv
│  │  ├──......
├──PPI
│  ├──Human
│  │  ├──ppi_database.npy
│  │  ├──map_gene_protein.csv
│  │  ├──map_gene_protein_full.csv
│  │  ├──symbol2id.csv
│  ├──Mus
│  │  ├──ppi_database.npy
│  │  ├──map_gene_protein.csv
│  │  ├──map_gene_protein_full.csv
│  │  ├──symbol2id.csv
├──......
```

## Running GGOT for datasets with different types

- running model for the chronic progressive non-critical disease, GSE48452

```bash
python run_model.py -d GSE48452 -s Human
```

- running model for the acute progressive non-critical disease, GSE2565

```bash
python run_model.py -d GSE2565 -s Mus
```

- running model for the chronic progressive critical disease, LUAD, COAD

```bash
python run_model.py -d LUAD -s Human
python run_model.py -d COAD -s Human
```

- running model for the acute progressive critical disease, XJTUSepsis, GSE154918

```bash
python run_model.py -d XJTUSepsis -s Human
python run_model.py -d GSE154918 -s Human
```

## Making the visualization

```bash
python visualization.py -d GSE2565
```

## Results (example on simulation dataset)

<img src="assets/Numsim.png" alt="Overview" style="zoom:25%;">

If you have any problems with this repository, please send emails to hwb0856@stu.xjtu.edu.cn for discussions.
