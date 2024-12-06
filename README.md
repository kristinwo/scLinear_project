# Evaluating scLinear's Performance in Cancerous Tissues

## Project Overview
This project investigates the performance of the **scLinear** model from the R package in predicting **antibody-derived tags (ADT)** from **gene expression (GEX)** data in cancerous tissues. While scLinear has demonstrated high accuracy in immune-related cells from healthy donors, this study aims to explore its prediction in kidney, breast and lung cancer.

The research is based on the work of **Daniel Hanhart**, whose method is outlined in the paper [ScLinear predicts protein abundance at single-cell resolution](https://www.nature.com/articles/s42003-024-05958-4). This study expands the application of scLinear beyond its initial focus on healthy hematopoietic cells.

## Supervisor
The project is conducted under the supervision of **Panagiotis Chouvardas**.

## About scLinear
The **scLinear** package provides a streamlined solution for predicting ADT data from single-cell RNA-sequencing (scRNA-seq) data. Equipped with pre-trained models and customizable options for training new models, scLinear simplifies preprocessing and modeling of ADT data, making it accessible to a wide range of users.

According to the package description:

> "The goal of scLinear is to predict antibody-derived tags (ADT) data from gene expression data in scRNA-seq data. It includes all the necessary pre-processing steps, comes equipped with pre-trained models and also allows the training of new models."

For more information, visit the scLinear GitHub page: [https://github.com/DanHanh/scLinear](https://github.com/DanHanh/scLinear)

## Research Aim
The project seeks to fill a gap in the current usage of scLinear, as the paper notes:

> "The measured proteins in the available CITE-seq datasets are, in their vast majority, immune-related. Therefore, although scLinear shows very accurate results in hematopoietic cells from healthy donors, its performance in pathological conditions and in modeling the abundance of non-immune protein markers remain to be explored."

Given that scLinear's efficacy in non-immune cells like cancer remains untested, the core objective is to determine whether scLinear can successfully predict ADT data for **cancerous tissues**.

## References
- **Daniel Hanhart**, scLinear: [GitHub Repository](https://github.com/DanHanh/scLinear)
- [Research Paper](https://www.nature.com/articles/s42003-024-05958-4)
