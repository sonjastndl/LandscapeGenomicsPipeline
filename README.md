# Landscape Genomics Pipeline 

This repository is supposed to hold documentation on the Pipeline used for the UC3 Project of FAIRiCUBE. 

## Objectives

This pipeline performs several landscape-genomics analysis on a SNP called data set. 

## Preliminary Results on Test Data

---

## Requirements
- VCFTools (Citation)
- Baypass2.3 (Citation)
- LEA R-Package (Citation)

## Workflow 

0) Download Data
1) Set Up Environment
2) Subsample VCF for chromosomal arms (computational reasons while pipeline development)
2.1) Remove Polyploides
2.2) Subsample 10k Variants
2.3) Convert to Allele Frequencies
3) Define Metadata / Covariates for analyses

4) Linear Model 

5) Latent Factor Mixed Model

6) Baypass Analysis

7) Result Comparison 

x) The basic workflow can be performed according to [main.sh](d/d/main.sh)

---

## How to execute

