# Overcoming selection bias in synthetic lethality prediction
[Link to paper](https://doi.org/10.1093/bioinformatics/btac523)


==============================

## Abstract

**Motivation:** Synthetic lethality (SL) between two genes occurs when simultaneous loss of function leads to cell death. This holds great promise for developing anti-cancer therapeutics that target synthetic lethal pairs of endogenously disrupted genes. Identifying novel SL relationships through exhaustive experimental screens is challenging, due to the vast number of candidate pairs. Computational SL prediction is therefore sought to identify promising SL gene pairs for further experimentation. However, current SL prediction methods lack consideration for generalizability in the presence of selection bias in SL data.

**Results:** We show that SL data exhibit considerable gene selection bias. Our experiments designed to assess the robustness of SL prediction reveal that models driven by the topology of known SL interactions (e.g. graph, matrix factorization) are especially sensitive to selection bias. We introduce selection bias-resilient synthetic lethality (SBSL) prediction using regularized logistic regression or random forests. Each gene pair is described by 27 molecular features derived from cancer cell line, cancer patient tissue and healthy donor tissue samples. SBSL models are built and tested using approximately 8000 experimentally derived SL pairs across breast, colon, lung and ovarian cancers. Compared to other SL prediction methods, SBSL showed higher predictive performance, better generalizability and robustness to selection bias. Gene dependency, quantifying the essentiality of a gene for cell survival, contributed most to SBSL predictions. Random forests were superior to linear models in the absence of dependency features, highlighting the relevance of mutual exclusivity of somatic mutations, co-expression in healthy tissue and differential expression in tumour samples.

This repository contains all the code used to process the data and generate the results.


Directory Structure
------------

    ├── LICENSE
    ├── r                                <- R scripts 
    │   └── data                         <- processed data files  
    │   └── experiments                  <- R scripts for experiment analysis
    │   │   ├── 0.Dataset Analysis       <- R scripts for analysing dataset labels
    │   │   ├── 1.Performance            <- R scripts for different training and testing scenarios
    │   │   │   ├── 1.1 Base Datasets            
    │   │   │   ├── 1.2 Per Cancer               
    │   │   │   ├── 1.3 Dataset Cross Comparison 
    │   │   │   ├── 1.4 Ensemble                 
    │   │   │   ├── 1.5 Embeddings               
    │   │   │   ├── 1.6 Training Curves 
    │   │   │   ├── 1.7 Cross Cancer 
    │   │   │   ├── 1.8 Dep Scores
    │   │   │   └── 1.9 Gene dropout  
    │   │   └── 2.Feature Analysis       <- R scripts for feature analysis
    │   │       ├── 2.Feature Importance 
    │   │       └── 2.Feature Selection  
    │   └── features                     <- R scripts for generating model features
    │   └── util                         <- R scripts for useful shared functions
    │
    ├── results                          <- data and scripts to produce paper figures
    │   └── figures                      <- paper figures
    │   └── scripts                      <- R scripts to produce the figures
    │   
    ├── raw_data                         <- raw data
        └── labels                       <- SL labels
--------
