# RBM47-driven Leukemic Cell-intrinsic Inflammation Shapes the Developmental Hierarchy of Human Acute Myeloid Leukemia

This repository contains analysis scripts and processed datasets corresponding to the primary figures of the research [RBM47-driven Leukemic Cell-intrinsic Inflammation Shapes the Developmental Hierarchy of Human Acute Myeloid Leukemia]()

The repository is structured to ensure transparency, reproducibility, and reusability of the computational analyses supporting this study.

## Overview

Acute myeloid leukemia (AML) exhibits a perturbed myeloid developmental hierarchy (DH), spanning from LSPCs to more differentiated monocytic states. In this work, we investigate how RBM47-mediated inflammatory signaling shapes this hierarchy, influences malignant cell-state composition, and modulates therapeutic response.

**This repository includes:**

1.  developmental trajectory modeling of AML differentiation

2.  inflammation-associated transcriptional profiling

3.  prioritization of RNA-binding protein (RBP) regulators

4.  eCLIP-defined RBM47 regulatory interactions

5.  drug-response screening

Scripts and resources are organized and documented by figure.

## Figure 1 — Inflammatory signaling progressively increases along the myeloid developmental hierarchy in AML

We conceptualize the leukemic LSPC → monocytic continuum as a pseudo-temporal developmental trajectory. Using the GSE232559 CITE-seq AML dataset, we applied Mfuzz clustering to identify distinct dynamic gene expression patterns along this hierarchy. Inflammatory pathway activity was then quantified across the trajectory.

Relevant resources Inflammation-related gene sets are provided in the repository: 
> ├─Data    
    └─GeneSet


To enable projection of bulk AML cohorts into developmental-hierarchy space, we generated a CIBERSORTx signature matrix derived from GSE232559. This matrix supports deconvolution of normalized RNA-seq datasets to infer leukemic developmental composition.

Bulk cohort deconvolution was performed using the [CIBERSORTx web portal](https://cibersortx.stanford.edu/)

with:
+ S-mode batch correction enabled
+ Absolute mode enabled
+ all other parameters at default settings

For downstream analyses, estimated abundances of the myeloid leukemic populations were normalized to sum to 1 per sample. Lymphoid and erythroid–megakaryocytic lineages were excluded.

Processed mixtures: 

[TCGA-LAML (with FAB annotation)](https://github.com/FionaMoon/Inflammation-Shapes-Developmental-Hierarchy-of-Human-Acute-Myeloid-Leukemia/blob/main/Data/CIBERSORTx/TCGA-LAML_Clinic_and_decon_mye.RData) 

[BEAT-AML (with FAB annotation)](https://github.com/FionaMoon/Inflammation-Shapes-Developmental-Hierarchy-of-Human-Acute-Myeloid-Leukemia/blob/main/Data/CIBERSORTx/BEAT-AML_decon_clinical_order.csv)

### Figure 2 — RBM47 emerges as a top-ranked RNA-binding protein associated with inflammation and AML hierarchy dynamics

Differential expression analyses were performed across high- versus low-inflammation AML groups in three independent datasets. RBM47 consistently ranked among the most strongly enriched inflammation-associated RBPs, identifying it as a putative master regulator of inflammatory and hierarchical cell-state programs.

### Figure 4 — RBM47 directly binds and regulates transcripts encoding inflammatory mediators

RBM47 target transcripts were defined using eCLIP-seq, processed end-to-end with the Skipper pipeline under default parameters. All preprocessing followed the standardized [YeoLab workflow](https://github.com/YeoLab/skipper); no additional custom preprocessing scripts are required.

### Figure 5 — TNFα / IL1β inflammatory signaling actively drives maturation of the AML developmental hierarchy

We evaluated the functional role of TNFα / IL1β in modulating hierarchy composition by CITE-seq. We combine 4 condition (CT, TNFα-Trt, IL1β-Trt and Combo) by HTO. Preprocessing procedures and analysis schemas are detailed in the manuscript. Representative downsampled datasets with annotations are provided [HRA015977_example]().

### Figure 7 — Targeting RBM47-dependent inflammatory signaling shifts AML toward a Venetoclax-sensitive primitive state

Drug-response screening demonstrated that: AML samples with high inflammatory signaling exhibit reduced Venetoclax sensitivity, and inhibition of signaling downstream of RBM47 reprograms hierarchical composition toward more primitive, Venetoclax-responsive cell states.

Together, these findings mechanistically link RBM47-dependent inflammation, AML hierarchical remodeling, and Venetoclax response.

**Notes** 
+ Additional documentation, extended analyses, and finalized identifiers will be released upon publication. 
+ Please cite the associated manuscript when using any component of this repository.
