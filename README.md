# IsoAnalyzer
Analysis pipeline for alternative splicing events using bulk RNA-seq data


# Overview

This pipeline performs comprehensive RNA-Seq analysis, focusing on alignment, quantification, and the identification of alternative splicing events with functional consequence prediction. It integrates standard tools with parallelized optimization for efficiency.

# Installation and Usage

## Docker 
The recommended method for execution is using the provided Dockerfile available in the repository. This ensures all dependencies and environment configurations are correct. Alternatively, individual tools can be installed manually.

## Execution 
To launch the pipeline, ensure the main R script is located in the same directory as the settings file and the sample table. These files define the analysis parameters and the samples to be processed. Templates for these files are provided in the repository.

# Pipeline Workflow

## Preprocessing

The pipeline handles raw read processing based on the specific requirements of the dataset. UMI Handling: If Unique Molecular Identifiers (UMIs) are present, they are extracted and moved to the read header to facilitate downstream deduplication. This mode uses automatic sequence and adapter detection. Standard Trimming: For datasets without UMIs, TrimGalore is used for adapter removal, utilizing sequences provided in the settings file.

## Alignment and Deduplication

Reads are aligned to the reference genome using STAR. If UMIs are present, deduplication is performed using UMI-tools. This pipeline implements a custom parallelized variant of UMI-tools to exponentially increase analysis speed compared to the standard single-threaded version.

## Quantification

Gene-level counting is performed using featureCounts.

## Alternative Splicing Detection

Splicing events are identified using rMATS Turbo and SplAdder. Users can choose to run either program individually or both simultaneously.

### GTF Addition: When SplAdder is enabled, the pipeline offers a "GTF Addition" option. This generates an updated GTF file incorporating novel splicing events detected by SplAdder, which is then used for subsequent analysis steps.

##Isoform Quantification

Isoform-level quantification is handled by RSEM. This step requires alignment to the transcriptome rather than the genome. The pipeline can either accept pre-existing transcriptomic BAM files or generate them directly based on the provided settings.

## Data Filtering

Outputs from the detection tools undergo a filtering phase to retain only statistically significant splicing events.

## Functional Consequence Analysis

The final stage utilizes IsoformSwitchAnalyzeR to evaluate isoform switches and predict their functional impact on proteins. This analysis is divided into three distinct steps involving two different environments.

### Step 1: 
Preliminary Analysis Executed within the main pipeline/container. This step filters data, performs statistical testing, and reconstructs amino acid sequences.

### Step 2: 
External Protein Analysis This step requires the External Protein Analysis Bash script. Note: This script must be executed using the second Dockerfile provided (the sequence analysis container). It runs the following external tools to annotate the amino acid sequences generated in Step 1:

    CPAT (Coding Potential)

    Pfam (Protein Domains)

    SignalP (Signal Peptides)

    IUPred2A (Disordered Regions)
  
###Step 3: 
Final Integration Executed via the IsoformSwitchAnalyzer_Final.R script. Note: This script must be run using the primary Dockerfile. It requires the output directories from the External Protein Analysis (Step 2) as arguments to integrate the functional annotations and generate final visualizations.
