# IsoAnalyzer
Analysis pipeline for alternative splicing events using bulk RNA-seq data


# Overview

This pipeline performs comprehensive RNA-Seq analysis, focusing on alignment, quantification, and the identification of alternative splicing events with functional consequence prediction. It integrates standard tools with parallelized optimization for efficiency.

## Table of contents

## Table of Contents

- [Overview](#overview)
- [Installation and Usage](#installation-and-usage)
  * [Docker](#docker)
  * [Execution](#execution)
- [Pipeline Workflow](#pipeline-workflow)
  * [Preprocessing](#preprocessing)
  * [Alignment and Deduplication](#alignment-and-deduplication)
  * [Quantification](#quantification)
  * [Alternative Splicing Detection](#alternative-splicing-detection)
  * [Isoform Quantification](#isoform-quantification)
  * [Data Filtering](#data-filtering)
- [Functional Consequence Analysis](#functional-consequence-analysis)
  * [Step 1: Preliminary Analysis](#step-1-preliminary-analysis)
  * [Step 2: External Protein Analysis](#step-2-external-protein-analysis)
  * [Step 3: Final Integration](#step-3-final-integration)
- [Configuration Settings](#configuration-settings)
  * [Preprocessing and QC](#preprocessing-and-qc)
  * [UMI Handling](#umi-handling)
  * [Genome Alignment (STAR)](#genome-alignment-star)
  * [Splicing Analysis](#splicing-analysis)
  * [Transcriptome Quantification (RSEM)](#transcriptome-quantification-rsem)
- [Bibliography](#Bibliography)
  * [References](#references)
  * [Acknowledgments](#Acknowledgments)

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

### GTF Addition: 
When SplAdder is enabled, the pipeline offers a "GTF Addition" option. This generates an updated GTF file incorporating novel splicing events detected by SplAdder, which is then used for subsequent analysis steps.

## Isoform Quantification

Isoform-level quantification is handled by RSEM. This step requires alignment to the transcriptome rather than the genome. The pipeline can either accept pre-existing transcriptomic BAM files or generate them directly based on the provided settings.

## Data Filtering

Outputs from the detection tools undergo a filtering phase to retain only statistically significant splicing events.

## Functional Consequence Analysis

The final stage utilizes IsoformSwitchAnalyzeR to evaluate isoform switches and predict their functional impact on proteins. This analysis is divided into three distinct steps involving two different environments.

### Step 1: 
Preliminary Analysis Executed within the main pipeline/container. This step filters data, performs statistical testing, and reconstructs amino acid sequences.

### Step 2: 
External Protein Analysis This step requires the External Protein Analysis Bash script. Note: This script must be executed using the second Dockerfile provided (the sequence analysis container). 

It runs the following external tools to annotate the amino acid sequences generated in Step 1:

   CPAT (Coding Potential)

   Pfam (Protein Domains)

   SignalP (Signal Peptides)

   IUPred2A (Disordered Regions)
  
### Step 3: 
Final Integration Executed via the IsoformSwitchAnalyzer_Final.R script. 
Note: This script must be run using the primary Dockerfile. It requires the output directories from the External Protein Analysis (Step 2) as arguments to integrate the functional annotations and generate final visualizations.



# Configuration Settings

The pipeline behavior is controlled via the `settings.xlsx` file. Below are the available parameters.

### Preprocessing and QC
* **`trimming`**: Set to `TRUE` to enable adapter and quality trimming.
* **`fastqc`**: Set to `TRUE` to enable FastQC analysis.
* **`mode`**: Sequencing mode: `single` or `paired` (Required if alignment or trimming is enabled).
* **`adapter`**, **`adapter2`**: Adapter sequences to trim for Read 1 and Read 2, respectively (Required if `trimming = TRUE`).
* **`clip_5_1`**, **`clip_5_2`**: Number of bases to clip from the 5' end of Read 1 and Read 2.
* **`clip_3_1`**, **`clip_3_2`**: Number of bases to clip from the 3' end of Read 1 and Read 2.

### UMI Handling
* **`UMI`**: Set to `TRUE` to perform UMI extraction and deduplication.
* **`UMI_dedup`**: Set to `TRUE` to perform deduplication on existing BAM files containing UMI tags.
* **`umi_len`**: Length of the UMI sequence (Required if `UMI = TRUE`).
* **`umi_loc`**: Location of the UMI (e.g., `per_read`, `read1`, `read2`, `per_index`).
* **`umi_prefix`**: Separator character between the UMI and the read ID.
* **`umi_config`**: Path to the UMI_parallel configuration file (Required if `UMI_dedup = TRUE`).

### Genome Alignment (STAR)
* **`indexing`**: Set to `TRUE` to generate a genome index.
* **`genome_dir`**: Path to the existing indexed genome directory (Required if `indexing = FALSE`).
* **`alignment`**: Set to `TRUE` to perform genome alignment.
* **`gtf`**: Path to the GTF annotation file.
* **`fasta`**: Path to the genomic FASTA file.
* **`CPU`**: Number of threads allocated for STAR indexing and alignment.
* **`gene_info`**: Path to the `geneInfo.tab` file generated by STAR (Required if `indexing = FALSE`).

### Splicing Analysis
* **`Rmats`**: Set to `TRUE` to run rMATS.
* **`rmats_path`**: Directory containing condition separation text files (Required if `RSEM = TRUE` but `alignment = FALSE`).
* **`Spladder`**: Set to `TRUE` to run SplAdder.
* **`txt_path`**: Directory containing SplAdder input lists (`all_bams.txt`, `conditionA.txt`, `conditionB.txt`).
* **`bam_dir`**: Directory containing input BAM files (Required for rMATS/SplAdder if `alignment = FALSE`).
* **`GTF_addition`**: Set to `TRUE` to update the original GTF with novel splicing events detected by SplAdder.

### Transcriptome Quantification (RSEM)
* **`RSEM`**: Set to `TRUE` to perform isoform quantification.
* **`transcriptome_indexing`**: Set to `TRUE` to index the transcriptome (Required if `GTF_addition = TRUE`).
* **`transcriptome_genome_dir`**: Path to the existing indexed transcriptome (Required if `RSEM = TRUE` and `transcriptome_indexing = FALSE`).
* **`transcriptome_alignment`**: Set to `TRUE` to align reads to the transcriptome.
* **`bam_dir_transcriptome`**: Directory containing transcriptomic BAM files (Required if `transcriptome_alignment = FALSE`).
* **`gene_length`**: Path to the gene length file generated by `gtftools -l` (Required if `indexing = FALSE`).

# Bibliography

## References

#### Preprocessing & Quality Control

  * **FastQC**: Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. Available online at: [http://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * **fastp**: Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884-i890. [https://doi.org/10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
  * **Trim Galore**: Krueger, F. (2015). Trim Galore\!: A wrapper around Cutadapt and FastQC to consistently apply adapter and quality trimming to FastQ files. Available online at: [https://github.com/FelixKrueger/TrimGalore](https://github.com/FelixKrueger/TrimGalore)

#### Alignment & Deduplication

  * **STAR**: Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., ... & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1), 15-21. [https://doi.org/10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)
  * **SAMtools**: Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., ... & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. [https://doi.org/10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)
  * **UMI-tools**: Smith, T., Heger, A., & Sudbery, I. (2017). UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. [cite\_start]*Genome Research*, 27(3), 491–499[cite: 28]. [cite\_start][https://doi.org/10.1101/gr.209601.116](https://doi.org/10.1101/gr.209601.116) [cite: 29]
  * **UMI\_parallel**: A parallel wrapper for UMI-tools. Available online at: [https://github.com/Milda85/UMI\_parallel](https://github.com/Milda85/UMI_parallel)

#### Quantification

  * **featureCounts**: Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7), 923-930. [https://doi.org/10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656)
  * **RSEM**: Li, B., & Dewey, C. N. (2011). RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. *BMC Bioinformatics*, 12(1), 323. [https://doi.org/10.1186/1471-2105-12-323](https://doi.org/10.1186/1471-2105-12-323)

#### Alternative Splicing Analysis

  * **rMATS**: Shen, S., Park, J. W., Lu, Z. X., Lin, L., Henry, M. D., Wu, Y. N., ... & Xing, Y. (2014). rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. *Proceedings of the National Academy of Sciences*, 111(51), E5593-E5601. [https://doi.org/10.1073/pnas.1419161111](https://doi.org/10.1073/pnas.1419161111)
  * **SplAdder**: Kahles, A., Ongen, H., Zhong, Y., & Rätsch, G. (2016). SplAdder: identification, quantification and testing of alternative splicing events from RNA-Seq data. *Bioinformatics*, 32(12), 1840-1847. [https://doi.org/10.1093/bioinformatics/btw076](https://doi.org/10.1093/bioinformatics/btw076)

#### Functional Consequence Analysis

  * **IsoformSwitchAnalyzeR**: Vitting-Seerup, K., & Sandelin, A. (2019). IsoformSwitchAnalyzeR: Analysis of changes in genome-wide patterns of alternative splicing and its functional consequences. [cite\_start]*Bioinformatics*, 35(21), 4469-4471[cite: 1722, 2179]. [https://doi.org/10.1093/bioinformatics/btz247](https://doi.org/10.1093/bioinformatics/btz247)

## Acknowledgments

The concept for the script structure originated from a collaboration with **Sara Bianchi** from the **HEDGe Lab** (University of Turin).

Original repository: [https://github.com/sara-bianchi/Bulk_RNA_seq_pipeline](https://github.com/sara-bianchi/Bulk_RNA_seq_pipeline)

