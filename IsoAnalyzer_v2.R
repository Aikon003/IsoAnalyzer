#------------------------------------------------------------------------------
# Isoform Analyzer (v2)
#------------------------------------------------------------------------------
#
# table.xlsx columns (accessed by NAME, order does not matter):
#   sample_ID | cell_line | condition | replicate | fastq1 | fastq2 | bam | bam_transcriptome | cfr_1vs2
#   - fastq1 / fastq2 are needed only if trimming, alignment or RSEM transcriptome
#     alignment are run (single-end: only fastq1; an old "fastq" column is still accepted)
#   - bam: full path of the existing genomic BAM (needed if alignment = FALSE and
#     Rmats or Spladder are used; ignored if alignment = TRUE)
#   - bam_transcriptome: full path of the existing transcriptome BAM (Aligned.toTranscriptome),
#     needed if RSEM = TRUE and transcriptome_alignment = FALSE
#   - cfr_1vs2 assigns each sample to one of the two groups compared WITHIN each cell line
#     (the lower value is group 1 = condition A / b1, the higher is group 2 = B / b2)
#
# settings.xlsx (name in column 1, values from column 3 on). Main entries:
#   trimming, fastqc, mode, adapter, adapter2, clip_5_1, clip_3_1, clip_5_2, clip_3_2
#   UMI, umi_len, umi_loc, umi_prefix, umi_config, UMI_dedup
#   indexing, alignment, CPU, gtf, fasta, genome_dir, gene_length, gene_info
#   Rmats (TRUE/FALSE), read_length (optional, default 151)
#   Spladder: False | All | Test and Viz
#   spladder_graph_dir : output dir of a previous SplAdder run (only for "Test and Viz")
#   GTF_addition
#   RSEM (run quantification), transcriptome_indexing, transcriptome_alignment,
#   transcriptome_genome_dir
#   IsoformSwitch (run import + satuRn + Part1; default = value of RSEM)
#   rsem_counts_dir (existing RSEM counts, needed if IsoformSwitch = TRUE and RSEM = FALSE)
#   gtf_isoformswitch_exists (TRUE/FALSE), gtf_isoformswitch (path of an already
#   cleaned GTF, compatible with RSEM and IsoformSwitchAnalyzeR)
#------------------------------------------------------------------------------

#packages
library(tidyr)
library(stringr)
library(grid)
library(readxl)
library(dplyr)
library(readr)
library(ggplot2)
library(parallel)


#set up environment

set.seed(1234)
current_time <- format(Sys.time(), "%y%m%d")
base_dir <- paste0("./Output_", current_time)
dir.create(base_dir, showWarnings = FALSE)


#------------------------------------------------------------------------------
# helpers
#------------------------------------------------------------------------------

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))

check_status <- function(status, what) {
  if (!is.null(status) && !is.na(status) && status != 0) {
    stop(what, " failed (exit code ", status, ")")
  }
  invisible(status)
}

needs_zcat <- function(file) {
  if (is.null(file) || !is.character(file) || length(file) == 0 || is.na(file) || nchar(file) == 0) return(FALSE)
  endsWith(file, ".gz")
}


#upload settings

print("uploading settings file")

settings <- as.data.frame(read_excel("./settings.xlsx", col_names = FALSE))
settings <- settings[!is.na(settings$...1) & settings$...1 != "", , drop = FALSE]
rownames(settings) <- settings$...1
settings <- settings %>% dplyr::select(-c("...1", "...2"))
settings <- as.data.frame(t(settings), stringsAsFactors = FALSE)

# unique, non-empty, non-"NA" values of a setting (character(0) if absent)
setting_raw <- function(name) {
  if (!name %in% colnames(settings)) return(character(0))
  v <- unique(na.omit(trimws(as.character(settings[[name]]))))
  v[v != "" & toupper(v) != "NA"]
}

get_setting <- function(name, required = TRUE) {
  v <- setting_raw(name)
  if (length(v) > 1) { stop("invalid argument for ", name, ": more than one value found") }
  if (length(v) == 0) {
    if (required) stop("invalid argument for ", name, ": missing")
    return(NA_character_)
  }
  v
}

get_logical <- function(name, required = TRUE, default = FALSE) {
  v <- get_setting(name, required)
  if (is.na(v)) return(default)
  out <- as.logical(v)
  if (is.na(out)) stop("invalid argument for ", name, ": must be TRUE or FALSE (found '", v, "')")
  out
}

get_int <- function(name, required = TRUE, default = 0L) {
  v <- get_setting(name, required)
  if (is.na(v)) return(default)
  out <- suppressWarnings(as.integer(v))
  if (is.na(out)) stop("invalid argument for ", name, ": must be an integer (found '", v, "')")
  out
}

get_choice <- function(name, choices, required = TRUE) {
  v <- get_setting(name, required)
  if (is.na(v)) return(NA_character_)
  if (!v %in% choices) stop("invalid argument for ", name, ": must be one of ", paste(choices, collapse = ", "))
  v
}

get_path <- function(name, type = c("file", "dir"), required = TRUE) {
  type <- match.arg(type)
  v <- get_setting(name, required)
  if (is.na(v)) return(NA_character_)
  ok <- if (type == "file") file.exists(v) else dir.exists(v)
  if (!ok) {
    if (required) stop("invalid argument for ", name, ": ", type, " not found: ", v)
    # not needed for this run (e.g. placeholder path left in settings.xlsx)
    warning("setting ", name, " ignored: ", type, " not found: ", v)
    return(NA_character_)
  }
  v
}

get_bases <- function(name) {
  v <- get_setting(name, required = FALSE)
  if (is.na(v)) return(NA_character_)
  v <- toupper(v)
  bases <- strsplit(v, "")[[1]]
  if (any(!bases %in% c("A", "T", "G", "C"))) stop("invalid argument for ", name)
  v
}


#------------------------------------------------------------------------------
# main switches
#------------------------------------------------------------------------------

trimming  <- get_logical("trimming")
indexing  <- get_logical("indexing")
alignment <- get_logical("alignment")
UMI       <- get_logical("UMI")
RSEM      <- get_logical("RSEM")
Rmats     <- get_logical("Rmats")
GTF_addition <- get_logical("GTF_addition")
CPU       <- get_int("CPU")

params <- list(
  UMI_dedup = get_logical("UMI_dedup", required = FALSE, default = FALSE)
)

# SplAdder: False | All | Test and Viz
spladder_raw <- tolower(gsub("[ _-]+", " ", trimws(get_setting("Spladder"))))
Spladder_mode <- switch(spladder_raw,
                        "false"        = "none",
                        "all"          = "all",
                        "true"         = "all",
                        "test and viz" = "test_viz",
                        stop("invalid argument for Spladder: must be False, All or Test and Viz"))
Spladder_any   <- Spladder_mode != "none"
Spladder_build <- Spladder_mode == "all"

IsoformSwitch <- get_logical("IsoformSwitch", required = FALSE, default = RSEM)

if (GTF_addition && !Spladder_any) {
  stop("GTF_addition = TRUE needs SplAdder results: set Spladder to 'All' or 'Test and Viz'")
}

# read mode (single / paired)
need_mode <- trimming || alignment || RSEM || Rmats
mode <- if (need_mode) get_choice("mode", c("single", "paired")) else NA_character_


#------------------------------------------------------------------------------
# trimming / UMI settings
#------------------------------------------------------------------------------

fastqc <- FALSE
fqc <- ""
adapter <- NA_character_
adapter2 <- NA_character_
clip_5_1 <- clip_3_1 <- clip_5_2 <- clip_3_2 <- 0L

if (trimming) {
  fastqc <- get_logical("fastqc")
  fqc <- if (fastqc) " --fastqc" else ""
  adapter <- get_bases("adapter")
  clip_5_1 <- get_int("clip_5_1", required = FALSE, default = 0L)
  clip_3_1 <- get_int("clip_3_1", required = FALSE, default = 0L)
  if (mode == "paired") {
    adapter2 <- get_bases("adapter2")
    clip_5_2 <- get_int("clip_5_2", required = FALSE, default = 0L)
    clip_3_2 <- get_int("clip_3_2", required = FALSE, default = 0L)
  }
}

if (UMI) {
  umi_len    <- get_int("umi_len")
  umi_loc    <- get_choice("umi_loc", c("per_read", "read1", "read2", "per_index"))
  umi_prefix <- get_setting("umi_prefix")
  if (params$UMI_dedup) {
    umi_config <- get_path("umi_config", "file")
  }
}


#------------------------------------------------------------------------------
# reference / index / RSEM / IsoformSwitch settings
#------------------------------------------------------------------------------

transcriptome_indexing  <- FALSE
transcriptome_alignment <- FALSE

if (RSEM) {
  transcriptome_indexing  <- get_logical("transcriptome_indexing")
  transcriptome_alignment <- get_logical("transcriptome_alignment")
}

# RSEM reference directory (existing, or the one that will be created in this run)
transcriptome_genome_dir <- get_path("transcriptome_genome_dir", "dir",
                                     required = (RSEM && !transcriptome_indexing) || (IsoformSwitch && !RSEM))
if (RSEM && transcriptome_indexing) {
  rsem_ref_dir <- file.path(base_dir, "rsem_reference", "GRCh38_RSEM")
} else {
  rsem_ref_dir <- transcriptome_genome_dir
}

# where the RSEM counts are (produced now, or already existing)
rsem_counts_dir <- NA_character_
gtf_isoformswitch_exists <- FALSE
gtf_isoformswitch <- NA_character_

if (IsoformSwitch) {
  rsem_counts_dir <- if (RSEM) file.path(base_dir, "Transcriptome", "counts") else get_path("rsem_counts_dir", "dir")
  gtf_isoformswitch_exists <- get_logical("gtf_isoformswitch_exists", required = FALSE, default = FALSE)
  if (gtf_isoformswitch_exists) {
    gtf_isoformswitch <- get_path("gtf_isoformswitch", "file")
  }
}

need_gtf <- indexing || alignment || Rmats || Spladder_build || GTF_addition ||
  (RSEM && transcriptome_indexing) || (IsoformSwitch && !gtf_isoformswitch_exists)
gtf <- get_path("gtf", "file", required = need_gtf)

fasta <- get_path("fasta", "file", required = indexing || (RSEM && transcriptome_indexing))

genome_dir <- NA_character_
if (alignment && !indexing) {
  genome_dir <- get_path("genome_dir", "dir")
}

# optional files (not used downstream, loaded only if present)
gene_length <- get_path("gene_length", "file", required = FALSE)
geneInfo <- NULL
gene_info_path <- get_path("gene_info", "file", required = FALSE)
if (!is.na(gene_info_path)) {
  geneInfo <- read.table(gene_info_path, quote = "\"", comment.char = "", skip = 1)
}

# SplAdder graphs for already existing data
spladder_graph_dir <- NA_character_
if (Spladder_mode == "test_viz") {
  spladder_graph_dir <- get_path("spladder_graph_dir", "dir")
}

read_length <- get_int("read_length", required = FALSE, default = 151L)

# executables
need_star <- indexing || alignment || (RSEM && (transcriptome_indexing || transcriptome_alignment))
Star_bin <- Sys.which("STAR")
if (need_star && nchar(Star_bin) == 0) stop("STAR not found in PATH")

if (alignment == TRUE) print("alignment = TRUE")


#------------------------------------------------------------------------------
# table
#------------------------------------------------------------------------------

print("uploading table file")

table <- read_excel("./table.xlsx")

if (!"sample_ID" %in% colnames(table)) stop("sample_ID column is missing in table.xlsx")
if (anyNA(table$sample_ID) || any(table$sample_ID == "")) stop("sample_ID contains missing or empty values.")
if (anyDuplicated(table$sample_ID)) stop("Duplicate sample_ID values found in table.xlsx.")

# fastq files
needs_fastq <- trimming || alignment || (RSEM && transcriptome_alignment)
if (needs_fastq) {
  # backward compatibility: single-end tables with a "fastq" column
  if (mode == "single" && !"fastq1" %in% colnames(table) && "fastq" %in% colnames(table)) {
    table <- dplyr::rename(table, fastq1 = fastq)
  }
  fq_cols <- if (mode == "single") "fastq1" else c("fastq1", "fastq2")
  missing_fq_cols <- setdiff(fq_cols, colnames(table))
  if (length(missing_fq_cols) > 0) {
    stop("Missing column(s) in table.xlsx: ", paste(missing_fq_cols, collapse = ", "))
  }
  for (col in fq_cols) {
    bad <- table[[col]][is.na(table[[col]]) | !file.exists(as.character(table[[col]]))]
    if (length(bad) > 0) stop("incorrect path in ", col, " column: ", paste(bad, collapse = ", "))
  }
}

# existing BAM files: full paths given directly in the table (no automatic search)
check_path_column <- function(col, why) {
  if (!col %in% colnames(table)) stop("Missing column in table.xlsx: ", col, " (", why, ")")
  paths <- as.character(table[[col]])
  bad <- table$sample_ID[is.na(paths) | paths == "" | !file.exists(paths)]
  if (length(bad) > 0) {
    stop("incorrect or missing path in column ", col, " for samples: ", paste(bad, collapse = ", "))
  }
  if (anyDuplicated(paths)) {
    stop("column ", col, " contains the same file for more than one sample: ",
         paste(unique(paths[duplicated(paths)]), collapse = ", "))
  }
}
if (!alignment && (Rmats || Spladder_any)) {
  check_path_column("bam", "needed for Rmats/Spladder when alignment = FALSE")
}
if (RSEM && !transcriptome_alignment) {
  check_path_column("bam_transcriptome", "needed for RSEM when transcriptome_alignment = FALSE")
}

# design / comparisons (used by rMATS, SplAdder and IsoformSwitchAnalyzeR)
cell_line_specs <- NULL
myDesign <- NULL
myComparisons <- NULL

if (Rmats || Spladder_any || IsoformSwitch) {

  required_columns <- c("sample_ID", "cell_line", "condition", "replicate", "cfr_1vs2")
  missing_columns <- setdiff(required_columns, colnames(table))
  if (length(missing_columns) > 0) {
    stop("Missing required columns in table.xlsx: ", paste(missing_columns, collapse = ", "))
  }
  for (col in c("cell_line", "condition", "replicate", "cfr_1vs2")) {
    if (anyNA(table[[col]]) || any(table[[col]] == "")) stop(col, " contains missing or empty values.")
  }

  table$cell_line <- as.character(table$cell_line)
  table$condition <- as.character(table$condition)
  table$cfr_1vs2  <- as.character(table$cfr_1vs2)
  table$combined_condition <- paste(table$cell_line, table$condition, sep = "_")

  comparison_map <- table %>%
    dplyr::distinct(cell_line, cfr_1vs2, condition, combined_condition) %>%
    dplyr::arrange(match(cell_line, unique(table$cell_line)), cfr_1vs2)

  # each cfr_1vs2 group must map to exactly one condition within each cell line
  mapping_check <- comparison_map %>%
    dplyr::group_by(cell_line, cfr_1vs2) %>%
    dplyr::summarise(n_conditions = dplyr::n_distinct(condition), .groups = "drop")
  if (any(mapping_check$n_conditions != 1)) {
    stop("Each cfr_1vs2 group must map to exactly one condition within each cell line.")
  }

  # exactly two groups per cell line
  group_check <- comparison_map %>%
    dplyr::distinct(cell_line, cfr_1vs2) %>%
    dplyr::count(cell_line, name = "n_groups")
  if (any(group_check$n_groups != 2)) {
    stop("Each cell line must have exactly two cfr_1vs2 groups. Found: ",
         paste(group_check$cell_line, group_check$n_groups, sep = "=", collapse = "; "))
  }

  # one spec per cell line: group 1 (lower cfr_1vs2) = A / b1, group 2 = B / b2
  cell_line_specs <- lapply(unique(table$cell_line), function(cl) {
    g <- comparison_map[comparison_map$cell_line == cl, , drop = FALSE]
    list(
      cell_line = cl,
      cond1 = g$condition[1], cond2 = g$condition[2],
      combined1 = g$combined_condition[1], combined2 = g$combined_condition[2],
      samples1 = table$sample_ID[table$cell_line == cl & table$cfr_1vs2 == g$cfr_1vs2[1]],
      samples2 = table$sample_ID[table$cell_line == cl & table$cfr_1vs2 == g$cfr_1vs2[2]]
    )
  })
  names(cell_line_specs) <- unique(table$cell_line)

  myComparisons <- dplyr::bind_rows(lapply(cell_line_specs, function(s) {
    data.frame(condition_1 = s$combined1, condition_2 = s$combined2, stringsAsFactors = FALSE)
  }))

  myDesign <- data.frame(
    sampleID = as.character(table$sample_ID),
    condition = as.character(table$combined_condition),
    stringsAsFactors = FALSE
  )
  if (length(unique(table$replicate)) > 1) {
    myDesign$replicate <- as.character(table$replicate)
  }

  cat("\nComparisons (condition_1 vs condition_2):\n")
  print(myComparisons)
}


#------------------------------------------------------------------------------
# trimming
#------------------------------------------------------------------------------

trim_dir <- file.path(base_dir, "trimmed")

if (trimming == TRUE) {
  print("starting trimming")
  dir.create(trim_dir, showWarnings = FALSE)

  trim_galore_bin <- if (nchar(Sys.which("trim_galore")) > 0) Sys.which("trim_galore") else "/usr/bin/trim_galore"
  no_adapter <- (is.na(adapter) || adapter == "")
  no_adapter2 <- (is.na(adapter2) || adapter2 == "")

  for (i in seq_len(nrow(table))) {
    sample_ID <- table$sample_ID[i]
    html <- file.path(trim_dir, paste0(sample_ID, "_fastp_report.html"))
    json <- file.path(trim_dir, paste0(sample_ID, "_fastp_report.json"))

    if (mode == "single") {

      fastq <- table$fastq1[i]

      if (UMI == TRUE) {
        # -- Single End + UMI (Fastp) --
        cleaned_fastq <- file.path(trim_dir, paste0(sample_ID, "_umi.fq.gz"))
        st <- system2("fastp", paste0(
          "-i ", fastq, " ",
          "-o ", cleaned_fastq, " ",
          "--umi ",
          "--umi_loc ", umi_loc, " ",
          "--umi_len ", umi_len, " ",
          "--dedup ",
          "--umi_prefix ", umi_prefix, " ",
          ifelse(no_adapter, "--detect_adapter_for_pe", paste0("--adapter_sequence ", adapter)),
          " --thread ", CPU,
          " --html ", html, " --json ", json
        ))
        check_status(st, paste("fastp (single, UMI) for", sample_ID))
        if (fastqc) check_status(system2("fastqc", paste0(cleaned_fastq, " -o ", trim_dir, " -t ", CPU)), "fastqc")
        fastq <- cleaned_fastq

      } else if (no_adapter) {
        # -- Single End + Auto Adapter (Fastp) --
        cleaned_fastq <- file.path(trim_dir, paste0(sample_ID, "_trimmed.fq.gz"))
        st <- system2("fastp", paste0(
          "-i ", fastq, " ",
          "-o ", cleaned_fastq, " ",
          "--detect_adapter_for_pe ",
          ifelse(clip_5_1 == 0, "", paste0("--trim_front1 ", clip_5_1, " ")),
          ifelse(clip_3_1 == 0, "", paste0("--trim_tail1 ", clip_3_1, " ")),
          "--thread ", CPU,
          " --html ", html, " --json ", json
        ))
        check_status(st, paste("fastp (single) for", sample_ID))
        if (fastqc) check_status(system2("fastqc", paste0(cleaned_fastq, " -o ", trim_dir, " -t ", CPU)), "fastqc")
        fastq <- cleaned_fastq

      } else {
        # -- Single End + Adapter (Trim Galore) --
        clipping_R1_3 <- if (clip_3_1 == 0) "" else paste0(" --three_prime_clip_R1 ", clip_3_1)
        clipping_R1_5 <- if (clip_5_1 == 0) "" else paste0(" --clip_R1 ", clip_5_1)
        st <- system2(trim_galore_bin, paste0(
          "--adapter ", adapter,
          clipping_R1_5, clipping_R1_3,
          " --gzip", fqc,
          " --basename ", sample_ID,
          " --output_dir ", trim_dir, " ",
          fastq
        ))
        check_status(st, paste("trim_galore (single) for", sample_ID))
        fastq <- file.path(trim_dir, paste0(sample_ID, "_trimmed.fq.gz"))
      }
      table$fastq1[i] <- fastq

    } else if (mode == "paired") {

      fastq1 <- table$fastq1[i]
      fastq2 <- table$fastq2[i]

      if (UMI == TRUE) {
        # -- Paired End + UMI (Fastp) --
        cleaned_fastq1 <- file.path(trim_dir, paste0(sample_ID, "_umi_R1.fq.gz"))
        cleaned_fastq2 <- file.path(trim_dir, paste0(sample_ID, "_umi_R2.fq.gz"))
        st <- system2("fastp", paste0(
          "-i ", fastq1, " ",
          "-I ", fastq2, " ",
          "-o ", cleaned_fastq1, " ",
          "-O ", cleaned_fastq2, " ",
          "--umi ",
          "--dedup ",
          "--umi_loc ", umi_loc, " ",
          "--umi_len ", umi_len, " ",
          "--umi_prefix ", umi_prefix, " ",
          ifelse(no_adapter && no_adapter2,
                 "--detect_adapter_for_pe",
                 paste0("--adapter_sequence ", adapter, " --adapter_sequence_r2 ", adapter2)),
          " --thread ", CPU,
          " --html ", html, " --json ", json
        ))
        check_status(st, paste("fastp (paired, UMI) for", sample_ID))
        if (fastqc) check_status(system2("fastqc", paste0(cleaned_fastq1, " ", cleaned_fastq2, " -o ", trim_dir, " -t ", CPU)), "fastqc")
        fastq1 <- cleaned_fastq1
        fastq2 <- cleaned_fastq2

      } else if (no_adapter && no_adapter2) {
        # -- Paired End + Auto Adapter (Fastp) --
        cleaned_fastq1 <- file.path(trim_dir, paste0(sample_ID, "_val_1.fq.gz"))
        cleaned_fastq2 <- file.path(trim_dir, paste0(sample_ID, "_val_2.fq.gz"))
        st <- system2("fastp", paste0(
          "-i ", fastq1, " ",
          "-I ", fastq2, " ",
          "-o ", cleaned_fastq1, " ",
          "-O ", cleaned_fastq2, " ",
          "--detect_adapter_for_pe ",
          ifelse(clip_5_1 == 0, "", paste0("--trim_front1 ", clip_5_1, " ")),
          ifelse(clip_3_1 == 0, "", paste0("--trim_tail1 ", clip_3_1, " ")),
          ifelse(clip_5_2 == 0, "", paste0("--trim_front2 ", clip_5_2, " ")),
          ifelse(clip_3_2 == 0, "", paste0("--trim_tail2 ", clip_3_2, " ")),
          "--thread ", CPU,
          " --html ", html, " --json ", json
        ))
        check_status(st, paste("fastp (paired) for", sample_ID))
        if (fastqc) check_status(system2("fastqc", paste0(cleaned_fastq1, " ", cleaned_fastq2, " -o ", trim_dir, " -t ", CPU)), "fastqc")
        fastq1 <- cleaned_fastq1
        fastq2 <- cleaned_fastq2

      } else {
        # -- Paired End + Adapter (Trim Galore) --
        clipping_R1_3 <- if (clip_3_1 == 0) "" else paste0(" --three_prime_clip_R1 ", clip_3_1)
        clipping_R1_5 <- if (clip_5_1 == 0) "" else paste0(" --clip_R1 ", clip_5_1)
        clipping_R2_3 <- if (clip_3_2 == 0) "" else paste0(" --three_prime_clip_R2 ", clip_3_2)
        clipping_R2_5 <- if (clip_5_2 == 0) "" else paste0(" --clip_R2 ", clip_5_2)
        st <- system2(trim_galore_bin, paste0(
          "--adapter ", adapter,
          " --adapter2 ", adapter2,
          clipping_R1_5, clipping_R1_3,
          clipping_R2_5, clipping_R2_3,
          " --gzip --paired", fqc,
          " --basename ", sample_ID,
          " --output_dir ", trim_dir, " ",
          fastq1, " ", fastq2
        ))
        check_status(st, paste("trim_galore (paired) for", sample_ID))
        fastq1 <- file.path(trim_dir, paste0(sample_ID, "_val_1.fq.gz"))
        fastq2 <- file.path(trim_dir, paste0(sample_ID, "_val_2.fq.gz"))
      }

      table$fastq1[i] <- fastq1
      table$fastq2[i] <- fastq2
    }
  }
  print("finishing trimming")
}


#------------------------------------------------------------------------------
# indexing
#------------------------------------------------------------------------------

if (indexing == TRUE) {
  print("starting indexing")

  system2("mkdir", "/home/genome")
  system2("chmod", "777 /home/genome")
  st <- system2(Star_bin, paste0("--runThreadN ", CPU, " --runMode genomeGenerate --limitGenomeGenerateRAM=150000000000 --genomeDir /home/genome --genomeFastaFiles ", fasta, " --sjdbGTFfile ", gtf))
  check_status(st, "STAR genomeGenerate")
  system2("rm", paste0("-r ", base_dir, "/genome"))
  system2("mv", paste0("/home/genome ", base_dir))
  genome_dir <- file.path(base_dir, "genome")
  system2("gtftools", paste0("-l ", base_dir, "/gene_length.txt ", gtf))
  gene_length <- file.path(base_dir, "gene_length.txt")
  geneInfo <- read.table(file.path(genome_dir, "geneInfo.tab"), quote = "\"", comment.char = "", skip = 1)

  print("finishing indexing")
}


#------------------------------------------------------------------------------
# alignment
#------------------------------------------------------------------------------

# fastq files to use for a sample (trimmed if trimming was run, otherwise raw)
get_fastq_files <- function(sample_ID, i) {
  if (mode == "single") {
    if (trimming) {
      suffix <- if (UMI) "_umi.fq.gz" else "_trimmed.fq.gz"
      return(file.path(trim_dir, paste0(sample_ID, suffix)))
    }
    return(table$fastq1[i])
  } else {
    if (trimming) {
      suffix1 <- if (UMI) "_umi_R1.fq.gz" else "_val_1.fq.gz"
      suffix2 <- if (UMI) "_umi_R2.fq.gz" else "_val_2.fq.gz"
      return(list(file.path(trim_dir, paste0(sample_ID, suffix1)),
                  file.path(trim_dir, paste0(sample_ID, suffix2))))
    }
    return(list(table$fastq1[i], table$fastq2[i]))
  }
}

if (alignment == TRUE) {
  print("starting alignment stage")

  bam_parent_dir <- file.path(base_dir, "BAM")
  counts_parent_dir <- file.path(base_dir, "counts")
  logs_parent_dir <- file.path(base_dir, "logs")
  dir.create(bam_parent_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(counts_parent_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(logs_parent_dir, showWarnings = FALSE, recursive = TRUE)

  # --- Loop 1: STAR Alignment and Indexing ---
  print("--- Starting Loop 1: STAR Alignment & Indexing ---")
  table$bam_aligned_genomic <- NA_character_

  for (i in seq_len(nrow(table))) {
    sample_ID <- table$sample_ID[i]
    cat("  Processing sample:", sample_ID, "(", i, "/", nrow(table), ")\n")

    sample_bam_dir <- file.path(bam_parent_dir, sample_ID)
    dir.create(sample_bam_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(logs_parent_dir, sample_ID), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(counts_parent_dir, sample_ID), showWarnings = FALSE, recursive = TRUE)

    star_prefix <- file.path(sample_bam_dir, sample_ID)
    fastq_files <- get_fastq_files(sample_ID, i)

    cat("    Running STAR alignment...\n")
    if (mode == "single") {
      reads_in <- shQuote(fastq_files)
      zcat_flag <- if (needs_zcat(fastq_files)) " --readFilesCommand zcat" else ""
    } else {
      reads_in <- paste(shQuote(fastq_files[[1]]), shQuote(fastq_files[[2]]))
      zcat_flag <- if (needs_zcat(fastq_files[[1]])) " --readFilesCommand zcat" else ""
    }

    status_star <- system2(Star_bin, paste0(
      "--runThreadN ", CPU,
      " --genomeDir ", shQuote(genome_dir),
      " --readFilesIn ", reads_in,
      " --outSAMtype BAM SortedByCoordinate",
      zcat_flag,
      " --outFileNamePrefix ", shQuote(star_prefix)
    ))
    if (status_star != 0) stop("STAR alignment failed for sample: ", sample_ID)

    bam_aligned_path <- paste0(star_prefix, "Aligned.sortedByCoord.out.bam")
    if (!file.exists(bam_aligned_path)) stop("Aligned BAM not found after STAR run for sample: ", sample_ID)

    cat("    Running samtools index...\n")
    status_index <- system2("samtools", paste("index", shQuote(bam_aligned_path)))
    if (status_index != 0) stop("Samtools index failed for sample: ", sample_ID)
    print(paste("  Finished STAR and indexing for", sample_ID))

    table$bam_aligned_genomic[i] <- bam_aligned_path
  }
  print("--- Finished Loop 1: STAR Alignment & Indexing ---")


  # --- UMI deduplication (optional) ---
  table$bam_final_genomic <- NA_character_

  if (UMI == TRUE && params$UMI_dedup == TRUE) {
    cat("\n--- Starting Parallel UMI Deduplication ---\n")

    dedup_input_dir <- bam_parent_dir
    dedup_output_dir <- file.path(base_dir, "BAM_dedup")
    dir.create(dedup_output_dir, showWarnings = FALSE, recursive = TRUE)

    umi_script <- "/opt/UMI_parallel/para_umi_dedup.sh"
    if (!file.exists(umi_script)) stop("UMI parallel script not found at expected path: ", umi_script)

    system(paste("chmod +x", shQuote(umi_script)))

    cmd <- paste(
      "bash", shQuote(umi_script),
      "-i", shQuote(dedup_input_dir),
      "-o", shQuote(dedup_output_dir),
      "-t", CPU,
      "-f", shQuote(umi_config)
    )
    cat("Running parallel deduplication command:", cmd, "\n")
    status_dedup <- system(cmd)
    if (status_dedup != 0) stop("Error executing para_umi_dedup.sh. Check logs or script output.")

    cat("Parallel UMI deduplication finished.\n")
    for (i in seq_len(nrow(table))) {
      sample_ID <- table$sample_ID[i]
      dedup_bam_path <- file.path(dedup_output_dir, paste0(sample_ID, "Aligned.sortedByCoord.out.dedup.bam"))
      if (!file.exists(dedup_bam_path)) {
        warning("Expected deduplicated BAM not found for sample ", sample_ID, ": ", dedup_bam_path)
      } else {
        table$bam_final_genomic[i] <- dedup_bam_path
      }
    }

    cat("Indexing deduplicated BAM files...\n")
    for (bam_path in na.omit(table$bam_final_genomic)) {
      cat("  Indexing:", basename(bam_path), "\n")
      status_index_dedup <- system2("samtools", paste("index", shQuote(bam_path)))
      if (status_index_dedup != 0) warning("Samtools index failed for deduplicated BAM: ", bam_path)
    }
    cat("Finished indexing deduplicated BAM files.\n")

  } else {
    print("UMI deduplication skipped (UMI is FALSE or UMI_dedup is FALSE).")
    table$bam_final_genomic <- table$bam_aligned_genomic
  }


  # --- Loop 2: featureCounts ---
  print("--- Starting Loop 2: featureCounts ---")
  for (i in seq_len(nrow(table))) {
    sample_ID <- table$sample_ID[i]
    bam_file_for_counting <- table$bam_final_genomic[i]
    if (is.na(bam_file_for_counting) || !file.exists(bam_file_for_counting)) {
      warning("Skipping featureCounts for ", sample_ID, ": final BAM path not found or file missing.")
      next
    }
    cat("  Counting features for sample:", sample_ID, "(", i, "/", nrow(table), ")\n")
    fc_out_file <- file.path(counts_parent_dir, sample_ID, paste0(sample_ID, "_counts.txt"))
    fc_log_file <- file.path(logs_parent_dir, sample_ID, paste0(sample_ID, "_featureCounts.log"))
    fc_cmd <- paste0(
      "featureCounts",
      " -a ", shQuote(gtf),
      " -o ", shQuote(fc_out_file),
      " -T ", CPU,
      ifelse(mode == "paired", " -p", ""),
      " ", shQuote(bam_file_for_counting),
      " > ", shQuote(fc_log_file), " 2>&1"
    )
    status_fc <- system(fc_cmd)
    if (status_fc != 0) {
      warning("featureCounts failed for sample: ", sample_ID, ". Check log: ", fc_log_file)
    } else {
      print(paste("  Finished featureCounts for", sample_ID))
    }
  }

  print("finishing alignment stage")
}


#------------------------------------------------------------------------------
# BAM paths for rMATS / SplAdder (fresh alignment, or "bam" column of the table)
#------------------------------------------------------------------------------

if (Rmats || Spladder_any) {
  table$bam_path <- if (alignment == TRUE) table$bam_final_genomic else as.character(table$bam)
  missing_bam <- is.na(table$bam_path) | !file.exists(table$bam_path)
  if (any(missing_bam)) {
    stop("BAM files not found for samples: ", paste(table$sample_ID[missing_bam], collapse = ", "))
  }
  table$bam_path <- normalizePath(table$bam_path)
}


#------------------------------------------------------------------------------
# rMATS (one run per cell line: group 1 vs group 2 of cfr_1vs2)
#------------------------------------------------------------------------------

if (Rmats == TRUE) {
  print("Preparing and running rMATS (one run per cell line)")

  rmats_base <- file.path(base_dir, "rmatsout")

  for (spec in cell_line_specs) {
    cl_safe <- safe_name(spec$cell_line)
    out_cl  <- file.path(rmats_base, cl_safe)
    out_OD  <- file.path(out_cl, "OD")
    out_TM  <- file.path(out_cl, "TM")
    out_txt <- file.path(out_cl, "txt")
    for (d in c(out_OD, out_TM, out_txt)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

    b1 <- table$bam_path[match(spec$samples1, table$sample_ID)]
    b2 <- table$bam_path[match(spec$samples2, table$sample_ID)]
    message("[", spec$cell_line, "] n B1: ", length(b1), " | n B2: ", length(b2))

    b1_path <- file.path(out_txt, "B1.txt")
    b2_path <- file.path(out_txt, "B2.txt")
    write(paste(b1, collapse = ","), file = b1_path)
    write(paste(b2, collapse = ","), file = b2_path)

    st <- system2("python3", args = c(
      "/rmats/rmats.py",
      "--b1", b1_path,
      "--b2", b2_path,
      "--gtf", gtf,
      "-t", mode,
      "--readLength", as.character(read_length),
      "--variable-read-length",
      "--nthread", as.character(CPU),
      "--od", out_OD,
      "--tmp", out_TM
    ))
    check_status(st, paste("rMATS for", spec$cell_line))
  }
  print("Finished rMATS")
}


#------------------------------------------------------------------------------
# SplAdder: False | All (build + test + viz) | Test and Viz (existing graphs)
#------------------------------------------------------------------------------

spladder_out <- NA_character_

if (Spladder_any) {

  spladder_bin <- Sys.which("spladder")
  if (nchar(spladder_bin) == 0) stop("spladder not found in PATH")

  mc.cores <- max(1, min(4, ceiling(detectCores() / 4)))

  spladder_out <- if (Spladder_build) file.path(base_dir, "spladder_out") else normalizePath(spladder_graph_dir)
  dir.create(spladder_out, recursive = TRUE, showWarnings = FALSE)
  spladder_dir <- file.path(base_dir, "spladder")
  spladder_txt <- file.path(spladder_dir, "txt")
  spladder_logs <- file.path(spladder_out, "logs")
  dir.create(spladder_txt, recursive = TRUE, showWarnings = FALSE)
  dir.create(spladder_logs, recursive = TRUE, showWarnings = FALSE)

  if (Spladder_mode == "test_viz") {
    message("SplAdder: using existing graphs in ", spladder_out,
            " (BAM file names must match those used when the graphs were built)")
  }

  # --- condition files, one pair per cell line ---
  spladder_specs <- lapply(cell_line_specs, function(spec) {
    cl_safe <- safe_name(spec$cell_line)
    labelA <- safe_name(paste0(spec$cell_line, "_", spec$cond1))
    labelB <- safe_name(paste0(spec$cell_line, "_", spec$cond2))
    fileA <- file.path(spladder_txt, paste0(cl_safe, "_conditionA.txt"))
    fileB <- file.path(spladder_txt, paste0(cl_safe, "_conditionB.txt"))
    writeLines(table$bam_path[match(spec$samples1, table$sample_ID)], fileA)
    writeLines(table$bam_path[match(spec$samples2, table$sample_ID)], fileB)
    list(cell_line = spec$cell_line, cl_safe = cl_safe,
         labelA = labelA, labelB = labelB, fileA = fileA, fileB = fileB,
         testcase = paste0("testing_", labelA, "_vs_", labelB))
  })

  # ------------------------------------------------------------------
  # Build (only for "All")
  # ------------------------------------------------------------------
  if (Spladder_build) {
    print("Starting SplAdder build")

    all_bams <- unique(table$bam_path)
    all_bams_file <- file.path(spladder_txt, "all_bams.txt")
    writeLines(all_bams, all_bams_file)

    # filtered GTF (protein_coding only)
    gtf_pc <- file.path(spladder_dir, "Homo_sapiens.GRCh38.113.prot.gtf")
    if (!file.exists(gtf_pc)) {
      message("Creating filtered GTF (protein_coding only)...")
      system(paste("grep 'protein_coding' ", shQuote(gtf), " > ", shQuote(gtf_pc), sep = ""))
      if (!file.exists(gtf_pc) || file.info(gtf_pc)$size == 0) {
        stop("Filtered GTF file was not created correctly. Check path or grep command.")
      }
    } else {
      message("Filtered GTF already exists, skipping creation.")
    }
    message("Filtered GTF created at: ", gtf_pc)

    log_step2  <- file.path(spladder_logs, "spladder_step2_merge_graphs.log")
    log_step3b <- file.path(spladder_logs, "spladder_step3b_quant_collect.log")
    log_step4  <- file.path(spladder_logs, "spladder_step4_event_extraction.log")

    run_parallel <- function(items, fun) {
      res <- mclapply(items, fun, mc.cores = mc.cores)
      failed <- vapply(res, function(r) inherits(r, "try-error") || !identical(r, 0L), logical(1))
      if (any(failed)) warning("SplAdder step failed for: ", paste(basename(items)[failed], collapse = ", "))
      invisible(res)
    }

    # 1) single graphs
    build_single_graph <- function(bam) {
      args <- c("build", "--outdir", spladder_out, "--annotation", gtf_pc,
                "--bams", bam, "--merge-strat", "single",
                "--no-extract-ase", "--set-mm-tag", "nM")
      message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Building single graph for: ", bam)
      lg <- file.path(spladder_logs, paste0("spladder_step1_", basename(bam), ".log"))
      as.integer(system2(spladder_bin, args = args, stdout = lg, stderr = lg))
    }
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Step 1: Building single graphs in parallel")
    run_parallel(all_bams, build_single_graph)

    # 2) merge graphs
    merge_args <- c("build", "--outdir", spladder_out, "--annotation", gtf_pc,
                    "--bams", all_bams_file, "--merge-strat", "merge_graphs",
                    "--no-extract-ase", "--set-mm-tag", "nM")
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Step 2: Merging graphs into a single merged graph")
    check_status(system2(spladder_bin, args = merge_args, stdout = log_step2, stderr = log_step2), "SplAdder merge graphs")

    # 3a) quantification per BAM
    quant_single <- function(bam) {
      args <- c("build", "--outdir", spladder_out, "--annotation", gtf_pc,
                "--bams", bam, "--merge-strat", "merge_graphs",
                "--no-extract-ase", "--quantify-graph", "--qmode", "single",
                "--set-mm-tag", "nM")
      message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Quantifying merged graph for: ", bam)
      lg <- file.path(spladder_logs, paste0("spladder_step3a_", basename(bam), ".log"))
      as.integer(system2(spladder_bin, args = args, stdout = lg, stderr = lg))
    }
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Step 3a: Quantification per BAM in parallel (qmode single)")
    run_parallel(all_bams, quant_single)

    # 3b) quantification collection
    collect_args <- c("build", "--outdir", spladder_out, "--annotation", gtf_pc,
                      "--bams", all_bams_file, "--merge-strat", "merge_graphs",
                      "--no-extract-ase", "--quantify-graph", "--qmode", "collect",
                      "--set-mm-tag", "nM")
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Step 3b: Collecting all quantifications (qmode collect)")
    check_status(system2(spladder_bin, args = collect_args, stdout = log_step3b, stderr = log_step3b), "SplAdder quantification collect")

    # 4) event extraction
    event_args <- c("build", "--outdir", spladder_out, "--annotation", gtf_pc,
                    "--bams", all_bams_file, "--merge-strat", "merge_graphs",
                    "--event-types", "exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons",
                    "--ase-edge-limit", "1000", "--set-mm-tag", "nM")
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Step 4: Extracting ASE events on the merged graph")
    check_status(system2(spladder_bin, args = event_args, stdout = log_step4, stderr = log_step4), "SplAdder event extraction")

    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            " | Build completed: single graphs -> merge -> quantification -> event extraction")
  }

  # ------------------------------------------------------------------
  # Test (All and Test and Viz): one test per cell line
  # ------------------------------------------------------------------
  print("Starting SplAdder test")
  for (sp in spladder_specs) {
    test_args <- c(
      "test",
      "--conditionA", sp$fileA,
      "--conditionB", sp$fileB,
      "--labelA", sp$labelA,
      "--labelB", sp$labelB,
      "--parallel", as.character(CPU),
      "--diagnose-plots",
      "--plot-format", "pdf",
      "--outdir", spladder_out
    )
    message("Running SplAdder test for ", sp$cell_line, " (", sp$labelA, " vs ", sp$labelB, ")")
    lg <- file.path(spladder_logs, paste0("spladder_test_", sp$cl_safe, ".log"))
    st <- system2(spladder_bin, args = test_args, stdout = lg, stderr = lg)
    if (st != 0) warning("SplAdder test failed for ", sp$cell_line, ". Check log: ", lg)
  }

  # ------------------------------------------------------------------
  # Viz (All and Test and Viz): top events per cell line
  # ------------------------------------------------------------------
  message("Starting SplAdder viz (parallel mode)")

  event_types <- c("exon_skip", "intron_retention", "alt_3prime",
                   "alt_5prime", "mutex_exons", "mult_exon_skip")
  n_jobs <- max(1, min(length(event_types), CPU))

  for (sp in spladder_specs) {
    test_dir <- file.path(spladder_out, sp$testcase)
    if (!dir.exists(test_dir)) {
      warning("SplAdder test directory not found for ", sp$cell_line, ": ", test_dir, ". Skipping viz.")
      next
    }

    files_A <- readLines(sp$fileA); files_A <- files_A[files_A != ""]
    files_B <- readLines(sp$fileB); files_B <- files_B[files_B != ""]
    stringa_A <- paste0("ConditionA:", paste(files_A, collapse = ","))
    stringa_B <- paste0("ConditionB:", paste(files_B, collapse = ","))

    plot_top_events <- function(event_type, top_n = 10, log_scale = FALSE) {
      args <- c("viz",
                "--test", sp$testcase, event_type, as.character(top_n),
                "--testdir", test_dir,
                "--outdir", spladder_out,
                "--confidence", "3",
                if (log_scale) "--log",
                "--track", "coverage", stringa_A, stringa_B,
                "--track", "splicegraph",
                "--outbase", paste0(sp$cl_safe, "_top", top_n, if (log_scale) "_log_scale" else "", "_"),
                "--format", "pdf")
      message("Plotting ", if (log_scale) "logarithmic " else "", "top ", top_n, " events for ", sp$cell_line, ": ", event_type)
      lg <- file.path(spladder_logs, paste0("spladder_viz_", sp$cl_safe, "_", event_type, if (log_scale) "_log" else "", ".log"))
      as.integer(system2(spladder_bin, args = args, stdout = lg, stderr = lg))
    }

    mclapply(event_types, plot_top_events, mc.cores = n_jobs)
    mclapply(event_types, plot_top_events, log_scale = TRUE, mc.cores = n_jobs)

    message("All top events plotted for ", sp$cell_line, " in ", test_dir)
  }
}


#------------------------------------------------------------------------------
# GTF addition (SplAdder events added to the reference GTF)
#------------------------------------------------------------------------------

gtf_merged <- NA_character_
spladder_combined_gff <- NA_character_

if (GTF_addition == TRUE) {

  if (!file.exists(gtf)) stop("Original GTF file not found.")
  if (!dir.exists(spladder_out)) stop("SplAdder output directory not found.")

  # GFF3 produced by SplAdder for each category
  gff3_files <- list.files(
    path = spladder_out,
    pattern = "merge_graphs_.*_C[0-9]+\\.confirmed\\.gff3$",
    full.names = TRUE
  )
  if (length(gff3_files) == 0) stop("No merge_graphs_*.confirmed.gff3 files found in ", spladder_out)

  # GFF3 merge
  spladder_combined_gff <- file.path(spladder_out, "spladder_combined.gff3")
  system(paste("cat", paste(shQuote(gff3_files), collapse = " "), ">", shQuote(spladder_combined_gff)))

  # GFF3 to GTF
  gff3_to_gtf <- Sys.which("rsem-gff3-to-gtf")
  if (nchar(gff3_to_gtf) == 0) stop("rsem-gff3-to-gtf not found in PATH")

  spladder_gtf <- file.path(spladder_out, "spladder_combined.gtf")
  check_status(system2("python3", args = c(gff3_to_gtf, spladder_combined_gff, spladder_gtf, "--RNA-patterns", "mRNA")),
               "rsem-gff3-to-gtf")

  # temporary files
  spladder_cleaned_gtf  <- file.path(spladder_out, "spladder_combined.cleaned.gtf")
  spladder_filtered_gtf <- file.path(spladder_out, "spladder_combined.filtered.gtf")
  spladder_prefixed_gtf <- file.path(spladder_out, "spladder_prefixed.gtf")
  merged_temp           <- file.path(spladder_out, "merged_temp.gtf")
  gtf_merged            <- file.path(spladder_out, "merged_sorted.gtf")

  # STEP 1: remove non-printable characters
  cat("Step 1/5: Cleaning non-printable characters (Perl)...\n")
  system(sprintf("perl -pe 's/[^[:print:]\\t\\n]//g' %s > %s", shQuote(spladder_gtf), shQuote(spladder_cleaned_gtf)))
  if (!file.exists(spladder_cleaned_gtf)) stop("Step 1 (Perl clean) failed.")

  # STEP 2: filter empty IDs
  cat("Step 2/5: Filtering empty IDs (Awk)...\n")
  system(sprintf("awk -F'\\t' 'BEGIN{OFS=\"\\t\"} $9 ~ /transcript_id \"\"/ || $9 ~ /gene_id \"\"/ {next} {print}' %s > %s",
                 shQuote(spladder_cleaned_gtf), shQuote(spladder_filtered_gtf)))
  if (!file.exists(spladder_filtered_gtf)) stop("Step 2 (Awk filter) failed.")

  # STEP 3: add 'spladder_' prefix
  cat("Step 3/5: Adding 'spladder_' prefix (Awk)...\n")
  system(sprintf(
    "awk -F'\\t' 'BEGIN{OFS=\"\\t\"} { if($3==\"transcript\" || $3==\"exon\"){ gsub(/transcript_id \"/, \"transcript_id \\\"spladder_\"); gsub(/gene_id \"/, \"gene_id \\\"spladder_\") } print }' %s > %s",
    shQuote(spladder_filtered_gtf), shQuote(spladder_prefixed_gtf)))
  if (!file.exists(spladder_prefixed_gtf)) stop("Step 3 (Awk prefix) failed.")

  # STEP 4: merge with the original GTF
  cat("Step 4/5: Merging with original GTF...\n")
  system(paste("cat", shQuote(gtf), shQuote(spladder_prefixed_gtf), ">", shQuote(merged_temp)))

  # STEP 5: sort
  cat("Step 5/5: Sorting final file...\n")
  system(sprintf("sort -k1,1 -k4,4n %s > %s", shQuote(merged_temp), shQuote(gtf_merged)))

  file.remove(spladder_gtf, spladder_cleaned_gtf, spladder_filtered_gtf, spladder_prefixed_gtf, merged_temp)

  cat("New GTF created successfully:", gtf_merged, "\n")
}


#------------------------------------------------------------------------------
# RSEM (reference, transcriptome alignment, quantification)
#------------------------------------------------------------------------------

if (RSEM == TRUE) {

  # --- RSEM reference ---
  if (transcriptome_indexing == TRUE) {

    print("Transcriptome indexing is TRUE. Starting RSEM reference creation...")

    gtf_file_path <- if (GTF_addition == TRUE) gtf_merged else gtf
    gtf_file_abs  <- tools::file_path_as_absolute(gtf_file_path)
    fasta_abs     <- tools::file_path_as_absolute(fasta)

    dir.create(rsem_ref_dir, recursive = TRUE, showWarnings = FALSE)
    rsem_prefix_basename <- basename(rsem_ref_dir)
    Star_dir <- dirname(Star_bin)

    original_wd <- getwd()
    setwd(rsem_ref_dir)
    print(paste("Creating RSEM reference inside:", getwd()))

    status <- system2("rsem-prepare-reference", args = c(
      "--gtf", gtf_file_abs,
      "--star",
      "--star-path", Star_dir,
      "--num-threads", CPU,
      fasta_abs,
      rsem_prefix_basename
    ))

    setwd(original_wd)
    print(paste("Returned to working directory:", getwd()))
    check_status(status, "rsem-prepare-reference")

  } else {
    print("Transcriptome indexing is FALSE. Skipping reference creation.")
    print(paste("Using pre-existing RSEM reference from:", rsem_ref_dir))
  }

  rsem_prefix_basename <- basename(rsem_ref_dir)
  rsem_prefix_full <- file.path(rsem_ref_dir, rsem_prefix_basename)

  if (!file.exists(paste0(rsem_prefix_full, ".grp"))) {
    warning("RSEM prefix not found at ", rsem_prefix_full,
            ". Please ensure the prefix name matches the directory name.")
  }

  transcriptome_out <- file.path(base_dir, "Transcriptome")
  table$bam_final_transcriptome <- NA_character_

  # --- Transcriptome alignment ---
  if (transcriptome_alignment == TRUE) {

    print("Transcriptome alignment is TRUE. Starting STAR alignment for RSEM...")

    for (i in seq_len(nrow(table))) {
      sample_ID <- table$sample_ID[i]
      fastq_files <- get_fastq_files(sample_ID, i)

      tr_bam_dir <- file.path(transcriptome_out, "BAM", sample_ID)
      dir.create(tr_bam_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(transcriptome_out, "logs", sample_ID), recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(transcriptome_out, "counts", sample_ID), recursive = TRUE, showWarnings = FALSE)

      star_prefix <- file.path(tr_bam_dir, sample_ID)

      if (mode == "single") {
        reads_in <- shQuote(fastq_files)
        zcat_flag <- if (needs_zcat(fastq_files)) " --readFilesCommand zcat" else ""
      } else {
        reads_in <- paste(shQuote(fastq_files[[1]]), shQuote(fastq_files[[2]]))
        zcat_flag <- if (needs_zcat(fastq_files[[1]])) " --readFilesCommand zcat" else ""
      }

      st <- system2(Star_bin, paste0(
        "--runThreadN ", CPU,
        " --genomeDir ", shQuote(rsem_ref_dir),
        " --readFilesIn ", reads_in,
        " --outSAMtype None",
        " --quantMode TranscriptomeSAM",
        zcat_flag,
        " --outFileNamePrefix ", shQuote(star_prefix)
      ))
      check_status(st, paste("STAR transcriptome alignment for", sample_ID))

      table$bam_final_transcriptome[i] <- paste0(star_prefix, "Aligned.toTranscriptome.out.bam")
    }

  } else {

    print("Transcriptome alignment is FALSE. Skipping alignment loop.")
    print("Using the bam_transcriptome column of the table.")

    table$bam_final_transcriptome <- normalizePath(as.character(table$bam_transcriptome))
    print("Transcriptome BAM file paths successfully located.")
  }

  # --- RSEM expression quantification ---
  for (i in seq_len(nrow(table))) {
    sample_ID <- table$sample_ID[i]
    counts_dir <- file.path(transcriptome_out, "counts", sample_ID)
    dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)

    rsem_args <- c(
      "--alignments",
      "--no-bam-output",
      "--num-threads", CPU,
      "--calc-pme",
      "--seed", 12345,
      table$bam_final_transcriptome[i],
      rsem_prefix_full,
      file.path(counts_dir, sample_ID)
    )
    if (mode == "paired") rsem_args <- c("--paired-end", rsem_args)

    check_status(system2("rsem-calculate-expression", args = rsem_args),
                 paste("rsem-calculate-expression for", sample_ID))
  }

  print("RSEM expression calculation finished for all samples")
}


#------------------------------------------------------------------------------
# IsoformSwitchAnalyzeR: import (RSEM counts), satuRn test, Part 1
#------------------------------------------------------------------------------

if (IsoformSwitch == TRUE) {

  library(IsoformSwitchAnalyzeR)
  if (!requireNamespace("satuRn", quietly = TRUE)) {
    stop("Package 'satuRn' is required for isoformSwitchTestSatuRn. Install it with BiocManager::install('satuRn').")
  }

  cat("\n--- Preparing data for IsoformSwitchAnalyzeR ---\n")

  isar_dir <- file.path(base_dir, "IsoformSwitchAnalyzer_Output")
  dir.create(isar_dir, showWarnings = FALSE, recursive = TRUE)

  if (!dir.exists(rsem_counts_dir)) stop("RSEM counts directory not found: ", rsem_counts_dir)
  cat("Using RSEM counts directory:", rsem_counts_dir, "\n")

  # --- 1. transcript FASTA of the RSEM reference ---
  if (is.na(rsem_ref_dir) || !dir.exists(rsem_ref_dir)) stop("RSEM reference directory not found: ", rsem_ref_dir)
  transcript_fasta_path <- file.path(rsem_ref_dir, paste0(basename(rsem_ref_dir), ".transcripts.fa"))
  cat("Using transcript FASTA:", transcript_fasta_path, "\n")
  if (!file.exists(transcript_fasta_path)) stop("Transcript FASTA file not found: ", transcript_fasta_path)

  # --- 2. GTF compatible with RSEM and IsoformSwitchAnalyzeR ---
  if (gtf_isoformswitch_exists == TRUE) {

    gtf_cleaned_path <- gtf_isoformswitch
    cat("Using existing GTF for IsoformSwitchAnalyzeR:", gtf_cleaned_path, "\n")

  } else {

    if (GTF_addition == TRUE) {
      gtf_source  <- gtf_merged
      gff3_source <- spladder_combined_gff
      cat("Using merged GTF from SplAdder:", gtf_source, "\n")

      gene_info_file <- file.path(rsem_ref_dir, "geneInfo.tab")
      cat("Using:", gene_info_file, "\n", gtf_source, "\n", gff3_source, "\n")

      if (file.exists(gene_info_file) && file.exists(gff3_source)) {

        gtf_obj <- rtracklayer::import(gtf_source)

        gene_map <- readr::read_tsv(gene_info_file, col_names = FALSE, show_col_types = FALSE, skip = 1) %>%
          dplyr::select(gene_id = 1, gene_symbol = 2)

        cat("  Parsing GFF3 to link gene name to splicing events...\n")
        gff_text <- readr::read_tsv(gff3_source, comment = "#", col_names = FALSE, show_col_types = FALSE)

        spladder_map <- gff_text %>%
          dplyr::select(X9) %>%
          dplyr::filter(stringr::str_detect(X9, "ID=") & stringr::str_detect(X9, "GeneName=")) %>%
          dplyr::mutate(
            raw_id = stringr::str_extract(X9, "ID=([^;]+)"),
            raw_id = stringr::str_remove(raw_id, "ID="),
            transcript_id = paste0("spladder_", raw_id),
            gene_id_raw = stringr::str_extract(X9, 'GeneName="([^"]+)"'),
            gene_id = stringr::str_remove_all(gene_id_raw, 'GeneName=|"')
          ) %>%
          dplyr::select(transcript_id, gene_id) %>%
          dplyr::distinct()

        complete_map <- spladder_map %>%
          dplyr::left_join(gene_map, by = "gene_id") %>%
          dplyr::filter(!is.na(gene_symbol)) %>%
          dplyr::mutate(new_gene_name = paste0("spladder_", gene_symbol)) %>%
          dplyr::select(transcript_id, new_gene_name)

        cat("Found", nrow(complete_map), "gene names for SplAdder events.\n")

        meta_df <- as.data.frame(S4Vectors::mcols(gtf_obj))
        if (!"gene_name" %in% colnames(meta_df)) meta_df$gene_name <- NA

        cat("  Adding gene names to the GTF...\n")
        meta_df <- meta_df %>%
          dplyr::left_join(complete_map, by = "transcript_id") %>%
          dplyr::mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", new_gene_name, gene_name)) %>%
          dplyr::select(-new_gene_name)

        S4Vectors::mcols(gtf_obj) <- meta_df

        gtf_enriched_path <- file.path(spladder_out, "merged_sorted_enriched.gtf")
        cat("  Exporting new GTF:", gtf_enriched_path, "\n")
        rtracklayer::export(gtf_obj, gtf_enriched_path, format = "GTF")
        gtf_per_isoformswitch <- gtf_enriched_path

      } else {
        warning("geneInfo.tab or spladder_combined.gff3 not found, skipping gene name addition.")
        gtf_per_isoformswitch <- gtf_source
      }

    } else {
      gtf_per_isoformswitch <- gtf
      cat("Using original GTF:", gtf_per_isoformswitch, "\n")
    }

    if (!file.exists(gtf_per_isoformswitch)) stop("GTF file for IsoformSwitchAnalyzeR not found: ", gtf_per_isoformswitch)

    # clean the GTF so that it is compatible with IsoformSwitchAnalyzeR
    out_gtf <- file.path(base_dir, "GTF_cleaned")
    dir.create(out_gtf, recursive = TRUE, showWarnings = FALSE)

    python_script_path <- "/opt/isoformSwitchAnalyzeR/modify_gtf.py"
    gtf_cleaned_path <- file.path(out_gtf, "isoform_switch_cleaned.gtf")
    if (!file.exists(python_script_path)) stop("Script Python 'modify_gtf.py' not found: ", python_script_path)

    cat("Running python script to clean GTF...\n")
    check_status(system2("python3", args = c(shQuote(python_script_path),
                                              shQuote(gtf_per_isoformswitch),
                                              shQuote(gtf_cleaned_path))),
                 "modify_gtf.py")
    if (!file.exists(gtf_cleaned_path)) stop("Cleaned GTF was not created by python script.")
    cat("GTF cleaned successfully:", gtf_cleaned_path, "\n")
  }

  # --- 3. import RSEM quantifications ---
  cat("\nImporting RSEM quantification...\n")
  data <- importIsoformExpression(parentDir = rsem_counts_dir)

  # importIsoformExpression() can return isoform_id as an ordinary first column
  # rather than as row names: convert it explicitly before subsetting samples.
  prepare_rsem_matrix <- function(x, matrix_name) {
    if (is.null(dim(x)) || length(dim(x)) != 2) {
      stop("RSEM ", matrix_name, " is not a two-dimensional matrix/data.frame.")
    }
    if ("isoform_id" %in% colnames(x)) {
      isoform_ids <- as.character(x[["isoform_id"]])
      if (anyNA(isoform_ids) || any(isoform_ids == "")) stop("Missing or empty isoform_id values in RSEM ", matrix_name, ".")
      if (anyDuplicated(isoform_ids)) stop("Duplicated isoform_id values in RSEM ", matrix_name, ".")
      x[["isoform_id"]] <- NULL
      x <- as.matrix(x)
      rownames(x) <- isoform_ids
    } else {
      rn <- rownames(x)
      if (is.null(rn) || identical(rn, as.character(seq_len(nrow(x))))) {
        stop("RSEM ", matrix_name, " has no isoform_id column and no meaningful isoform row names.")
      }
      x <- as.matrix(x)
    }
    x
  }

  data$abundance <- prepare_rsem_matrix(data$abundance, "abundance")
  data$counts    <- prepare_rsem_matrix(data$counts, "counts")

  if (!setequal(rownames(data$abundance), rownames(data$counts))) {
    stop("RSEM abundance and counts matrices contain different isoform IDs.")
  }
  data$counts <- data$counts[rownames(data$abundance), , drop = FALSE]

  valid_samples <- myDesign$sampleID
  rsem_samples <- colnames(data$abundance)
  if (is.null(rsem_samples) || is.null(colnames(data$counts))) {
    stop("Could not identify sample names in RSEM abundance/count matrices.")
  }
  missing_samples <- setdiff(valid_samples, rsem_samples)
  extra_samples <- setdiff(rsem_samples, valid_samples)
  if (length(missing_samples) > 0) {
    stop("Samples listed in table.xlsx but missing from RSEM quantifications: ", paste(missing_samples, collapse = ", "))
  }
  if (!setequal(colnames(data$abundance), colnames(data$counts))) {
    stop("RSEM abundance and counts matrices do not contain the same sample columns.")
  }
  if (length(extra_samples) > 0) {
    cat("\nThe following RSEM samples are not present in table.xlsx and will be ignored:\n",
        paste(extra_samples, collapse = ", "), "\n")
  }

  data$abundance <- data$abundance[, valid_samples, drop = FALSE]
  data$counts    <- data$counts[, valid_samples, drop = FALSE]
  stopifnot(identical(colnames(data$abundance), myDesign$sampleID))
  stopifnot(identical(colnames(data$counts), myDesign$sampleID))

  cat("Number of quantified isoforms:", nrow(data$abundance), "\n")
  cat("Example isoform IDs:", paste(head(rownames(data$abundance), 5), collapse = ", "), "\n")

  # --- 4. switchAnalyzeRlist ---
  cat("\nCreating switchAnalyzeRlist with the comparisons defined by cfr_1vs2...\n")
  switchList <- IsoformSwitchAnalyzeR::importRdata(
    isoformCountMatrix   = data$counts,
    isoformRepExpression = data$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = gtf_cleaned_path,
    isoformNtFasta       = transcript_fasta_path,
    comparisonsToMake    = myComparisons
  )

  cat("\nswitchAnalyzeRlist created successfully.\n")
  cat("Number of isoforms imported:", nrow(switchList$isoformFeatures), "\n")
  cat("Number of genes imported:", length(unique(switchList$isoformFeatures$gene_id)), "\n")

  output_rds_imported <- file.path(isar_dir, "switchAnalyzeRlist_imported.rds")
  saveRDS(switchList, file = output_rds_imported)
  cat("Imported object saved to:", output_rds_imported, "\n")

  # --- 5. differential isoform usage ---

  cat("\nRunning isoformSwitchTestDEXSeq...\n")
    switchList_tested <- IsoformSwitchAnalyzeR::isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = switchList,
    alpha = 0.05,
    reduceToSwitchingGenes = TRUE
  )



  output_rds_tested <- file.path(isar_dir, "switchAnalyzeRlist_tested.rds")
  saveRDS(switchList_tested, file = output_rds_tested)
  cat("Statistically tested object saved to:", output_rds_tested, "\n")

  if ("isoform_switch_q_value" %in% colnames(switchList_tested$isoformFeatures)) {
    cat("Non-NA isoform switching q-values:", sum(!is.na(switchList_tested$isoformFeatures$isoform_switch_q_value)), "\n")
  }
  if ("gene_switch_q_value" %in% colnames(switchList_tested$isoformFeatures)) {
    cat("Non-NA gene switching q-values:", sum(!is.na(switchList_tested$isoformFeatures$gene_switch_q_value)), "\n")
  }

  # --- 6. Part 1: sequence extraction and consequence preparation ---
  if (nrow(switchList_tested$isoformFeatures) == 0) {
    warning("No isoforms left after isoformSwitchTestSatuRn (no switching genes): skipping Part 1.")
  } else {
    cat("\n--- Starting IsoformSwitchAnalyzeR Part 1 ---\n")
    library(BSgenome.Hsapiens.UCSC.hg38)

    fasta_output_dir <- file.path(base_dir, "Analysis", "fasta")
    dir.create(fasta_output_dir, recursive = TRUE, showWarnings = FALSE)

    switchList_part1 <- isoformSwitchAnalysisPart1(
      switchAnalyzeRlist = switchList_tested,
      pathToGTF = gtf_cleaned_path,
      genomeObject = NULL,
      outputSequences = TRUE,
      pathToOutput = fasta_output_dir,
      prepareForWebServers = FALSE
    )

    output_rds_file_p1 <- file.path(isar_dir, "switchList_part1.rds")
    saveRDS(switchList_part1, file = output_rds_file_p1)
    cat("switchList_part1 saved to:", output_rds_file_p1, "\n")
  }

  cat("\n--- IsoformSwitchAnalyzeR workflow completed ---\n")

} else {
  cat("\nIsoformSwitch set to FALSE. Skipping IsoformSwitchAnalyzeR import, testing and Part 1.\n")
}


#------------------------------------------------------------------------------
# Gene table output filtering
#------------------------------------------------------------------------------

if (Rmats || Spladder_any || RSEM) {

  cat("\n--- Starting Output Filtering ---\n")

  # common thresholds
  fdr_cutoff <- 0.05
  pval_cutoff <- 0.01
  psi_diff_cutoff <- 0.05
  jc_min <- 5
  tpm_cutoff <- 1
  count_cutoff <- 10

  # -----------------------------------------------------
  # 1. rMATS filtering (per cell line)
  # -----------------------------------------------------
  if (Rmats == TRUE) {
    cat("Processing rMATS output...\n")

    event_types_rmats <- c("SE", "A3SS", "A5SS", "MXE", "RI")

    get_JCcounts <- function(x) {
      vals <- suppressWarnings(as.numeric(strsplit(as.character(x), ",")[[1]]))
      sum(vals[!is.na(vals)])
    }

    for (spec in cell_line_specs) {
      cl_safe <- safe_name(spec$cell_line)
      input_dir_rmats <- file.path(base_dir, "rmatsout", cl_safe, "OD")
      output_dir_rmats <- file.path(base_dir, "rmatsout", cl_safe, "Filtered_rMATS")
      dir.create(output_dir_rmats, showWarnings = FALSE, recursive = TRUE)

      if (!dir.exists(input_dir_rmats)) {
        warning("rMATS output directory not found: ", input_dir_rmats, ". Skipping rMATS filtering.")
        next
      }

      for (etype in event_types_rmats) {
        file_path <- file.path(input_dir_rmats, paste0(etype, ".MATS.JC.txt"))
        if (!file.exists(file_path)) {
          cat("  [", spec$cell_line, "] File not found:", basename(file_path), "\n")
          next
        }
        cat("  [", spec$cell_line, "] Filtering:", basename(file_path), "\n")

        tryCatch({
          df <- read_tsv(file_path, show_col_types = FALSE)

          required_cols_rmats <- c("GeneID", "geneSymbol", "IJC_SAMPLE_1", "SJC_SAMPLE_1",
                                   "IJC_SAMPLE_2", "SJC_SAMPLE_2", "PValue", "FDR", "IncLevelDifference")
          if (!all(required_cols_rmats %in% colnames(df))) {
            stop("One or more required rMATS columns missing. Found: ", paste(colnames(df), collapse = ", "))
          }

          # column 12 is the event ID (the file has two columns named "ID")
          df_filtered <- df %>%
            dplyr::select(GeneID, geneSymbol, ID = 12,
                          IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2,
                          PValue, FDR, IncLevelDifference) %>%
            dplyr::mutate(
              IncCount1 = sapply(IJC_SAMPLE_1, get_JCcounts),
              SkipCount1 = sapply(SJC_SAMPLE_1, get_JCcounts),
              IncCount2 = sapply(IJC_SAMPLE_2, get_JCcounts),
              SkipCount2 = sapply(SJC_SAMPLE_2, get_JCcounts),
              TotalCount = IncCount1 + SkipCount1 + IncCount2 + SkipCount2
            ) %>%
            dplyr::filter(
              !is.na(FDR) & FDR < fdr_cutoff,
              !is.na(IncLevelDifference) & abs(IncLevelDifference) > psi_diff_cutoff,
              !is.na(TotalCount) & TotalCount > jc_min,
              !is.na(PValue) & PValue < pval_cutoff
            ) %>%
            dplyr::select(GeneID, geneSymbol, ID, PValue, FDR, IncLevelDifference, TotalCount)

          if (nrow(df_filtered) > 0) {
            write_tsv(df_filtered, file.path(output_dir_rmats, paste0(etype, "_filtered.tsv")))
            cat("    Filtered results saved. Kept", nrow(df_filtered), "events.\n")
          } else {
            cat("    No events passed filtering.\n")
          }
        }, error = function(e) {
          cat("    Error processing file:", basename(file_path), "-", e$message, "\n")
        })
      }
    }
    cat("rMATS filtering complete.\n")
  } else {
    cat("\nRmats set to FALSE. Skipping rMATS filtering.\n")
  }

  # -----------------------------------------------------
  # 2. SplAdder filtering (per cell line)
  # -----------------------------------------------------
  if (Spladder_any) {
    cat("Processing SplAdder output...\n")

    for (sp in spladder_specs) {
      input_dir_spladder <- file.path(spladder_out, sp$testcase)
      output_dir_spladder <- file.path(spladder_out, "Filtered_SplAdder", sp$cl_safe)
      dir.create(output_dir_spladder, showWarnings = FALSE, recursive = TRUE)

      if (!dir.exists(input_dir_spladder)) {
        warning("SplAdder output directory not found: ", input_dir_spladder, ". Skipping SplAdder filtering.")
        next
      }

      spladder_files <- list.files(input_dir_spladder, pattern = "^test_results_C3_.*\\.tsv$", full.names = TRUE)
      spladder_files <- spladder_files[!grepl("\\.gene_unique\\.", spladder_files)]
      if (length(spladder_files) == 0) {
        cat("  [", sp$cell_line, "] No SplAdder result files (*_C3_*.tsv) found in", input_dir_spladder, "\n")
      }

      for (file_path in spladder_files) {
        cat("  [", sp$cell_line, "] Filtering:", basename(file_path), "\n")

        event_type_spladder <- str_match(basename(file_path), "test_results_C3_(.+)\\.tsv")[, 2]
        if (is.na(event_type_spladder)) event_type_spladder <- "unknown"

        tryCatch({
          df <- read_tsv(file_path, show_col_types = FALSE)

          required_cols <- c("p_val_adj", "dPSI", "mean_event_count_A", "mean_event_count_B", "p_val")
          if (!all(required_cols %in% colnames(df))) {
            warning("    Skipping file ", basename(file_path), ": missing one or more required columns (",
                    paste(required_cols, collapse = ", "), "). Available columns: ", paste(colnames(df), collapse = ", "))
            next
          }

          df_filtered <- df %>%
            dplyr::mutate(TotalMeanCount = mean_event_count_A + mean_event_count_B) %>%
            dplyr::filter(!is.na(p_val) & p_val < pval_cutoff) %>%
            dplyr::select(any_of(c("gene_id", "gene_name", "event_id", "chrm", "exon_pos",
                                   "p_val", "p_val_adj", "dPSI", "TotalMeanCount")))

          if (nrow(df_filtered) > 0) {
            write_tsv(df_filtered, file.path(output_dir_spladder, paste0(event_type_spladder, "_filtered.tsv")))
            cat("    Filtered results saved. Kept", nrow(df_filtered), "events.\n")
          } else {
            cat("    No events passed filtering.\n")
          }
        }, error = function(e) {
          cat("    Error processing file:", basename(file_path), "-", e$message, "\n")
        })
      }
    }
    cat("SplAdder filtering complete.\n")

    # move folder if needed

    # Test and Viz: the results were written inside the graph folder, move them to this run
    if (Spladder_mode == "test_viz") {
      dest <- file.path(base_dir, "spladder")
      dir.create(file.path(dest, "Filtered_SplAdder"), recursive = TRUE, showWarnings = FALSE)
      for (sp in spladder_specs) {
        # test results + viz plots (the plot folder is inside the testing_* folder)
        src <- file.path(spladder_out, sp$testcase)
        if (dir.exists(src)) check_status(system2("mv", c(shQuote(src), shQuote(dest))), paste("mv", sp$testcase))
        # filtered tables
        src <- file.path(spladder_out, "Filtered_SplAdder", sp$cl_safe)
        if (dir.exists(src)) check_status(system2("mv", c(shQuote(src), shQuote(file.path(dest, "Filtered_SplAdder")))), paste("mv filtered", sp$cl_safe))
      }
    }




  } else {
    cat("\nSpladder set to False. Skipping SplAdder filtering.\n")
  }

  # -----------------------------------------------------
  # 3. RSEM filtering
  # -----------------------------------------------------
  rsem_filter_input <- if (RSEM) file.path(base_dir, "Transcriptome", "counts") else rsem_counts_dir

  if ((RSEM || IsoformSwitch) && !is.na(rsem_filter_input) && dir.exists(rsem_filter_input)) {
    cat("Processing RSEM output...\n")

    output_dir_rsem <- file.path(base_dir, "Transcriptome", "Filtered_RSEM")
    dir.create(output_dir_rsem, showWarnings = FALSE, recursive = TRUE)

    sample_dirs <- list.dirs(rsem_filter_input, full.names = TRUE, recursive = FALSE)

    if (length(sample_dirs) == 0) {
      warning("No sample directories found in: ", rsem_filter_input, ". Skipping RSEM filtering.")
    } else {
      filter_rsem_file <- function(file, out_file, label) {
        if (!file.exists(file)) { cat("    File not found:", basename(file), "\n"); return(invisible(NULL)) }
        tryCatch({
          df <- read_tsv(file, show_col_types = FALSE)
          if (!all(c("expected_count", "TPM") %in% colnames(df))) {
            stop("Columns 'expected_count' or 'TPM' missing in ", label, " results.")
          }
          df_f <- df %>%
            dplyr::filter(!is.na(expected_count) & expected_count > count_cutoff,
                          !is.na(TPM) & TPM > tpm_cutoff)
          if (nrow(df_f) > 0) {
            write_tsv(df_f, out_file)
            cat("    Filtered", label, "results saved. Kept", nrow(df_f), label, "\n")
          } else {
            cat("    No", label, "passed filtering.\n")
          }
        }, error = function(e) {
          cat("    Error processing file:", basename(file), "-", e$message, "\n")
        })
      }

      for (sample_dir in sample_dirs) {
        sample_id <- basename(sample_dir)
        cat("  Filtering RSEM results for sample:", sample_id, "\n")
        sample_output_dir <- file.path(output_dir_rsem, sample_id)
        dir.create(sample_output_dir, showWarnings = FALSE, recursive = TRUE)

        filter_rsem_file(file.path(sample_dir, paste0(sample_id, ".genes.results")),
                         file.path(sample_output_dir, paste0(sample_id, ".genes.results_filtered.tsv")), "genes")
        filter_rsem_file(file.path(sample_dir, paste0(sample_id, ".isoforms.results")),
                         file.path(sample_output_dir, paste0(sample_id, ".isoforms.results_filtered.tsv")), "isoforms")
      }
    }
    cat("RSEM filtering complete.\n")
  } else {
    cat("\nNo RSEM output to filter. Skipping RSEM filtering.\n")
  }

  cat("\n--- All selected filtering steps finished! ---\n")

} else {
  cat("\n--- No filtering steps selected (Rmats, Spladder, RSEM are all off) ---\n")
}
