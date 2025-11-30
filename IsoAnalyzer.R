#------------------------------------------------------------------------------
#BULK RNA SEQ
#------------------------------------------------------------------------------

#packages
library(tidyr)
library(stringr)
library(grid)
library(readxl)
library(dplyr)
library(readr)
library(ggplot2)




#set up environment
theme_set(theme_bw(12) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face = "bold", margin = margin(10, 0, 10, 0)),
                  axis.text.x = element_text(angle = 45, hjust = 1)))
options(stringsAsFactors = FALSE)
update_geom_defaults("point", aes(size = 4))
set.seed(1234)

pal <- c("#3283FE", "#FA0087", "#009E73", "#FBE426", "#56B4E9", "#FEAF16", "#DEA0FD", "#1CBE4F", "#F6222E", "#1CFFCE", "#325A9B", "#AA0DFE","#D55E00", "#2ED9FF", "#f0E442", "#1C8356", "#0072B2", "#CC79A7")

current_time <- "251106"  #format(Sys.time(), "%d%m%y_%H%M%S")
dir.create(paste0("./Output_", current_time))


#upload settings

print("uploading settings file")

settings <- as.data.frame(read_excel("./settings.xlsx", col_names = FALSE))
rownames(settings) <- settings$...1
settings <- settings %>% dplyr::select(-c("...1", "...2"))
settings <- as.data.frame(t(settings))


if(length(unique(na.omit(settings$trimming))) != 1 || is.na(as.logical(settings$trimming[1]))){
  print("invalid argument for trimming"); stop()
} else {
  trimming = as.logical(settings$trimming[1])
}

if(trimming == TRUE){
  
  # fastqc
  if(length(unique(na.omit(settings$fastqc))) != 1 || is.na(as.logical(settings$fastqc[1]))){
    print("invalid argument for fastqc"); stop()
  } else {
    fastqc = as.logical(settings$fastqc[1])
    fqc = if (fastqc) " --fastqc" else ""
  }

  # mode
  if(length(unique(na.omit(settings$mode))) != 1 || !as.character(settings$mode[1]) %in% c("single", "paired")){
    print("invalid argument for mode"); stop()
  } else {
    mode = as.character(settings$mode[1])
  }

  # adapter 1
  adapter = as.character(settings$adapter[1])
  if(!is.na(adapter) && nchar(adapter) > 0){
    bases <- strsplit(adapter, "")[[1]]
    invalid_bases = bases[!bases %in% c("A", "T", "G", "C")]
    if(length(unique(na.omit(settings$adapter))) != 1 || length(invalid_bases) != 0) {
      print("invalid argument for adapter"); stop()
    }
  }

  # clip R1
  if(is.null(settings$clip_5_1[1]) || is.na(settings$clip_5_1[1]) || settings$clip_5_1[1] == ""){
    clip_5_1 = 0
  } else {
    clip_5_1 = as.integer(settings$clip_5_1[1])
    if(is.na(clip_5_1)){ print("invalid argument for clip_5_1"); stop() }
  }

  if(is.null(settings$clip_3_1[1]) || is.na(settings$clip_3_1[1]) || settings$clip_3_1[1] == ""){
    clip_3_1 = 0
  } else {
    clip_3_1 = as.integer(settings$clip_3_1[1])
    if(is.na(clip_3_1)){ print("invalid argument for clip_3_1"); stop() }
  }

  # paired-end specific
  if(mode == "paired"){
    
    # adapter2
    adapter2 = as.character(settings$adapter2[1])
    if(!is.na(adapter2) && nchar(adapter2) > 0){
      bases = strsplit(adapter2, "")[[1]]
      invalid_bases = bases[!bases %in% c("A", "T", "G", "C")]
      if(length(unique(na.omit(settings$adapter2))) != 1 || length(invalid_bases) != 0){
        print("invalid argument for adapter2"); stop()
      }
    }

    # clip R2
    if(is.null(settings$clip_5_2[1]) || is.na(settings$clip_5_2[1]) || settings$clip_5_2[1] == ""){
      clip_5_2 = 0
    } else {
      clip_5_2 = as.integer(settings$clip_5_2[1])
      if(is.na(clip_5_2)){ print("invalid argument for clip_5_2"); stop() }
    }

    if(is.null(settings$clip_3_2[1]) || is.na(settings$clip_3_2[1]) || settings$clip_3_2[1] == ""){
      clip_3_2 = 0
    } else {
      clip_3_2 = as.integer(settings$clip_3_2[1])
      if(is.na(clip_3_2)){ print("invalid argument for clip_3_2"); stop() }
    }
  }
}


if(length(unique(na.omit(settings$UMI))) != 1 || !as.logical(settings$UMI[1]) != "NA"){
  print("invalid argument for UMI"); stop()
} else{UMI = as.logical(settings$UMI[1])}

if(UMI == "TRUE"){
  if(length(unique(na.omit(settings$umi_len))) != 1 || as.integer(settings$umi_len[1]) == "NA"){
    print("invalid argument for umi_len"); stop()
  } else{umi_len = as.integer(settings$umi_len[1])}
  
  if(length(unique(na.omit(settings$umi_loc))) != 1 || !as.character(settings$umi_loc[1]) %in% c("per_read", "read1", "read2", "per_index")){
    print("invalid argument for umi_loc"); stop()
  } else{umi_loc = as.character(settings$umi_loc[1])}
  
  if(length(unique(na.omit(settings$umi_prefix))) != 1 || as.character(settings$umi_prefix[1]) == "NA"){
    print("invalid argument for umi_prefix"); stop()
  } else{umi_prefix = as.character(settings$umi_prefix[1])}

  if(length(unique(na.omit(settings$umi_config))) != 1 || is.na(settings$umi_config[1]) || settings$umi_config[1] == ""){
    print("invalid argument for umi_config"); stop()
  } else {
    umi_config = as.character(settings$umi_config[1])
  }
}


if(length(unique(na.omit(settings$indexing))) != 1 || as.logical(settings$indexing[1]) == "NA"){
  print("invalid argument for indexing"); stop()
} else{indexing = as.logical(settings$indexing[1])}

if(length(unique(na.omit(settings$alignment))) != 1 || as.logical(settings$alignment[1]) == "NA"){
  print("invalid argument for alignment"); stop()
} else{alignment = as.logical(settings$alignment[1])}


if (length(unique(na.omit(settings$gtf))) != 1 || !file.exists(settings$gtf[1])) {
    print("invalid argument for gtf"); stop()
  } else {
    gtf <- as.character(settings$gtf[1])
  }


# --- RSEM settings ---
if(length(unique(na.omit(settings$RSEM))) != 1 || is.na(settings$RSEM[1])){
  print("invalid argument for RSEM"); stop()
} else {
  RSEM <- as.logical(settings$RSEM[1])
}


if(length(unique(na.omit(settings$transcriptome_indexing))) != 1 || is.na(as.logical(settings$transcriptome_indexing[1]))){
  print("invalid argument for transcriptome_indexing"); stop()
} else {
  transcriptome_indexing = as.logical(settings$transcriptome_indexing[1])
}

if(transcriptome_indexing == FALSE) {
 
  if (length(unique(na.omit(settings$transcriptome_genome_dir))) != 1 || !dir.exists(settings$transcriptome_genome_dir[1])) {
    print("invalid argument for transcriptome_genome_dir (required when transcriptome_indexing is FALSE)"); stop()
    } else {
    transcriptome_genome_dir <- as.character(settings$transcriptome_genome_dir[1])
    }
} else {
  transcriptome_genome_dir <- NA 
}

# --- RSEM transcriptome alignment settings ---
if(length(unique(na.omit(settings$transcriptome_alignment))) != 1 || is.na(as.logical(settings$transcriptome_alignment[1]))){
  print("invalid argument for transcriptome_alignment"); stop()
} else {
  transcriptome_alignment = as.logical(settings$transcriptome_alignment[1])
}

if(transcriptome_alignment == FALSE) {
  
  if (length(unique(na.omit(settings$bam_dir_transcriptome))) != 1 || !dir.exists(settings$bam_dir_transcriptome[1])) {
    print("invalid argument for bam_dir_transcriptome (required when transcriptome_alignment is FALSE)"); stop()
    } else {
      bam_dir_transcriptome <- as.character(settings$bam_dir_transcriptome[1])
    }
} else {
  bam_dir_transcriptome <- NA 
}



if (indexing == TRUE || RSEM == TRUE) {
  
  if (length(unique(na.omit(settings$fasta))) != 1 || !file.exists(settings$fasta[1])) {
    print("invalid argument for fasta"); stop()
  } else {
    fasta <- as.character(settings$fasta[1])
  }

} else if (alignment == TRUE) {

  if (length(unique(na.omit(settings$genome_dir))) != 1 || !dir.exists(settings$genome_dir[1])) {
    print("invalid argument for genome_dir"); stop()
  } else {
    genome_dir <- as.character(settings$genome_dir[1])
  }

  if (length(unique(na.omit(settings$gtf))) != 1 || !file.exists(settings$gtf[1])) {
    print("invalid argument for gtf"); stop()
  } else {
    gtf <- as.character(settings$gtf[1])
  }

}

if(length(unique(na.omit(settings$CPU))) != 1 || as.integer(settings$CPU[1]) == "NA"){
    print("invalid argument for CPU"); stop()
  } else{CPU = as.integer(settings$CPU[1])}



if (alignment == TRUE ) {

  if (length(unique(na.omit(settings$gene_length))) != 1 || !file.exists(settings$gene_length[1])) {
    print("invalid argument for gene_length"); stop()
  } else {
    gene_length <- as.character(settings$gene_length[1])
  }

  if (length(unique(na.omit(settings$gene_info))) != 1 || !file.exists(settings$gene_info[1])) {
    print("invalid argument for gene_info"); stop()
  } else {
    geneInfo <- as.character(settings$gene_info[1])
    geneInfo <- read.table(geneInfo, quote = "\"", comment.char = "", skip = 1)
  }

}



if(alignment == TRUE || RSEM == TRUE){if(length(unique(na.omit(settings$mode))) != 1 || !as.character(settings$mode[1]) %in% c("single", "paired")){
  print("invalid argument for mode"); stop()
} else{mode = as.character(settings$mode[1])}}


if(length(unique(na.omit(settings$Rmats))) != 1 || !as.logical(settings$Rmats[1]) != "NA"){
  print("invalid argument for Rmats"); stop()
} else{Rmats = as.logical(settings$Rmats[1])}

if (Rmats == TRUE && alignment == FALSE){
   if(length(unique(na.omit(settings$rmats_path))) != 1 || as.character(settings$rmats_path[1]) == "NA"){
     print("invalid argument for rmats_path"); stop()
   } else{rmats_path = as.character(settings$rmats_path[1])}

   if (length(unique(na.omit(settings$bam_dir))) != 1 || !dir.exists(settings$bam_dir[1])) {
     stop("Invalid argument for bam_dir")
   } else {bam_dir <- as.character(settings$bam_dir[1])}

}


if(length(unique(na.omit(settings$Spladder))) != 1 || !as.logical(settings$Spladder[1]) != "NA"){
  print("invalid argument for Spladder"); stop()
} else{Spladder = as.logical(settings$Spladder[1])}

if (Spladder == TRUE && alignment == FALSE){ 
   if (length(unique(na.omit(settings$txt_path))) != 1 || !file.exists(settings$txt_path[1])) {
     print("invalid argument for txt_path"); stop()
   } else { txt_path <- as.character(settings$txt_path[1]) }

   if (length(unique(na.omit(settings$bam_dir))) != 1 || !dir.exists(settings$bam_dir[1])) {
    stop("Invalid argument for bam_dir")
   } else { bam_dir <- as.character(settings$bam_dir[1])
}

}




if(length(unique(na.omit(settings$GTF_addition))) != 1 || is.na(as.logical(settings$GTF_addition[1]))){
  print("invalid argument for GTF_addition"); stop()
} else {
  GTF_addition = as.logical(settings$GTF_addition[1])
}







#Analysis parameters

params <- list(
  
  UMI_dedup = if("UMI_dedup" %in% colnames(settings)) {
                val = settings$UMI_dedup[1]
                ifelse(is.na(val), FALSE, as.logical(val))
              } else {
                FALSE
              }
)

if(alignment == TRUE){ 
  print("alignment = TRUE")
}

#table

print("uploading table file")

table = read_excel("./table.xlsx")
if(alignment == TRUE){
  if(mode == "single"){
    if(!"fastq" %in% colnames(table)){
      print("fastq column is missing"); stop()
    } else{
      for(i in table$fastq){
        if(file.exists(i) == FALSE){
          print("incorrect path in fastq column"); stop()
        }
      }
    }
  } else if(mode == "paired"){
    if(!"fastq1" %in% colnames(table) || !"fastq2" %in% colnames(table)){
      print("fastq1 or fastq2 columns are missing"); stop()
    } else{
      for(i in table$fastq1){if(file.exists(i) == FALSE){
        print("incorrect path in fastq1 column"); stop()}
      }
      for(i in table$fastq2){if(file.exists(i) == FALSE){
        print("incorrect path in fastq2 column"); stop()}
      }
    }
  }
} 


#trimming


if(trimming == TRUE){
  print("starting trimming")
  dir.create(paste0("./Output_", current_time, "/trimmed"))
  out_dir = paste0("./Output_", current_time, "/trimmed")
  
  if(mode == "single"){
    for(i in 1:length(rownames(table))){
      sample_ID = table$sample_ID[i]
      fastq = table$fastq[i]
      
      if (UMI == TRUE){
        cleaned_fastq = paste0(out_dir, "/", sample_ID, "_umi.fq.gz")
        system2("fastp", paste0(
          "-i ", fastq, " ",
          "-o ", cleaned_fastq, " ",
          "--umi ",
          "--umi_loc ", umi_loc, " ",
          "--umi_len ", umi_len, " ",
          "--dedup ",
          "--umi_prefix ", umi_prefix, " ",
          ifelse(is.na(adapter) || adapter == "", "--detect_adapter_for_pe", paste0("--adapter_sequence ", adapter)),
          " --thread ", CPU,
          " --html ", out_dir, "/", sample_ID, "_fastp_report.html",
          " --json ", out_dir, "/", sample_ID, "_fastp_report.json"
        ))
        fastq = cleaned_fastq
      } else {
        if(is.na(adapter) || adapter == ""){
          cleaned_fastq = paste0(out_dir, "/", sample_ID, "_trimmed.fq.gz")
          system2("fastp", paste0(
            "-i ", fastq, " ",
            "-o ", cleaned_fastq, " ",
            "--detect_adapter_for_pe ",
            ifelse(clip_5_1 == 0, "", paste0("--trim_front1 ", clip_5_1, " ")),
            ifelse(clip_3_1 == 0, "", paste0("--trim_tail1 ", clip_3_1, " ")),
            "--thread ", CPU,
            " --html ", out_dir, "/", sample_ID, "_fastp_report.html",
            " --json ", out_dir, "/", sample_ID, "_fastp_report.json"
          ))
          fastq = cleaned_fastq
        } else {
          if(clip_3_1 == 0){ clipping_R1_3 = "" } else { clipping_R1_3 = paste0(" --three_prime_clip_R1 ", clip_3_1) }
          if(clip_5_1 == 0){ clipping_R1_5 = "" } else { clipping_R1_5 = paste0(" --clip_R1 ", clip_5_1) }
          
          system2("/usr/bin/trim_galore", paste0(
            "--adapter ", adapter,
            clipping_R1_5, clipping_R1_3,
            " --gzip", fqc,
            " --basename ", sample_ID,
            " --output_dir ", out_dir, " ",
            fastq
          ))
        }
      }
    }
  } else if(mode == "paired"){
    for(i in 1:length(rownames(table))){
      sample_ID = table$sample_ID[i]
      fastq1 = table$fastq1[i]
      fastq2 = table$fastq2[i]
      
      if (UMI == TRUE){
        cleaned_fastq1 = paste0(out_dir, "/", sample_ID, "_umi_R1.fq.gz")
        cleaned_fastq2 = paste0(out_dir, "/", sample_ID, "_umi_R2.fq.gz")
        
        system2("fastp", paste0(
          "-i ", fastq1, " ",
          "-I ", fastq2, " ",
          "-o ", cleaned_fastq1, " ",
          "-O ", cleaned_fastq2, " ",
          "--umi ",
          "--dedup ",
          "--umi_loc ", umi_loc, " ",
          "--umi_len ", umi_len, " ",
          "--umi_prefix ", umi_prefix, " ",
          ifelse((is.na(adapter) || adapter == "") && (is.na(adapter2) || adapter2 == ""), 
                 "--detect_adapter_for_pe", 
                 paste0("--adapter_sequence ", adapter, " --adapter_sequence_r2 ", adapter2)),
          " --thread ", CPU,
          " --html ", out_dir, "/", sample_ID, "_fastp_report.html",
          " --json ", out_dir, "/", sample_ID, "_fastp_report.json"
        ))
        fastq1 = cleaned_fastq1
        fastq2 = cleaned_fastq2
      } else {
        if((is.na(adapter) || adapter == "") && (is.na(adapter2) || adapter2 == "")){
          cleaned_fastq1 = paste0(out_dir, "/", sample_ID, "_val_1.fq.gz")
          cleaned_fastq2 = paste0(out_dir, "/", sample_ID, "_val_2.fq.gz")
          
          system2("fastp", paste0(
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
            " --html ", out_dir, "/", sample_ID, "_fastp_report.html",
            " --json ", out_dir, "/", sample_ID, "_fastp_report.json"
          ))
          fastq1 = cleaned_fastq1
          fastq2 = cleaned_fastq2
        } else {
          if(clip_3_1 == 0){ clipping_R1_3 = "" } else { clipping_R1_3 = paste0(" --three_prime_clip_R1 ", clip_3_1) }
          if(clip_5_1 == 0){ clipping_R1_5 = "" } else { clipping_R1_5 = paste0(" --clip_R1 ", clip_5_1) }
          if(clip_3_2 == 0){ clipping_R2_3 = "" } else { clipping_R2_3 = paste0(" --three_prime_clip_R2 ", clip_3_2) }
          if(clip_5_2 == 0){ clipping_R2_5 = "" } else { clipping_R2_5 = paste0(" --clip_R2 ", clip_5_2) }
          
          system2("/usr/bin/trim_galore", paste0(
            "--adapter ", adapter,
            " --adapter2 ", adapter2,
            clipping_R1_5, clipping_R1_3,
            clipping_R2_5, clipping_R2_3,
            " --gzip --paired", fqc,
            " --basename ", sample_ID,
            " --output_dir ", out_dir, " ",
            fastq1, " ", fastq2
          ))
        }
      }
    }
  }
  print("finishing trimming")
}



#indexing

Star_bin <- Sys.which("STAR")
  if(nchar(Star_bin) == 0) stop("STAR non trovato nel PATH")


if(indexing == TRUE){
  print("starting indexing")
  
  system2("mkdir", "/home/genome")
  system2("chmod", "777 /home/genome")
  system2(Star_bin, paste0("--runThreadN ", CPU, " --runMode genomeGenerate --limitGenomeGenerateRAM=150000000000 --genomeDir /home/genome --genomeFastaFiles ", fasta, " --sjdbGTFfile ", gtf));
  system2("rm", paste0("-r ./Output_", current_time, "/genome"));
  system2("mv", paste0("/home/genome ./Output_", current_time));
  genome_dir = paste0("./Output_", current_time, "/genome");
  system2("gtftools", paste0("-l ./Output_", current_time, "/gene_length.txt ", gtf));
  gene_length = paste0("./Output_", current_time, "/gene_length.txt");
  geneInfo = read.table(paste0(genome_dir, "/geneInfo.tab"), quote = "\"", comment.char = "", skip = 1)
  
  print("finishing indexing")

}



#alignment



get_fastq_files <- function(sample_ID, mode, out_dir, trimming, UMI, adapter = NA, adapter2 = NA, i = NULL, table = NULL) {
  
  if (mode == "single") {
    
    if (UMI) {
      suffix <- "_umi.fq.gz"
    } else {
      if (is.na(adapter) || adapter == "") {
        suffix <- "_trimmed.fq.gz"
      } else {
        suffix <- "_trimmed.fq.gz" 
      }
    }
    
    
    return(if (trimming) paste0(out_dir, "/", sample_ID, suffix) else table$fastq[i])
    
  } else if (mode == "paired") {
    
    if (UMI) {
      suffix1 <- "_umi_R1.fq.gz"
      suffix2 <- "_umi_R2.fq.gz"
    } else {
      suffix1 <- "_val_1.fq.gz"
      suffix2 <- "_val_2.fq.gz"
    }
    
    return(list(
      if (trimming) paste0(out_dir, "/", sample_ID, suffix1) else table$fastq1[i],
      if (trimming) paste0(out_dir, "/", sample_ID, suffix2) else table$fastq2[i]
    ))
  }
}

if (alignment == TRUE) {
  print("starting alignment stage")

  
  bam_parent_dir <- paste0("./Output_", current_time, "/BAM")
  counts_parent_dir <- paste0("./Output_", current_time, "/counts")
  logs_parent_dir <- paste0("./Output_", current_time, "/logs")
  dir.create(bam_parent_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(counts_parent_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(logs_parent_dir, showWarnings = FALSE, recursive = TRUE)

  
  # --- Loop 1: STAR Alignment and Indexing ---
  print("--- Starting Loop 1: STAR Alignment & Indexing ---")
  table$bam_aligned_genomic <- NA 

  needs_zcat <- function(file) {
      
      if (is.null(file) || !is.character(file) || length(file) == 0 || is.na(file) || nchar(file) == 0) return(FALSE)
      endsWith(file, ".gz")
  }

  for (i in 1:nrow(table)) {
    sample_ID <- table$sample_ID[i]
    cat("  Processing sample:", sample_ID, "(", i, "/", nrow(table), ")\n")

    
    bam_dir <- file.path(bam_parent_dir, sample_ID) 
    log_dir <- file.path(logs_parent_dir, sample_ID)
    counts_dir <- file.path(counts_parent_dir, sample_ID) 
    dir.create(bam_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(counts_dir, showWarnings = FALSE, recursive = TRUE) 

   
    star_prefix <- file.path(bam_dir, sample_ID)

   
    fastq_files_list <- get_fastq_files(sample_ID, mode, out_dir, trimming, UMI, i = i, table = table)

    
    cat("    Running STAR alignment...\n")
    if (mode == "single") {
        input_file <- fastq_files_list 
        zcat_flag <- if (needs_zcat(input_file)) " --readFilesCommand zcat" else ""

        status_star <- system2(Star_bin, paste0(
                "--runThreadN ", CPU,
                " --genomeDir ", shQuote(genome_dir),
                " --readFilesIn ", shQuote(input_file),
                " --outSAMtype BAM SortedByCoordinate",
                zcat_flag,
                " --outFileNamePrefix ", shQuote(star_prefix)
        ))
    } else { 
        input_file1 <- fastq_files_list[[1]]
        input_file2 <- fastq_files_list[[2]]
        zcat_flag <- if (needs_zcat(input_file1)) " --readFilesCommand zcat" else ""

        status_star <- system2(Star_bin, paste0(
                "--runThreadN ", CPU,
                " --genomeDir ", shQuote(genome_dir),
                " --readFilesIn ", shQuote(input_file1), " ", shQuote(input_file2),
                " --outSAMtype BAM SortedByCoordinate",
                zcat_flag,
                " --outFileNamePrefix ", shQuote(star_prefix)
        ))
    }
    

    if(status_star != 0) stop("STAR alignment failed for sample: ", sample_ID)

   
    bam_aligned_path <- paste0(star_prefix, "Aligned.sortedByCoord.out.bam")
    if(!file.exists(bam_aligned_path)) stop("Aligned BAM not found after STAR run for sample: ", sample_ID)

    cat("    Running samtools index...\n")
    status_index <- system2("samtools", paste("index", shQuote(bam_aligned_path)))
    if(status_index != 0) stop("Samtools index failed for sample: ", sample_ID)
    print(paste("  Finished STAR and indexing for", sample_ID))

    
    table$bam_aligned_genomic[i] <- bam_aligned_path
  }
  print("--- Finished Loop 1: STAR Alignment & Indexing ---")


  
  table$bam_final_genomic <- NA 

  if (UMI == TRUE && params$UMI_dedup == TRUE) { 
    cat("\n--- Starting Parallel UMI Deduplication ---\n")

    dedup_input_dir <- bam_parent_dir
    dedup_output_dir <- paste0("./Output_", current_time, "/BAM_dedup")
    dir.create(dedup_output_dir, showWarnings = FALSE, recursive = TRUE)

    umi_script <- "/opt/UMI_parallel/para_umi_dedup.sh"

    if (!file.exists(umi_script)) stop("UMI parallel script not found at expected path: ", umi_script)

    system(paste("chmod +x", shQuote(umi_script)))
    if (!exists("umi_config") || !file.exists(umi_config)) stop("UMI config file path 'umi_config' is not defined or file does not exist.")
    if (!dir.exists(dedup_input_dir)) stop("Input directory for parallel deduplication not found: ", dedup_input_dir)

    cmd <- paste(
        "bash", shQuote(umi_script),
        "-i", shQuote(dedup_input_dir),
        "-o", shQuote(dedup_output_dir),
        "-t", CPU,
        "-f", shQuote(umi_config)
    )
    cat("Running parallel deduplication command:", cmd, "\n")
    status_dedup <- system(cmd)

    if(status_dedup != 0){
        stop("Error executing para_umi_dedup.sh. Check logs or script output.")
    } else {
        cat("Parallel UMI deduplication finished.\n")
        for (i in 1:nrow(table)) {
            sample_ID <- table$sample_ID[i]
            dedup_bam_path <- file.path(dedup_output_dir, paste0(sample_ID, "Aligned.sortedByCoord.out.dedup.bam"))
            if (!file.exists(dedup_bam_path)) {
                warning("Expected deduplicated BAM not found for sample ", sample_ID, ": ", dedup_bam_path)
                 table$bam_final_genomic[i] <- NA
            } else {
                table$bam_final_genomic[i] <- dedup_bam_path
            }
        }
        cat("Indexing deduplicated BAM files...\n")

        for (bam_path in na.omit(table$bam_final_genomic)) {
            cat("  Indexing:", basename(bam_path), "\n")
            status_index_dedup <- system2("samtools", paste("index", shQuote(bam_path)))
            if(status_index_dedup != 0) warning("Samtools index failed for deduplicated BAM: ", bam_path)
        }
        cat("Finished indexing deduplicated BAM files.\n")
    }
  } else {
    print("\nUMI deduplication skipped (UMI is FALSE or params$UMI_dedup is FALSE).")
    table$bam_final_genomic <- table$bam_aligned_genomic
  }


  # --- Loop 2: featureCounts ---
  
  print("\n--- Starting Loop 2: featureCounts ---")
  for (i in 1:nrow(table)) {
    sample_ID <- table$sample_ID[i]
    bam_file_for_counting <- table$bam_final_genomic[i]
    if (is.na(bam_file_for_counting) || !file.exists(bam_file_for_counting)) {
        warning("Skipping featureCounts for ", sample_ID, ": Final BAM path not found or file missing.")
        next
    }
    cat("  Counting features for sample:", sample_ID, "(", i, "/", nrow(table), ")\n")
    counts_dir <- file.path(counts_parent_dir, sample_ID)
    log_dir <- file.path(logs_parent_dir, sample_ID)
    fc_out_file <- file.path(counts_dir, paste0(sample_ID, "_counts.txt"))
    fc_log_file <- file.path(log_dir, paste0(sample_ID, "_featureCounts.log"))
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
    if(status_fc != 0) {
        warning("featureCounts failed for sample: ", sample_ID, ". Check log: ", fc_log_file)
    } else {
        print(paste("  Finished featureCounts for", sample_ID))
    }
  }
  

  print("finishing alignment stage")

} 





#RSEM




if(RSEM == TRUE){

  # --- RSEM Reference Handling ---
  if(transcriptome_indexing == TRUE) {
    
    print("Transcriptome indexing is TRUE. Starting RSEM reference creation...")

    spladder_out <- "/archive/home/paolo.zarcone/Bulk-RNA/IsoAnalyzer/Complete_RUN/Nodo1_rMATS_SplAdder_allignment/Output_251106/spladder_out"
    
    # --- GTF selection ---
    if(GTF_addition == TRUE){
      gtf_file_path <- file.path(spladder_out, "merged_sorted.gtf")
    } else {
      gtf_file_path <- gtf
    }
    
   
    gtf_file_abs <- tools::file_path_as_absolute(gtf_file_path)
    fasta_abs <- tools::file_path_as_absolute(fasta)

    
    rsem_ref_dir <- "./rsem_reference/GRCh38_RSEM" 
    dir.create(rsem_ref_dir, recursive = TRUE, showWarnings = FALSE)
    

    rsem_prefix_basename <- basename(rsem_ref_dir)
    
    Star_bin <- Sys.which("STAR")
    if(nchar(Star_bin) == 0) stop("STAR not found in PATH")
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
    
    if (status != 0) {
      stop(paste("rsem-prepare-reference failed with exit code:", status))
    }
   

    
    rsem_prefix_full <- file.path(rsem_ref_dir, rsem_prefix_basename)
  
  } else {
    
    print("Transcriptome indexing is FALSE. Skipping reference creation.")
    rsem_ref_dir <- transcriptome_genome_dir 
    print(paste("Using pre-existing RSEM reference from:", rsem_ref_dir))
    
    
    rsem_prefix_basename <- basename(rsem_ref_dir)
    rsem_prefix_full <- file.path(rsem_ref_dir, rsem_prefix_basename)
    
  
    if(!file.exists(paste0(rsem_prefix_full, ".grp"))){
        warning(paste("Warning: RSEM prefix not found at", rsem_prefix_full))
        warning("Please ensure the prefix name matches the directory name.")
    }
  }


  # --- RSEM Transcriptome Alignment ---
  if (transcriptome_alignment == TRUE) {
  
    print("Transcriptome alignment is TRUE. Starting STAR alignment for RSEM...")

    # --- STAR alignment per sample ---
    for(i in 1:nrow(table)){
      sample_ID <- table$sample_ID[i]
  
      
      fastq_files <- get_fastq_files(sample_ID, mode, out_dir, trimming, UMI)
  
      # --- If trimming=FALSE, create temporary trimmed files ---
      if(trimming == FALSE){
        tmp_out_dir <- paste0("./Output_", current_time, "/tmp_trimmed")
        dir.create(tmp_out_dir, recursive = TRUE, showWarnings = FALSE)
  
        if(mode == "single"){
          fastq <- table$fastq[i]
          if(UMI){
            tmp_fastq <- paste0(tmp_out_dir, "/", sample_ID, "_umi.fq.gz")
            system2("fastp", paste0(
              "-i ", fastq, " -o ", tmp_fastq, " --umi --umi_loc ", umi_loc,
              " --umi_len ", umi_len, " --dedup --umi_prefix ", umi_prefix,
              " --thread ", CPU,
              " --html ", tmp_out_dir, "/", sample_ID, "_fastp_report.html",
              " --json ", tmp_out_dir, "/", sample_ID, "_fastp_report.json"
            ))
          } else {
            tmp_fastq <- paste0(tmp_out_dir, "/", sample_ID, "_trimmed.fq.gz")
            system2("fastp", paste0(
              "-i ", fastq, " -o ", tmp_fastq,
              " --detect_adapter_for_pe ",
              ifelse(clip_5_1==0,"",paste0("--trim_front1 ",clip_5_1," ")),
              ifelse(clip_3_1==0,"",paste0("--trim_tail1 ",clip_3_1," ")),
              " --thread ", CPU,
              " --html ", tmp_out_dir, "/", sample_ID, "_fastp_report.html",
              " --json ", tmp_out_dir, "/", sample_ID, "_fastp_report.json"
            ))
          }
          fastq_files <- tmp_fastq
        } else { # paired-end
          fastq1 <- table$fastq1[i]
          fastq2 <- table$fastq2[i]
          if(UMI){
            tmp_fastq1 <- paste0(tmp_out_dir, "/", sample_ID, "_umi_R1.fq.gz")
            tmp_fastq2 <- paste0(tmp_out_dir, "/", sample_ID, "_umi_R2.fq.gz")
            system2("fastp", paste0(
              "-i ", fastq1, " -I ", fastq2,
              " -o ", tmp_fastq1, " -O ", tmp_fastq2,
              " --umi --dedup --umi_loc ", umi_loc, " --umi_len ", umi_len, " --umi_prefix ", umi_prefix,
              " --thread ", CPU,
              " --html ", tmp_out_dir, "/", sample_ID, "_fastp_report.html",
              " --json ", tmp_out_dir, "/", sample_ID, "_fastp_report.json"
            ))
          } else {
            tmp_fastq1 <- paste0(tmp_out_dir, "/", sample_ID, "_val_1.fq.gz")
            tmp_fastq2 <- paste0(tmp_out_dir, "/", sample_ID, "_val_2.fq.gz")
            system2("fastp", paste0(
              "-i ", fastq1, " -I ", fastq2,
              " -o ", tmp_fastq1, " -O ", tmp_fastq2,
              " --detect_adapter_for_pe ",
              ifelse(clip_5_1==0,"",paste0("--trim_front1 ",clip_5_1," ")),
              ifelse(clip_3_1==0,"",paste0("--trim_tail1 ",clip_3_1," ")),
              ifelse(clip_5_2==0,"",paste0("--trim_front2 ",clip_5_2," ")),
              ifelse(clip_3_2==0,"",paste0("--trim_tail2 ",clip_3_2," ")),
              " --thread ", CPU,
              " --html ", tmp_out_dir, "/", sample_ID, "_fastp_report.html",
              " --json ", tmp_out_dir, "/", sample_ID, "_fastp_report.json"
            ))
          }
          fastq_files <- list(tmp_fastq1, tmp_fastq2)
        }
      }
  
      # --- Output directories ---
      bam_dir <- paste0("./Output_", current_time, "/Transcriptome/BAM/", sample_ID, "/")
      log_dir <- paste0("./Output_", current_time, "/Transcriptome/logs/", sample_ID, "/")
      counts_dir <- paste0("./Output_", current_time, "/Transcriptome/counts/", sample_ID, "/")
      dir.create(bam_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
  
      star_prefix <- paste0(bam_dir, sample_ID)
      
      needs_zcat <- function(file) {
      if (is.null(file)) return(FALSE)
      endsWith(file, ".gz")
      }
  
      # --- Run STAR alignment ---
      
      if(mode == "single"){
        
        zcat_flag <- if (needs_zcat(fastq_files)) " --readFilesCommand zcat" else ""

        system2(Star_bin, args = c(
          "--runThreadN", CPU,
          "--genomeDir", rsem_ref_dir,
          "--readFilesIn", fastq_files,
          "--outSAMtype", "None",
          "--quantMode TranscriptomeSAM", zcat_flag,
          "--outFileNamePrefix", star_prefix
        ))
      } else {

        zcat_flag <- if (needs_zcat(fastq_files[[1]])) " --readFilesCommand zcat" else ""

        system2(Star_bin, args = c(
          "--runThreadN", CPU,
          "--genomeDir", rsem_ref_dir,
          "--readFilesIn", fastq_files[[1]], fastq_files[[2]],
          "--outSAMtype", "None",
          "--quantMode TranscriptomeSAM", zcat_flag,
          "--outFileNamePrefix", star_prefix
        ))
      }

      bam_file_per_rsem <- paste0(star_prefix, "Aligned.toTranscriptome.out.bam")
  
      table$bam_final[i] <- bam_file_per_rsem 
    }
    } else {
      
    print("Transcriptome alignment is FALSE. Skipping alignment loop.")
    print(paste("Locating pre-existing BAM files in:", bam_dir_transcriptome))
    
    for(i in 1:nrow(table)){
      sample_ID <- table$sample_ID[i]
      
      
      bam_file <- file.path(bam_dir_transcriptome, paste0(sample_ID, "Aligned.toTranscriptome.out.bam"))

      if(!file.exists(bam_file)){
        stop(paste("Could not find expected *Transcriptome* BAM file:", bam_file))
      }
      
      table$bam_final[i] <- bam_file
    }
    print("Transcriptome BAM file paths successfully located.")
  }
  

}


print("Pipeline rmats + Spaldder + Rsem allignment finished.")

