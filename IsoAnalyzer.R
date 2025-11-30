#------------------------------------------------------------------------------
#Isoform Analyzer
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

set.seed(1234)
current_time <- format(Sys.time(), "%y%m%d")
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
  dir.create(paste0("./Output_", current_time, "/trimmed"), showWarnings = FALSE)
  out_dir = paste0("./Output_", current_time, "/trimmed")
  
  if(mode == "single"){
    
    for(i in 1:length(rownames(table))){
      sample_ID = table$sample_ID[i]
      fastq = table$fastq[i]
      
      if (UMI == TRUE){
        # -- Single End + UMI (Fastp) --
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
        
        if(fastqc == TRUE){
          system2("fastqc", paste0(cleaned_fastq, " -o ", out_dir, " -t ", CPU))
        }
        
        fastq = cleaned_fastq
        
      } else {
        if(is.na(adapter) || adapter == ""){
          # -- Single End + Auto Adapter (Fastp) --
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
          
          if(fastqc == TRUE){
            system2("fastqc", paste0(cleaned_fastq, " -o ", out_dir, " -t ", CPU))
          }
          
          fastq = cleaned_fastq
          
        } else {
          # -- Single End + Adapter (Trim Galore) --
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
      table$fastq[i] <- fastq
    }
    
  } else if(mode == "paired"){ 
    
    for(i in 1:length(rownames(table))){
      sample_ID = table$sample_ID[i]
      fastq1 = table$fastq1[i]
      fastq2 = table$fastq2[i]
      
      if (UMI == TRUE){
        # -- Paired End + UMI (Fastp) --
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
        
        if(fastqc == TRUE){
          system2("fastqc", paste0(cleaned_fastq1, " ", cleaned_fastq2, " -o ", out_dir, " -t ", CPU))
        }
        
        fastq1 = cleaned_fastq1
        fastq2 = cleaned_fastq2
        
      } else {
        if((is.na(adapter) || adapter == "") && (is.na(adapter2) || adapter2 == "")){
          # -- Paired End + Auto Adapter (Fastp) --
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
          
          if(fastqc == TRUE){
            system2("fastqc", paste0(cleaned_fastq1, " ", cleaned_fastq2, " -o ", out_dir, " -t ", CPU))
          }
          
          fastq1 = cleaned_fastq1
          fastq2 = cleaned_fastq2
          
        } else {
          # -- Paired End + Adapter (Trim Galore) --
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

      table$fastq1[i] <- fastq1
      table$fastq2[i] <- fastq2
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



#Rmats
if (Rmats == TRUE && alignment == TRUE) {
  print("Preparing rMATS input files")

  
  out_dir <- paste0("./Output_", current_time, "/rmatsout")
  out_OD  <- file.path(out_dir, "OD")
  out_TM  <- file.path(out_dir, "TM")
  out_txt  <- file.path(out_dir, "txt")
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_OD, showWarnings = FALSE)
  dir.create(out_TM, showWarnings = FALSE)
  dir.create(out_txt, showWarnings = FALSE)

  # sample_ID and conditions 
  conditions <- as.character(table[[6]])
  sample_IDs <- as.character(table$sample_ID)

  # creation of B1.txt and B2.txt
  bam_list <- list("1" = character(), "2" = character())


  for (i in seq_along(sample_IDs)) {
    sample_ID <- sample_IDs[i]
    bam_dir   <- paste0("./Output_", current_time, "/BAM/", sample_ID)


    # different suffix if UMI deduplication was performed
    

    suffix <- if (params$UMI_dedup == TRUE) "dedup.bam" else "Aligned.sortedByCoord.out.bam"

    bam_file <- file.path(bam_dir, paste0(sample_ID, suffix))



    # Bam check
    if (file.exists(bam_file)) {
      # normalizzo il gruppo
      group_val <- suppressWarnings(as.integer(as.character(conditions[i])))
      if (is.na(group_val)) {
        warning(sprintf("Invalid/NA group for sample %s (groups[i] = %s). Skipping.", sample_ID, conditions[i]))
        next
      }
      if (group_val == 1L) {
        bam_list[["1"]] <- c(bam_list[["1"]], bam_file)
      } else if (group_val == 2L) {
        bam_list[["2"]] <- c(bam_list[["2"]], bam_file)
      } else {
        warning(sprintf("Group %s for sample %s is not 1 or 2. Skipping.", group_val, sample_ID))
      }
    } else {
      warning(sprintf("No BAM found for sample %s. Expected: %s", sample_ID, bam_file))
    }
  }

  # remove duplicates
  bam_list[["1"]] <- unique(bam_list[["1"]])
  bam_list[["2"]] <- unique(bam_list[["2"]])

  message("n B1: ", length(bam_list[["1"]]))
  message("n B2: ", length(bam_list[["2"]]))

  
  if (length(bam_list[["1"]]) > 0) {
    write(paste(bam_list[["1"]], collapse=","), file = file.path(out_txt, "B1.txt"))
    message("Wrote B1.txt")
  }
  if (length(bam_list[["2"]]) > 0) {
    write(paste(bam_list[["2"]], collapse=","), file = file.path(out_txt, "B2.txt"))
    message("Wrote B2.txt")
  }

  print(warnings())
  print("Finished preparing rMATS input files (B1.txt and B2.txt)")
}


# rMATS execution



if (Rmats == TRUE && alignment == TRUE) {
  print("starting rMATS")

  b1_path <- file.path(out_txt, "B1.txt")
  b2_path <- file.path(out_txt, "B2.txt")


  system2("python3", args = c(
    "/rmats/rmats.py",
    "--b1", b1_path,
    "--b2", b2_path,
    "--gtf", gtf,
    "-t", "paired",
    "--readLength", "151",
    "--variable-read-length",
    "--nthread", as.character(CPU),
    "--od", out_OD,
    "--tmp", out_TM
  ))
}

if (Rmats == TRUE && alignment == FALSE ) {
  print("starting rMATS")

  out_dir <- paste0("./Output_", current_time, "/rmatsout")
  out_OD  <- file.path(out_dir, "OD")
  out_TM  <- file.path(out_dir, "TM")

  dir.create(out_dir, recursive = TRUE)
  dir.create(out_OD)
  dir.create(out_TM)

  system2("python3", args = c(
    "/rmats/rmats.py",
    "--b1", file.path(rmats_path, "b1.txt"),
    "--b2", file.path(rmats_path, "b2.txt"),
    "--gtf", gtf,
    "-t", "paired",
    "--readLength", "151",
    "--variable-read-length",
    "--nthread", as.character(CPU),
    "--od", out_OD,
    "--tmp", out_TM
  ))
}







# SplAdder




# SplAdder input files preparation

if (Spladder == TRUE) {
  
  if (alignment == FALSE) {
    # --- Case 1: pre-existing BAMs and txt files ---
    if (length(unique(na.omit(settings$txt_path))) != 1 || !file.exists(settings$txt_path[1])) {
      stop("Invalid argument for txt_path")
    } else {
      txt_path <- as.character(settings$txt_path[1])
    }
    
    if (length(unique(na.omit(settings$bam_dir))) != 1 || !dir.exists(settings$bam_dir[1])) {
      stop("Invalid argument for bam_dir")
    } else {
      bam_dir <- as.character(settings$bam_dir[1])
    }

    all_bams_file   <- file.path(txt_path, "all_bams.txt")
    conditionA_file <- file.path(txt_path, "conditionA.txt")
    conditionB_file <- file.path(txt_path, "conditionB.txt")

    if (!file.exists(all_bams_file))   stop("all_bams.txt not found in txt_path")
    if (!file.exists(conditionA_file)) stop("conditionA.txt not found in txt_path")
    if (!file.exists(conditionB_file)) stop("conditionB.txt not found in txt_path")

    bams_A <- readLines(conditionA_file)
    bams_B <- readLines(conditionB_file)
    all_bams <- readLines(all_bams_file)

    message("Using pre-existing SplAdder input files from: ", txt_path)
    message("BAM directory (from settings): ", bam_dir)
    message("n condA: ", length(bams_A), " | n condB: ", length(bams_B))

    # --- BAM indexing ---
    bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
    for (bam in bam_files) {
      bai <- paste0(bam, ".bai")
      if (!file.exists(bai)) {
        message("Indexing missing BAM: ", bam)
        system2("samtools", c("index", bam))
      } else {
        message("Index already present: ", bam)
      }
    }
    print("Finished indexing all BAMs in directory")

  } else if (alignment == TRUE) {
    # --- Case 2: freshly aligned BAMs ---
    print("Preparing SplAdder input files from freshly aligned BAMs")
    
    out_dir <- paste0("./Output_", current_time, "/spladder")
    out_txt <- file.path(out_dir, "txt")

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_txt, showWarnings = FALSE)

    # 1) BAM list for build function
    all_bams <- unique(c(bam_list[["1"]], bam_list[["2"]]))
    if (length(all_bams) > 0) {
      writeLines(all_bams, file.path(out_txt, "all_bams.txt"))
      message("Wrote all_bams.txt for SplAdder build")
    } else {
      warning("No BAMs found to write all_bams.txt")
    }

    # 2) Condition A file
    if (length(bam_list[["1"]]) > 0) {
      writeLines(bam_list[["1"]], file.path(out_txt, "conditionA.txt"))
      message("Wrote conditionA.txt")
    } else {
      warning("No BAMs found for condition A")
    }

    # 3) Condition B file
    if (length(bam_list[["2"]]) > 0) {
      writeLines(bam_list[["2"]], file.path(out_txt, "conditionB.txt"))
      message("Wrote conditionB.txt")
    } else {
      warning("No BAMs found for condition B")
    }

    
    all_bams_file <- file.path(out_txt, "all_bams.txt")
    all_bams <- readLines(all_bams_file)
    conditionA_file <- file.path(out_txt, "conditionA.txt")
    conditionB_file <- file.path(out_txt, "conditionB.txt")

    

    print("Finished preparing SplAdder input files (all_bams.txt, conditionA.txt, conditionB.txt)")
  }
}



library(parallel)


spladder_bin <- Sys.which("spladder")
if (nchar(spladder_bin) == 0) stop("spladder non trovato nel PATH")


spladder_out <- paste0("./Output_", current_time, "/spladder_out")
dir.create(spladder_out, recursive = TRUE, showWarnings = FALSE)
outdir <- spladder_out

mc.cores <- min(4, ceiling(detectCores() / 4))



if (Spladder == TRUE ) {
  print("Starting SplAdder build")

  library(parallel)

  spladder_bin <- Sys.which("spladder")
  if (nchar(spladder_bin) == 0) stop("spladder non trovato nel PATH")


  spladder_out <- paste0("./Output_", current_time, "/spladder_out")
  dir.create(spladder_out, recursive = TRUE, showWarnings = FALSE)

  outdir <- spladder_out

  mc.cores <- min(4, ceiling(detectCores() / 4))



  # Define filtered GTF path
  spladder_dir <- paste0("./Output_", current_time, "/spladder")
  dir.create(spladder_dir, recursive = TRUE, showWarnings = FALSE)

  gtf_pc <- file.path(spladder_dir, "Homo_sapiens.GRCh38.113.prot.gtf")

  # Filter only if not already existing
  if (!file.exists(gtf_pc)) {
    message("Creating filtered GTF (protein_coding only)...")
    
    cmd <- paste("grep 'protein_coding' ", shQuote(gtf), " > ", shQuote(gtf_pc), sep = "")
    system(cmd)
    
    # check that the file was actually created
    if (!file.exists(gtf_pc) || file.info(gtf_pc)$size == 0) {
      stop("Filtered GTF file was not created correctly. Check path or grep command.")
    }
    
  } else {
    message("Filtered GTF already exists, skipping creation.")
  }

  message("Filtered GTF created at: ", gtf_pc)



  # ----- Define log files -----
  log_step1 <- file.path(outdir, "spladder_step1_single_graphs.log")
  log_step2 <- file.path(outdir, "spladder_step2_merge_graphs.log")
  log_step3a <- file.path(outdir, "spladder_step3a_quant_single.log")
  log_step3b <- file.path(outdir, "spladder_step3b_quant_collect.log")
  log_step4 <- file.path(outdir, "spladder_step4_event_extraction.log")

  # ----- 1) Single graphs for BAM -----
  build_single_graph <- function(bam) {
    args <- c("build",
              "--outdir", outdir,
              "--annotation", gtf_pc,
              "--bams", bam,
              "--merge-strat", "single",
              "--no-extract-ase",
              "--set-mm-tag", "nM")
    
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Building single graph for: ", bam)
    system2(spladder_bin, args=args, stdout=log_step1, stderr=log_step1)
    return(bam)
  }

  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Step 1: Building single graphs in parallel")
  mclapply(all_bams, build_single_graph, mc.cores = mc.cores)


  # ----- 2) Merge graphs -----
  merge_args <- c("build",
                  "--outdir", outdir,
                  "--annotation", gtf_pc,
                  "--bams", all_bams_file,
                  "--merge-strat", "merge_graphs",
                  "--no-extract-ase",
                  "--set-mm-tag", "nM")

  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Step 2: Merging graphs into a single merged graph")
  system2(spladder_bin, args=merge_args, stdout=log_step2, stderr=log_step2)



  # ----- 3a) Quantification per BAM -----
  quant_single <- function(bam) {
    args <- c("build",
              "--outdir", outdir,
              "--annotation", gtf_pc,
              "--bams", bam,
              "--merge-strat", "merge_graphs",
              "--no-extract-ase",
              "--quantify-graph",
              "--qmode", "single",
              "--set-mm-tag", "nM")
    
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Quantifying merged graph for: ", bam)
    system2(spladder_bin, args=args, stdout=log_step3a, stderr=log_step3a)
    return(bam)
  }

  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Step 3a: Quantification per BAM in parallel (qmode single)")
  mclapply(all_bams, quant_single, mc.cores = mc.cores)


  # ----- 3b) Quantification collection -----
  collect_args <- c("build",
                    "--outdir", outdir,
                    "--annotation", gtf_pc,
                    "--bams", all_bams_file,
                    "--merge-strat", "merge_graphs",
                    "--no-extract-ase",
                    "--quantify-graph",
                    "--qmode", "collect",
                    "--set-mm-tag", "nM")

  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Step 3b: Collecting all quantifications (qmode collect)")
  system2(spladder_bin, args=collect_args, stdout=log_step3b, stderr=log_step3b)


  # ----- 4) Final event extraction -----
  event_args <- c("build",
                  "--outdir", outdir,
                  "--annotation", gtf_pc,
                  "--bams", all_bams_file,
                  "--merge-strat", "merge_graphs",
                  "--event-types", "exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons",
                  "--ase-edge-limit", "1000",
                  "--set-mm-tag", "nM")

  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | Step 4: Extracting ASE events on the merged graph")
  system2(spladder_bin, args=event_args, stdout=log_step4, stderr=log_step4)


  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          " | Pipeline completed: single graphs -> merge -> quantification -> event extraction")


}

#SplAdder test



if (Spladder == TRUE ) {

  outdir <- spladder_out
  print("Starting SplAdder test")
  test_args <- c(
    "test",
    "--conditionA", conditionA_file,
    "--conditionB", conditionB_file,
    "--labelA", "treated",
    "--labelB", "Wildtype",
    "--parallel", CPU,
    "--diagnose-plots",
    "--plot-format", "pdf",
    "--outdir", outdir
  )


  message("Running SplAdder test step")
  system2(spladder_bin, args=test_args, stdout=TRUE, stderr=TRUE)
}

#Spladder viz


if (Spladder == TRUE ) {

  message("Starting SplAdder viz (parallel mode)")

  test_dir <- file.path(spladder_out, "testing_treated_vs_Wildtype")
  testcase <- "testing_treated_vs_Wildtype"


  files_A <- readLines(conditionA_file)
  files_A <- files_A[files_A != ""]
  stringa_A <- paste0("ConditionA:", paste(files_A, collapse = ","))

  files_B <- readLines(conditionB_file)
  files_B <- files_B[files_B != ""]
  stringa_B <- paste0("ConditionB:", paste(files_B, collapse = ","))




  # ---- events type ----
  event_types <- c("exon_skip", "intron_retention", "alt_3prime",
                   "alt_5prime", "mutex_exons", "mult_exon_skip")
 
  plot_top_events <- function(event_type, top_n = 10) {
    args <- c("viz",
              "--test", testcase, event_type, as.character(top_n),
              "--testdir", test_dir,
              "--outdir", spladder_out,
              "--confidence", "3",
              "--track", "coverage", stringa_A, stringa_B,
              "--track", "splicegraph",          
              "--outbase", paste0("top", top_n, "_"),
              "--format", "pdf")
  

    message("Plotting top ", top_n, " events for: ", event_type)
    system2(spladder_bin, args = args, stdout = TRUE, stderr = TRUE)
  }


  plot_top_events_log <- function(event_type, top_n = 10) {
    args <- c("viz",
              "--test", testcase, event_type, as.character(top_n),
              "--testdir", test_dir,
              "--outdir", spladder_out,
              "--confidence", "3",
              "--log",
              "--track", "coverage", stringa_A, stringa_B,
              "--track", "splicegraph",          
              "--outbase", paste0("top", top_n, "_log_scale_"),
              "--format", "pdf")
  

    message("Plotting logarithmic version of top ", top_n, " events for: ", event_type)
    system2(spladder_bin, args = args, stdout = TRUE, stderr = TRUE)
  }


  n_jobs <- length(event_types)

  mclapply(event_types, plot_top_events, mc.cores = n_jobs)
  mclapply(event_types, plot_top_events_log, mc.cores = n_jobs)

  message("All top events plotted in ", file.path(test_dir, "plots"))
}




#Updated GTF

if (GTF_addition == TRUE ) {
  
  gtf_original <- as.character(settings$gtf[1])

  if (!file.exists(gtf_original)) stop("Original GTF file not found.")
  if (!dir.exists(spladder_out)) stop("SplAdder output directory not found.")
  
  # --- GFF3 produced by spladder for each category ---
  gff3_files <- list.files(
    path = spladder_out,
    pattern = "merge_graphs_.*_C[0-9]+\\.confirmed\\.gff3$",
    full.names = TRUE
  )
  
  if (length(gff3_files) == 0) stop("No merge_graphs_*.confirmed.gff3 files found.")
  
  # --- GFF3 merge ---
  spladder_combined_gff <- file.path(spladder_out, "spladder_combined.gff3")
  system(paste("cat", paste(shQuote(gff3_files), collapse = " "), ">", shQuote(spladder_combined_gff)))
  
  # --- GFF3 to GTF ---
  gff3_to_gtf <- Sys.which("rsem-gff3-to-gtf")
  if (nchar(gff3_to_gtf) == 0) stop("rsem-gff3-to-gtf not found in PATH")

  spladder_gtf <- file.path(spladder_out, "spladder_combined.gtf") # "Dirty" file
  system2("python3", args = c(gff3_to_gtf, spladder_combined_gff, spladder_gtf, "--RNA-patterns", "mRNA"))
  
  
  # --- New GTF filtering ---
  
  # temporary files
  spladder_cleaned_gtf <- file.path(spladder_out, "spladder_combined.cleaned.gtf")
  spladder_filtered_gtf <- file.path(spladder_out, "spladder_combined.filtered.gtf")
  spladder_prefixed_gtf <- file.path(spladder_out, "spladder_prefixed.gtf")
  merged_temp <- file.path(spladder_out, "merged_temp.gtf")
  merged_sorted <- file.path(spladder_out, "merged_sorted.gtf") 


  # Remove non-printable characters (\u0001)
  cat("Step 1/5: Cleaning non-printable characters (Perl)...\n")
  clean_cmd <- sprintf("perl -pe 's/[^[:print:]\\t\\n]//g' %s > %s", 
                       shQuote(spladder_gtf), 
                       shQuote(spladder_cleaned_gtf))
  system(clean_cmd)
  if (!file.exists(spladder_cleaned_gtf)) stop("Step 1 (Perl clean) failed.")

  
  # STEP 2: Filter (Awk) 
  cat("Step 2/5: Filtering empty IDs (Awk)...\n")
  filter_cmd <- sprintf("awk -F'\\t' 'BEGIN{OFS=\"\\t\"} $9 ~ /transcript_id \"\"/ || $9 ~ /gene_id \"\"/ {next} {print}' %s > %s", 
                        shQuote(spladder_cleaned_gtf), 
                        shQuote(spladder_filtered_gtf))
  system(filter_cmd)
  if (!file.exists(spladder_filtered_gtf)) stop("Step 2 (Awk filter) failed.")

  
  # STEP 3: Prefix (Awk) - Add 'spladder_' 
  cat("Step 3/5: Adding 'spladder_' prefix (Awk)...\n")
  awk_cmd <- sprintf(
    "awk -F'\\t' 'BEGIN{OFS=\"\\t\"} { if($3==\"transcript\" || $3==\"exon\"){ gsub(/transcript_id \"/, \"transcript_id \\\"spladder_\"); gsub(/gene_id \"/, \"gene_id \\\"spladder_\") } print }' %s > %s",
    shQuote(spladder_filtered_gtf), 
    shQuote(spladder_prefixed_gtf)
  )
  system(awk_cmd)
  if (!file.exists(spladder_prefixed_gtf)) stop("Step 3 (Awk prefix) failed.")

  
  # STEP 4: Merge (Cat) 
  cat("Step 4/5: Merging with original GTF...\n")
  system(paste("cat", shQuote(gtf_original), shQuote(spladder_prefixed_gtf), ">", shQuote(merged_temp)))

  
  # STEP 5: Sort (Sort) - Sort the final file
  cat("Step 5/5: Sorting final file...\n")
  sort_cmd <- sprintf("sort -k1,1 -k4,4n %s > %s", 
                      shQuote(merged_temp), 
                      shQuote(merged_sorted))
  system(sort_cmd)
  
  
  # Clean up temporary files
  file.remove(spladder_gtf) 
  file.remove(spladder_cleaned_gtf)
  file.remove(spladder_filtered_gtf)
  file.remove(spladder_prefixed_gtf)
  file.remove(merged_temp)

  cat("New GTF created successfully:", merged_sorted, "\n")
}



#RSEM




if(RSEM == TRUE){

  # --- RSEM Reference Handling ---
  if(transcriptome_indexing == TRUE) {
    
    print("Transcriptome indexing is TRUE. Starting RSEM reference creation...")
    
    # --- GTF selection ---
    if(GTF_addition == TRUE){
      gtf_file_path <- file.path(spladder_out, "merged_sorted.gtf")
    } else {
      gtf_file_path <- gtf
    }
    
   
    gtf_file_abs <- tools::file_path_as_absolute(gtf_file_path)
    fasta_abs <- tools::file_path_as_absolute(fasta)

    
    rsem_ref_dir <- paste0("./Output_", current_time, "rsem_reference/GRCh38_RSEM") 
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
  

  # --- RSEM expression quantification ---
  for(i in 1:nrow(table)){
    sample_ID <- table$sample_ID[i]
    bam_file <- table$bam_final[i]
    counts_dir <- paste0("./Output_", current_time, "/Transcriptome/counts/", sample_ID, "/")
    dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)

    
    rsem_args <- c(
      "--alignments",      
      "--no-bam-output",
      "--num-threads", CPU,
      "--calc-pme",
      "--seed", 12345,
      bam_file,            
      rsem_prefix_full, 
      paste0(counts_dir, sample_ID)
    )

    if(mode == "paired") {
      rsem_args <- c("--paired-end", rsem_args)
    }

    system2("rsem-calculate-expression", args = rsem_args)
  }

  print("RSEM expression calculation finished for all samples")
}






#Data import

if(RSEM == TRUE) {

  library(IsoformSwitchAnalyzeR)
 

  cat("\n--- Preparing data for IsoformSwitchAnalyzeR ---\n")

  # --- 1. Define RSEM base directory ---
  base_dir <- paste0("./Output_", current_time) 
  rsem_base_dir <- file.path(base_dir, "Transcriptome", "counts")

  if(!dir.exists(rsem_base_dir)){
    warning("RSEM base directory not found: ", rsem_base_dir, ". Skipping IsoformSwitchAnalyzeR import.")
  } else {
    cat("Using RSEM base directory:", rsem_base_dir, "\n")

    
    # --- 2. Design matrix building ---
    cat("Building design matrix...\n")
    myDesign <- data.frame(
      sampleID  = table$sample_ID, 
      condition = table$condition  
    )

    if("replicate" %in% colnames(table)) {
      if(length(unique(na.omit(table$replicate))) > 1) {
        myDesign$replicate <- table$replicate
        cat("  Replicate information included (column is not constant).\n")
      } else {
        cat("  'replicate' column found but ignored (all values are identical).\n")
      }
    }

    if(GTF_addition == TRUE){ 
          spladder_out_dir <- file.path(base_dir, "spladder_out") 
          
        
          gtf_source <- file.path(spladder_out_dir, "merged_sorted.gtf")
          gff3_source <- file.path(spladder_out_dir, "spladder_combined.gff3") 
          
          cat("Using potentially merged GTF from SplAdder:", gtf_source, "\n")
          

          gene_info_file <- NA  
          if(transcriptome_indexing == TRUE) {
              gene_info_file <- file.path(base_dir, "rsem_reference/GRCh38_RSEM/geneInfo.tab")
          } else {
              gene_info_file <- file.path(transcriptome_genome_dir, "geneInfo.tab")
          }

          cat( "using: " , gene_info_file, "\n" , gtf_source , "\n" , gff3_source , "\n" )
          
          if (!is.na(gene_info_file) && file.exists(gene_info_file) && file.exists(gff3_source)) {
              
              gtf_obj <- rtracklayer::import(gtf_source)
            
              gene_map <- readr::read_tsv(gene_info_file, col_names = FALSE, show_col_types = FALSE, skip = 1) %>%
                  dplyr::select(gene_id = 1, gene_symbol = 2)
              
              cat("  Parsing GFF3 to link gene name to splicing events...\n")
              gff_text <- readr::read_tsv(gff3_source, comment = "#", col_names = FALSE, show_col_types = FALSE)
              
              spladder_map <- gff_text %>%
                  dplyr::select(X9) %>%
                  dplyr::filter(stringr::str_detect(X9, "ID=") & stringr::str_detect(X9, "GeneName=")) %>%
                  dplyr::mutate(
                      # 1. Id extraction from GFF3 (es. "alt_3prime.123")
                      raw_id = stringr::str_extract(X9, "ID=([^;]+)"),
                      raw_id = stringr::str_remove(raw_id, "ID="),
                      
                      
                      transcript_id = paste0("spladder_", raw_id),
                      
                    
                      gene_id_raw = stringr::str_extract(X9, 'GeneName="([^"]+)"'),
                      gene_id = stringr::str_remove_all(gene_id_raw, 'GeneName=|"')
                  ) %>%
                  dplyr::select(transcript_id, gene_id) %>%
                  dplyr::distinct()
              
              # Maps Merging
              complete_map <- spladder_map %>%
                  dplyr::left_join(gene_map, by = "gene_id") %>%
                  dplyr::filter(!is.na(gene_symbol)) %>%
                  dplyr::mutate(new_gene_name = paste0("spladder_", gene_symbol)) %>%
                  dplyr::select(transcript_id, new_gene_name)
              
              cat(paste("Found", nrow(complete_map), "gene name for SplAdder events.\n"))

              
              meta_df <- as.data.frame(S4Vectors::mcols(gtf_obj))
              if(!"gene_name" %in% colnames(meta_df)) meta_df$gene_name <- NA
              
              cat(" Performing name addition to the GTF...\n")
              meta_df <- meta_df %>%
                  dplyr::left_join(complete_map, by = "transcript_id") %>%
                  dplyr::mutate(
                      gene_name = ifelse(is.na(gene_name) | gene_name == "", new_gene_name, gene_name)
                  ) %>%
                  dplyr::select(-new_gene_name)
              
              S4Vectors::mcols(gtf_obj) <- meta_df
              
            
              gtf_enriched_path <- file.path(spladder_out_dir, "merged_sorted_enriched.gtf")
              cat("  Esporting new gtf:", gtf_enriched_path, "\n")
              rtracklayer::export(gtf_obj, gtf_enriched_path, format="GTF")
              
              gtf_per_isoformswitch <- gtf_enriched_path
              
          } else {
              warning("geneInfo.tab or spladder_combined.gff3 not founded, skipping addition.")
              gtf_per_isoformswitch <- gtf_source
          }

    } else {
        gtf_per_isoformswitch <- gtf 
        cat("Using original GTF:", gtf_per_isoformswitch, "\n")
    }
      
    
    if(!file.exists(gtf_per_isoformswitch)) stop("GTF file for IsoformSwitchAnalyzeR not found:", gtf_per_isoformswitch)

    local_rsem_ref_dir <- NA_character_ 
    if (transcriptome_indexing == FALSE) { 
        local_rsem_ref_dir <- transcriptome_genome_dir 
    } else {
        local_rsem_ref_dir <- file.path(base_dir "/rsem_reference/GRCh38_RSEM") 
        warning("transcriptome_indexing=TRUE in filter-only script. Assuming default ref path: ", local_rsem_ref_dir)
    }
    if (is.na(local_rsem_ref_dir) || !dir.exists(local_rsem_ref_dir)) {
        stop("RSEM reference directory not found.")
    }
    local_rsem_prefix_basename <- basename(local_rsem_ref_dir)
    transcript_fasta_path <- file.path(local_rsem_ref_dir, paste0(local_rsem_prefix_basename, ".transcripts.fa"))
    cat("Using transcript FASTA:", transcript_fasta_path, "\n")
    if(!file.exists(transcript_fasta_path)) stop("Transcript FASTA file not found:", transcript_fasta_path)


    #Prepare GTF 
    
    out_gtf  <- file.path(base_dir, "GTF_cleaned") 
    dir.create(out_gtf, recursive = TRUE, showWarnings = FALSE) 

    python_script_path <- "/opt/isoformSwitchAnalyzeR/modify_gtf.py"
    gtf_cleaned_path <- file.path(out_gtf, "isoform_switch_cleaned.gtf")

    if(!file.exists(python_script_path)) stop("Script Python 'modify_gtf.py' not found.")

    cat("Running python script to clean GTF...\n")
    system2("python3", args = c(
      shQuote(python_script_path),
      shQuote(gtf_per_isoformswitch),
      shQuote(gtf_cleaned_path)
    ))
    if(!file.exists(gtf_cleaned_path)) stop("Cleaned GTF was not created by python script.")
    cat("GTF cleaned successfully.\n")


    #Load data into IsoformSwitchAnalyzeR

    cat("Importing RSEM quantification...\n")
    data <- importIsoformExpression(parentDir = rsem_base_dir)  
  
    cat("Creating switchAnalyzeRlist object...\n")
    switchList <- tryCatch({
        IsoformSwitchAnalyzeR::importRdata(
        isoformCountMatrix   = data$counts,  
        isoformRepExpression = data$abundance,  
        designMatrix         = myDesign,  
        isoformExonAnnoation = gtf_cleaned_path,
        isoformNtFasta       = transcript_fasta_path
        )
    }, error = function(e){
        warning("Error creating switchAnalyzeRlist: ", e$message, ". Skipping.")
        NULL 
    })

    # Proceed only if switchList was created
    if (!is.null(switchList)) {
      cat("\nswitchAnalyzeRlist created successfully.\n")
      cat("Number of isoforms imported:", nrow(switchList$isoformFeatures), "\n")
      cat("Number of genes imported:", length(unique(switchList$isoformFeatures$gene_id)), "\n")

      # --- 6. Save the result ---
      output_dir_isoformswitch <- file.path(base_dir, "IsoformSwitchAnalyzer_Output")
      dir.create(output_dir_isoformswitch, showWarnings = FALSE, recursive = TRUE)
      output_rds_file <- file.path(output_dir_isoformswitch, "switchAnalyzeRlist_imported.rds")

      cat("Saving switchAnalyzeRlist object to:", output_rds_file, "\n")
      saveRDS(switchList, file = output_rds_file)
    } 
    
  }  

  cat("--- Data import for IsoformSwitchAnalyzeR finished ---\n")

} else {
  cat("\nRSEM set to FALSE. Skipping IsoformSwitchAnalyzeR data import.\n")
}


#IsoformSwitchAnalyzeR Part1
cat("\n--- Starting IsoformSwitchAnalyzeR Part 1 ---\n")

library(BSgenome.Hsapiens.UCSC.hg38)
fasta_output_dir <- file.path(base_dir, "Analysis", "fasta")
dir.create(fasta_output_dir, recursive = TRUE, showWarnings = FALSE)

switchList <- readRDS(output_rds_file) 

switchList_part1 <- isoformSwitchAnalysisPart1(
    switchAnalyzeRlist = switchList,
    pathToGTF = gtf_cleaned_path,
    genomeObject = NULL, 
    outputSequences = TRUE, 
    pathToOutput = fasta_output_dir,
    prepareForWebServers = FALSE
 
)

output_rds_file_p1 <- file.path(output_dir_isoformswitch, "switchList_part1.rds")
cat("Saving switchList_part1 object to:",output_rds_file_p1, "\n")
saveRDS(switchList_part1, file = output_rds_file_p1)

#Gene table output filtering


if(Rmats == TRUE || Spladder == TRUE || RSEM == TRUE) {
    library(dplyr)
    library(readr)
    library(stringr)
    cat("\n--- Starting Output Filtering ---\n")
} else {
    cat("\n--- No filtering steps selected (Rmats, Spladder, RSEM are FALSE) ---\n")
}

# ==== Common Settings (defined only if needed) ====
if(Rmats == TRUE || SplAdder == TRUE || RSEM == TRUE) {

    # Common thresholds
    fdr_cutoff <- 0.05
    pval_cutoff <- 0.01 
    psi_diff_cutoff <- 0.05
    jc_min <- 5 
    tpm_cutoff <- 1
    count_cutoff <- 10
}

# -----------------------------------------------------
# 1. rMATS FILTERING
# -----------------------------------------------------
if(Rmats == TRUE) {
    cat("Processing rMATS output...\n")

    input_dir_rmats <- file.path(base_dir, "rmatsout", "OD")
    output_dir_rmats <- file.path(base_dir, "rmatsout", "Filtered_rMATS")
    dir.create(output_dir_rmats, showWarnings = FALSE, recursive = TRUE)

    if (!dir.exists(input_dir_rmats)) {
      warning("rMATS output directory not found: ", input_dir_rmats, ". Skipping rMATS filtering.")
    } else {
      event_types_rmats <- c("SE", "A3SS", "A5SS", "MXE", "RI")

      get_JCcounts <- function(x) {
        vals <- suppressWarnings(as.numeric(strsplit(as.character(x), ",")[[1]]))
        sum(vals[!is.na(vals)])
      }

      for (etype in event_types_rmats) {
        
        file_path <- file.path(input_dir_rmats, paste0(etype, ".MATS.JC.txt"))

        if (!file.exists(file_path)) {
          cat("  File not found:", basename(file_path), "\n")
          next
        }
        cat("  Filtering:", basename(file_path), "\n")

        tryCatch({
          df <- read_tsv(file_path, show_col_types = FALSE)

          
          required_cols_rmats <- c("GeneID", "geneSymbol", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "PValue", "FDR", "IncLevelDifference")
          if (!all(required_cols_rmats %in% colnames(df))) {
                 # Il tuo stop() ora cercherà solo le altre colonne
                 stop(paste0("One or more required rMATS columns missing. Found: ", paste(colnames(df), collapse=", ")))
          }
          
          df_filtered <- df %>%
           
            dplyr::select(
                GeneID, geneSymbol, 
                ID = 12, 
                IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2, PValue, FDR, IncLevelDifference
             ) %>%
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

          if(nrow(df_filtered) > 0){
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
} # End if(Rmats == TRUE)

# -----------------------------------------------------
# 2. SplAdder FILTERING
# -----------------------------------------------------
if(Spladder == TRUE) {
    cat("Processing SplAdder output...\n")

    input_dir_spladder <- file.path(base_dir, "spladder_out", "testing_treated_vs_Wildtype") 
    output_dir_spladder <- file.path(base_dir, "spladder_out", "Filtered_SplAdder")
    dir.create(output_dir_spladder, showWarnings = FALSE, recursive = TRUE)

    if (!dir.exists(input_dir_spladder)) {
      warning("SplAdder output directory not found: ", input_dir_spladder, ". Skipping SplAdder filtering.")
    } else {
      spladder_files <- list.files(input_dir_spladder, pattern = "^test_results_C3_.*\\.tsv$", full.names = TRUE)
      spladder_files <- spladder_files[!grepl("\\.gene_unique\\.", spladder_files)]

      if (length(spladder_files) == 0) {
          cat("  No SplAdder result files (*_C3_*.tsv) found in", input_dir_spladder, "\n")
      }

      for (file_path in spladder_files) {
        cat("  Filtering:", basename(file_path), "\n")

        event_type_spladder <- str_match(basename(file_path), "test_results_C3_(.+)\\.tsv")[,2]
        if (is.na(event_type_spladder)) event_type_spladder <- "unknown"

        tryCatch({
          df <- read_tsv(file_path, show_col_types = FALSE)

          
          required_cols <- c("p_val_adj", "dPSI", "mean_event_count_A", "mean_event_count_B", "p_val")
          if (!all(required_cols %in% colnames(df))) {
              warning("    Skipping file ", basename(file_path), ": missing one or more required columns (p_val_adj, dPSI, mean_event_count_A, mean_event_count_B, p_val). Available columns: ", paste(colnames(df), collapse=", "))
              next
          }

          df_filtered <- df %>%
            
            dplyr::mutate(
                TotalMeanCount = mean_event_count_A + mean_event_count_B
            ) %>%
            dplyr::filter(
              !is.na(p_val) & p_val < pval_cutoff
            ) %>%
            dplyr::select(any_of(c("gene_id", "gene_name", "event_id", "chrm", "exon_pos",
                            "p_val", "p_val_adj", "dPSI", "TotalMeanCount")))

          if(nrow(df_filtered) > 0) {
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
    cat("✔ SplAdder filtering complete.\n")
} else {
    cat("\nSpladder set to FALSE. Skipping SplAdder filtering.\n")
} 


# -----------------------------------------------------
# 3. RSEM FILTERING
# -----------------------------------------------------
if(RSEM == TRUE) {
    cat("Processing RSEM output...\n")

    input_dir_rsem_base <- file.path(base_dir, "Transcriptome", "counts")
    output_dir_rsem <- file.path(base_dir, "Transcriptome", "Filtered_RSEM")
    dir.create(output_dir_rsem, showWarnings = FALSE, recursive = TRUE)

    sample_dirs <- list.dirs(input_dir_rsem_base, full.names = TRUE, recursive = FALSE)

    if (length(sample_dirs) == 0) {
      warning("No sample directories found in: ", input_dir_rsem_base, ". Skipping RSEM filtering.")
    } else {

      for (sample_dir in sample_dirs) {
        sample_id <- basename(sample_dir)
        cat("  Filtering RSEM results for sample:", sample_id, "\n")

        sample_output_dir <- file.path(output_dir_rsem, sample_id)
        dir.create(sample_output_dir, showWarnings = FALSE, recursive = TRUE)

        # --- Filter GENES results ---
        genes_file <- file.path(sample_dir, paste0(sample_id, ".genes.results"))
        if (file.exists(genes_file)) {
          tryCatch({
              df_genes <- read_tsv(genes_file, show_col_types = FALSE)
              if(!all(c("expected_count", "TPM") %in% colnames(df_genes))) {
                  stop("Columns 'expected_count' or 'TPM' missing in genes results.")
              }

              df_genes_filtered <- df_genes %>%
               
                dplyr::filter(
                  !is.na(expected_count) & expected_count > count_cutoff,
                  !is.na(TPM) & TPM > tpm_cutoff
                )

              if(nrow(df_genes_filtered) > 0) {
                  out_name_genes <- paste0(sample_id, ".genes.results_filtered.tsv")
                  write_tsv(df_genes_filtered, file.path(sample_output_dir, out_name_genes))
                  cat("    Filtered genes results saved. Kept", nrow(df_genes_filtered), "genes.\n")
              } else {
                  cat("    No genes passed filtering.\n")
              }
          }, error = function(e) {
              cat("    Error processing file:", basename(genes_file), "-", e$message, "\n")
          })
        } else {
            cat("    File not found:", basename(genes_file), "\n")
        }

        # --- Filter ISOFORMS results ---
        isoforms_file <- file.path(sample_dir, paste0(sample_id, ".isoforms.results"))
        if (file.exists(isoforms_file)) {
          tryCatch({
              df_isoforms <- read_tsv(isoforms_file, show_col_types = FALSE)
              if(!all(c("expected_count", "TPM") %in% colnames(df_isoforms))) {
                  stop("Columns 'expected_count' or 'TPM' missing in isoforms results.")
              }

              df_isoforms_filtered <- df_isoforms %>%
                
                dplyr::filter(
                  !is.na(expected_count) & expected_count > count_cutoff,
                  !is.na(TPM) & TPM > tpm_cutoff
                )

              if(nrow(df_isoforms_filtered) > 0) {
                  out_name_isoforms <- paste0(sample_id, ".isoforms.results_filtered.tsv")
                  write_tsv(df_isoforms_filtered, file.path(sample_output_dir, out_name_isoforms))
                  cat("    Filtered isoforms results saved. Kept", nrow(df_isoforms_filtered), "isoforms.\n")
              } else {
                  cat("    No isoforms passed filtering.\n")
              }
          }, error = function(e) {
              cat("    Error processing file:", basename(isoforms_file), "-", e$message, "\n")
          })
        } else {
            cat("    File not found:", basename(isoforms_file), "\n")
        }
      } 
    } 
    cat("RSEM filtering complete.\n")
} else {
    cat("\nRSEM set to FALSE. Skipping RSEM filtering.\n")
}


# -----------------------------------------------------
# Final message
# -----------------------------------------------------
if(Rmats == TRUE || Spladder == TRUE || RSEM == TRUE) {
    cat("\n--- All selected filtering steps finished! ---\n")
}