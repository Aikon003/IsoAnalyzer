library(IsoformSwitchAnalyzeR)

library(BSgenome.Hsapiens.UCSC.hg38) 

library(ggplot2) 
# ---------------------------



args <- commandArgs(trailingOnly = TRUE)


if (length(args) == 0) {
  stop("Error: No path provided.\nUsage: Rscript ./isoanalyzer.R <output_folder_path>", call. = FALSE)
} else if (length(args) > 1) {
  stop("Error: Too many arguments provided. Please provide only the folder path.", call. = FALSE)
}


base_dir <- args[1]


if (!dir.exists(base_dir)) {
  stop(paste("Error: The specified directory does not exist:", base_dir), call. = FALSE)
}


cat("Base directory set to:", base_dir, "\n")


output_rds_file <- file.path(base_dir,"IsoformSwitchAnalyzer_Output", "switchList_part1.rds")
switchList <- readRDS(output_rds_file) 
analysis_dir <- file.path(base_dir, "Analysis")
protein_input_dir <- file.path(analysis_dir, "Protein")
plot_output_dir <- file.path(analysis_dir, "Plots_Finali")
dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)


cat("Filtering significant isoform switches...\n")

asl_analyzed_selected=extractTopSwitches(
     switchList,
     filterForConsequences=FALSE,
     extractGenes=FALSE,
     alpha=0.05,
     dIFcutoff = 0.1,
     n=Inf,
     inEachComparison=FALSE,
     sortByQvals=TRUE
 )
 isoform_txt <- file.path(base_dir, "Analysis", "isoforms_DTU.txt")
 write.table(asl_analyzed_selected, 'isoforms_DTU.txt', quote=F,sep='\t',row.names=F,col.names=T)

cat("Importing protein analysis...\n")

# Add CPAT analysis
asl_data = analyzeCPAT(
switchAnalyzeRlist   = switchList,
pathToCPATresultFile = file.path(protein_input_dir, "cpat_output.txt"),
codingCutoff         = 0.725, 
removeNoncodinORFs   = TRUE   
)

# Add PFAM analysis
asl_data = analyzePFAM(
switchAnalyzeRlist   = asl_data,
pathToPFAMresultFile = file.path(protein_input_dir, "pfam_output.txt"),
showProgress=FALSE
)

# Add SignalP analysis
asl_data = analyzeSignalP(
switchAnalyzeRlist       = asl_data,
pathToSignalPresultFile  = file.path(protein_input_dir, "signalp_output_summary.signalp5")
)

# Add IUPred2A analysis
asl_data = analyzeIUPred2A(
switchAnalyzeRlist        = asl_data,
pathToIUPred2AresultFile = file.path(protein_input_dir, "iupred2a_output.txt"),
showProgress = FALSE
)

cat("Predict Alternative Splicing...\n")
asl_data_altspl = analyzeAlternativeSplicing(switchAnalyzeRlist = asl_data, quiet=TRUE)

cat("Analyze consequences...\n")
asl_data_analyzed = analyzeSwitchConsequences( asl_data_altspl )

cat("Plotting...\n")

output_rds_for_plotting <- file.path(base_dir, "IsoformSwitchAnalyzer_Output", "switchList_for_plotting.rds")
cat("Saving switchList_for_plotting object to:",output_rds_for_plotting, "\n")
saveRDS(asl_data_analyzed, file = output_rds_for_plotting)



# 2. Extract top 10 significant genes
top_switches <- extractTopSwitches(asl_data_analyzed, n = 10, filterForConsequences = TRUE)


if(nrow(top_switches) > 0) {
  
  cat(paste("  Plotting the", nrow(top_switches), "most significant genes in separate PDF files...\n"))


  for (i in 1:nrow(top_switches)) {
    
    
    top_gene_name <- top_switches$gene_name[i]
    top_gene_id <- top_switches$gene_id[i]
    cond1 <- top_switches$condition_1[i]
    cond2 <- top_switches$condition_2[i]

    
    if (!is.na(top_gene_name) && top_gene_name != "") {
        gene_da_plottare <- top_gene_name
        cat("using gene name for plotting:", gene_da_plottare, "\n")
    } else {
        gene_da_plottare <- top_gene_id
        cat("using gene ID for plotting:", gene_da_plottare, "\n")
        
    }   


    file_name <- paste0('switchplot_', i, '_', top_gene_id, '.pdf')
    file_path <- file.path(plot_output_dir, file_name)

    cat(paste("    (", i, "/", nrow(top_switches), ") Plotting:", gene_da_plottare, "->", file_name, "\n"))

    # 4. Create plot for THAT gene
    pdf(file=file_path, width=20, height=15)
    switchPlot(
        asl_data_analyzed,
        gene = gene_da_plottare, 
        condition1 = cond1,
        condition2 = cond2,
        plotTopology = FALSE, 
        localTheme = theme_bw(base_size = 13)
    )
    dev.off()
    
  } 
  
} else {
  cat("  No significant switches with consequences found, skipping single gene plots.\n")
}

# Overview plot of q.value versus dIF
cat("  Plotting overview 1 (dIF vs q-value)...\n")
pdf(file=file.path(plot_output_dir, 'overview1_dIF_vs_qvalue.pdf'), width=6, height=6)
print( 
  ggplot(data=asl_data_analyzed$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
      aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
      size=1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_1) +
  scale_color_manual('Significant\nIsoform Switch', values = c('black','red')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()
)
dev.off()

# Overview plot of dIF versus gene fc
cat("  Plotting overview 2 (dIF vs Gene FC)...\n")
pdf(file=file.path(plot_output_dir, 'overview2_dIF_vs_GeneFC.pdf'), width=6, height=6)
print(
  ggplot(data=asl_data_analyzed$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(
      aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
      size=1
  ) + 
  facet_wrap(~ condition_1) +
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual('Significant\nIsoform Switch', values = c('black','red')) +
  labs(x='Gene log2 fold change', y='dIF') +
  theme_bw()
)
dev.off()




# Consequence summary
cat("  Plotting consequence summary...\n")
pdf(file=file.path(plot_output_dir, 'consequence_summary.pdf'), width=12, height=5)
extractConsequenceSummary(
  asl_data_analyzed,
  consequencesToAnalyze='all',
  plotGenes = FALSE,
  asFractionTotal = FALSE,
  returnResult=TRUE
)
dev.off()

# Consequence enrichment
cat("  Plotting consequence enrichment...\n")
pdf(file=file.path(plot_output_dir, 'consequence_enrichment.pdf'), width=12, height=6)
extractConsequenceEnrichment(
  asl_data_analyzed,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE,
  returnResult = TRUE
)
dev.off()




# Splicing enrichment
cat("  Plotting splicing enrichment...\n")
pdf(file=file.path(plot_output_dir, 'splicing_enrichment.pdf'), width=12, height=6)
extractSplicingEnrichment(
  asl_data_analyzed,
  returnResult = FALSE
)
dev.off()

# Splicing enrichment comparison
cat("  Plotting splicing enrichment comparison...\n")
pdf(file=file.path(plot_output_dir, 'splicing_enrichment_comparison.pdf'), width=12, height=6)
extractSplicingEnrichmentComparison(
  asl_data_analyzed,
  splicingToAnalyze = c('A3','MES','ATSS','ATTS'),
  returnResult = TRUE
)


# Consequence enrichment comparison
cat("  Plotting consequence enrichment comparison...\n")
pdf(file=file.path(plot_output_dir, 'consequence_enrichment_comparison.pdf'), width=12, height=6)
extractConsequenceEnrichmentComparison(
  asl_data_analyzed,
  consequencesToAnalyze=c('domains_identified','intron_retention','coding_potential'),
  analysisOppositeConsequence = TRUE,
  returnResult = TRUE
)

# Compare DTUs between comparisons (Venn Diagram)
cat("  Plotting overlap (Venn)...\n")

pdf(file=file.path(plot_output_dir, 'overlap_venn_diagram.pdf'))
extractSwitchOverlap(
  asl_data_analyzed,
  filterForConsequences=TRUE
)
dev.off()


cat("--- Plotting completed successfully! ---\n")
