#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {
  stop("At least three arguments must be supplied", call.=FALSE)
} else if (length(args) >= 4){
  #input_annotation = args[1]
  input_annotation = "test_data/reference/mgi_lncRNA_annotation.MF.tsv"
  #output_dir = args[2]
  output_dir = "output/semantic"
  #reference = args[3]
  reference = "test_data/annotation/mgi_genes_annotation.MF.tsv"
  #ontology_type = args[4]
  ontology_type = "molecular_function"
  ontology_type_short = "MF"
  ontology_name = "org.Mm.eg.db"
  if(length(args) >= 5){
    ontology_name = args[5]
  }
}
output_dir
reference
ontology_type
ontology_name
library("GOSemSim")
#########
#########

#########
#########
#hsGO <- godata('org.Mm.eg.db', ont="MF")
hsGO <- godata(ontology_name, ont=ontology_type_short)
id2go_path = reference
id2go = read.table(id2go_path, sep='\t', head=FALSE)
ontology_type
#head(id2go[(id2go[ ,3] == ontology_type), ])
#length(id2go[ ,3])
"Listing tests"
#tests <- list.files(predictions_dir, full.names = TRUE)
#tests

prediction_similarity <- function(gene_df,m,c,onto,ref_dfs){
  gene_name <- as.vector(gene_df$V1[[1]])
  #reference_df <- id2go[(id2go[ ,1] == gene_name), ]
  
  #reference_gos <-  as.vector(reference_df$go_id)
  prediction_gos <-  as.vector(gene_df$V2)
  #mgoSim(prediction_gos, reference_gos, semData=hsGO, measure=m, combine=c)
  custom_mgoSim <- function(custom_ref){
    ref_name <- as.vector(custom_ref$V1[1])
    reference_gos <- as.vector(custom_ref$V2)
    c(gene_name,"bab", mgoSim(prediction_gos, custom_ref, semData=hsGO, measure=m, combine=c))
  }
  sapply(ref_dfs, custom_mgoSim)
}

prediction_file = input_annotation
m = "Wang"
c = "BMA"
prediction_df <- read.table(prediction_file, header=FALSE,sep='\t')
subset_df <-function(id){
  prediction_df[(prediction_df[ ,1] == id), ]
}

unique_names <- as.vector(unique(prediction_df$V1))
dfs <- lapply(unique_names,subset_df)

reference_df <- id2go
onto <- ontology_type
reference_df <- reference_df[(reference_df[ ,3] == onto), ]
subset_df_ref <-function(id){
  reference_df[(reference_df[ ,1] == id), ]
}
unique_names_ref <- as.vector(unique(reference_df$V1))
ref_dfs <- lapply(unique_names,subset_df_ref)

custom_prediction_similarity <- function(custom_df){
  prediction_similarity(custom_df,m,c,onto,ref_dfs)
}

similarity <- sapply(dfs, custom_prediction_similarity)

sims <- as.vector(similarity)
out_path <- paste(output_dir,
                  basename(input_annotation),
                  sep="/")
fileConn <-file(out_path)
writeLines(paste(sims[!is.na(sims)],sep='\n'), fileConn)
close(fileConn)