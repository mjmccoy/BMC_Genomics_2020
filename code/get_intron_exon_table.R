# Function get_intron-exon_table uses function get_attributes to pull multiple attributes for each species.
get_intron_exon_table <- function(species, clade, values = "protein_coding", filters = "biotype", output_dir = "./"){
  species <- gsub("\\_.*", "", species)
  species <- paste(
    species,
    ifelse(
      clade == "chordata",
      "_gene_ensembl",
      "_eg_gene"
    ),
    sep = ""
  )
  filename <- paste(
    output_dir,
    species,
    "_intron_exon",
    ".RData", 
    sep = ""
  )
  if(!file.exists(filename)){
    data <- get_attributes(
      species = species,
      clade = clade,
      filters = filters,
      values = values,
      attributes = c(
        "ensembl_gene_id",
        "ensembl_transcript_id",
        "chromosome_name",
        "start_position",
        "end_position",
        "exon_chrom_start",
        "exon_chrom_end",
        "rank",
        "strand"
      )
    )
    data$gene_length <- data$end_position - data$start_position
    data$exon_length <- data$exon_chrom_end - data$exon_chrom_start
    data$intron_length <- NA
    data <- data[order(data$ensembl_transcript_id, data$exon_chrom_start),]
    for(i in 1:nrow(data)){
      if(data$strand[i] > 0){
        if(i != nrow(data)){
          if(data$ensembl_transcript_id[i] == data$ensembl_transcript_id[i + 1]){
            data$intron_length[i] <- data$exon_chrom_start[i + 1] - data$exon_chrom_end[i]
          }
        }
      } else {
        if(i != 1){
          if(data$ensembl_transcript_id[i] == data$ensembl_transcript_id[i - 1]){
            data$intron_length[i] <- data$exon_chrom_start[i] - data$exon_chrom_end[i - 1]
          }
        }
      }
    }
    saveRDS(data, filename)
    return(data)
  } else {
    data <- readRDS(filename)
    return(data)
  }
}