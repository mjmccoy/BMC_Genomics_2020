# Function get_gene_table uses function get_attributes to pull multiple attributes for each species.
# If simplified is TRUE, adds exonic/intronic length and exon number to the
# specified species gene table.
get_gene_table <- function(species, clade, values = "protein_coding", filters = "biotype", simplified = T, output_dir = "./"){
  species <- gsub("\\_.*", "", species)
  species <- paste(
    species,
    ifelse(
      clade == "chordata",
      "_gene_ensembl",
      "_eg_gene"),
    sep = "")
  filename <- paste(output_dir, species, ".RData", sep = "")
  if(!file.exists(filename)){
    genes <- get_attributes(
      species = species,
      clade = clade,
      filters = filters,
      values = values,
      attributes = c(
        "ensembl_gene_id",
        "chromosome_name",
        "start_position",
        "end_position")
      )
    exons <- get_attributes(
      species = species,
      clade = clade,
      filters = filters,
      values = values,
      attributes = c(
        "ensembl_gene_id",
        "rank",
        "exon_chrom_start",
        "exon_chrom_end")
      )
    if(simplified == F){
      genes <- merge(
        genes,
        exons,
        by = "ensembl_gene_id")
      saveRDS(genes, file = filename)
      return(genes)
    } else{
      # Calculate exon numbers and exonic length
      exons$ensembl_gene_id <- gsub(".", "_",exons$ensembl_gene_id, fixed = TRUE)
      exonic_lengths <- split(exons[,3:4],f=exons$ensembl_gene_id) # split up exon positions by gene
      exonic_lengths <- lapply(exonic_lengths, FUN = function(x) {interval_union(Intervals(x))}) # Find union of exons
      exonic_lengths <- lapply(exonic_lengths, FUN = as.data.frame)
      exonic_lengths <- do.call(rbind,exonic_lengths)
      exonic_lengths$ensembl_gene_id <- sub('\\..*', '', rownames(exonic_lengths))
      exon_numbers <- as.data.frame(table(exonic_lengths$ensembl_gene_id)) # calculate exon numbers
      exonic_lengths$exonic_length <- exonic_lengths$V2 - exonic_lengths$V1 # calculate exon lengths
      exonic_lengths <- aggregate(exonic_lengths$exonic_length, by = list(exonic_lengths$ensembl_gene_id), FUN = sum) # calculate total exonic length for each gene
      colnames(exonic_lengths) <- c("ensembl_gene_id", "exonic_length")
      # Calculate gene lengths
      agg.genes <- aggregate(genes[,3:4], by = list(genes$ensembl_gene_id), FUN = function(x){mean(x, na.rm = TRUE)}) # if multiple gene start and stop positions, find the average for each gene
      colnames(agg.genes)[1] <- "ensembl_gene_id"
      agg.genes$gene_length <- agg.genes$end_position - agg.genes$start_position
      # Add exonic and intornic length to main dataframe
      agg.genes$ensembl_gene_id <- gsub(".", "_",agg.genes$ensembl_gene_id, fixed = TRUE)
      agg.genes$exonic_length <- exonic_lengths$exonic_length[match(agg.genes$ensembl_gene_id, exonic_lengths$ensembl_gene_id)]
      agg.genes$intronic_length <- agg.genes$gene_length - agg.genes$exonic_length
      agg.genes$number_of_exons <- exon_numbers$Freq[match(agg.genes$ensembl_gene_id, exon_numbers$Var1)]
      saveRDS(agg.genes, filename)
      return(agg.genes)
    }
  } else{
    genes <- readRDS(filename)
    return(genes)
  }
}