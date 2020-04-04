# function LGE_data takes expression data and associates it with gene length data.
get_LGE <- function(expression_file, species, clade, output_dir = "./", quantile_normalized = F){
  data.df <- read.table(file = expression_file, sep = "\t", header = TRUE)
  data.df[is.na(data.df)] <- 0
  if (quantile_normalized == T){
    data.df[,-c(1:2)] <- preprocessCore::normalize.quantiles(as.matrix(data.df[,-c(1:2)]))
  }
  filters <- c("biotype", "ensembl_gene_id")
  values <- list("protein_coding", data.df$Gene.ID)
  gene_table <- get_gene_table(
    species = species,
    clade = clade,
    filters = filters,
    values = values,
    output_dir = output_dir)
  gene_table$gene_length <- gene_table$end_position - gene_table$start_position
  data.df$gene_length <- gene_table$gene_length[match(data.df$Gene.ID, gene_table$ensembl_gene_id)]
  data.df <- na.omit(data.df)
  # melt data for plotting
  melt.data <- melt(data.df, id.vars = c("Gene.ID", "Gene.Name", "gene_length"))
  melt.data$species <- species
  return(melt.data)
}