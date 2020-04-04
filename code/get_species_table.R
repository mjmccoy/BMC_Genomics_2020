get_species_table <- function(clade.datasets){
  clade <- gsub("\\..*", "",  deparse(substitute(clade.datasets)))
  # Make empty dataframe
  clade.species.lengths <- data.frame(
    species = sapply(1:nrow(clade.datasets), FUN = function(x){
      ifelse(clade == "chordata", clade.datasets$species[x], clade.datasets$description[x])
    }),
    median_gene_length = NA,
    median_long_gene_length = NA,
    median_exonic_length = NA,
    median_long_gene_exonic_length = NA,
    median_intronic_length = NA,
    median_long_gene_intronic_length = NA,
    number_of_exons = NA,
    number_of_long_gene_exons = NA,
    clade = rep(
      if(clade == "plants"){
        "Plantae"
      } else if(clade == "protists"){
        "Protista"
      } else if(clade == "metazoa"){
        "Metazoa"
      } else if(clade == "fungi"){
        "Fungi"
      } else{
        "Chordata"
      }, nrow(clade.datasets)),
    Taxon.ID = clade.datasets$Taxon.ID)
  # Loop through all species in clade to obtain feature estimates
  for (i in 1:nrow(clade.datasets)) {
    species <- clade.datasets[i,1]
    gene_table <- get_gene_table(species, clade, output_dir = "../data/gene_tables/")
    # Add individual species estimates to clade dataframe
    clade.species.lengths[i,2] <- median(gene_table$gene_length)
    clade.species.lengths[i,3] <- median(gene_table$gene_length[gene_table$gene_length > quantile(gene_table$gene_length, probs = 0.90)])
    clade.species.lengths[i,4] <- median(gene_table$exonic_length)
    clade.species.lengths[i,5] <- median(gene_table$exonic_length[gene_table$gene_length > quantile(gene_table$gene_length, probs = 0.90)])
    clade.species.lengths[i,6] <- median(gene_table$intronic_length)
    clade.species.lengths[i,7] <- median(gene_table$intronic_length[gene_table$gene_length > quantile(gene_table$gene_length, probs = 0.90)])
    clade.species.lengths[i,8] <- median(gene_table$number_of_exons)
    clade.species.lengths[i,9] <- median(gene_table$number_of_exons[gene_table$gene_length > quantile(gene_table$gene_length, probs = 0.90)])
    print(paste(i,"/", nrow(clade.datasets), sep = " "))
  }
  saveRDS(clade.species.lengths, file = paste("../data/clade_tables/", clade, "table.RDS"))
  return(clade.species.lengths)
}