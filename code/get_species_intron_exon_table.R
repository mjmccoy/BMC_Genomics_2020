get_species_intron_exon_table <- function(clade.datasets){
  clade <- gsub("\\..*", "",  deparse(substitute(clade.datasets)))
  # Make empty dataframe
  clade.species.lengths <- data.frame(
    species = sapply(1:nrow(clade.datasets), FUN = function(x){
      ifelse(clade == "chordata", clade.datasets$species[x], clade.datasets$description[x])
    }),
    median_gene_length = NA,
    median_exon_length = NA,
    median_intron_length = NA,
    median_number_of_exons = NA,
    number_of_genes = NA,
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
    species <- clade.datasets$dataset[i]
    tryCatch({
      intron_exon_table <- get_intron_exon_table(species, clade, output_dir = "./intron_exon_tables/")
      # Add individual species estimates to clade dataframe
      clade.species.lengths$median_gene_length[i] <- median(intron_exon_table$gene_length, na.rm = T)
      clade.species.lengths$median_exon_length[i] <- median(intron_exon_table$exon_length, na.rm = T)
      clade.species.lengths$median_intron_length[i] <- median(intron_exon_table$intron_length, na.rm = T)
      clade.species.lengths$median_number_of_exons[i] <- intron_exon_table %>%
        group_by(ensembl_gene_id) %>% 
        summarise(number_of_exons = max(rank, na.rm = T)) %>% 
        summarise(median_number_of_exons = median(number_of_exons, na.rm = T)) %>% 
        as.numeric
      clade.species.lengths$number_of_genes[i] <- length(unique(intron_exon_table$ensembl_gene_id))
      print(paste(i,"/", nrow(clade.datasets), sep = " "))
    }, error = function(e){cat(paste("ERROR :", i, "/", nrow(clade.datasets), species, sep = " "), "\n")})
  }
  saveRDS(clade.species.lengths, file = paste("./intron_exon_tables/clade_tables/", clade, "_intron_exon_table.RDS", sep = ""))
  return(clade.species.lengths)
}