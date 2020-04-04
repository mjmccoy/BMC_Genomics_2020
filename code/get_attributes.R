# Function get_attributes pulls specified attributes (i.e. gene and exon coordinates)
# for specified species.
get_attributes <- function(species, clade, filters, values, attributes){
  getBM(
    attributes = attributes,
    filters = filters,
    values = values,
    mart = useMart(
      biomart = ifelse(
        clade == "chordata",
        "ensembl",
        paste(clade, "_mart", sep = "")
      ),
      host = ifelse(
        clade == "chordata",
        "ensembl.org",
        paste(clade, ".ensembl.org", sep = "")
      ),
      dataset = species)
    )
}
