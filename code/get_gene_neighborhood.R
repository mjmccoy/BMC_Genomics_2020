# Function get_gene_neighborhood extracts sequences flanking gene of interest.
get_gene_neighborhood <- function(gene_table, ensembl_gene_id, neighborhood_length){
  gene <- gene_table[gene_table$ensembl_gene_id == ensembl_gene_id,]
  gene_length <- mean(gene$end_position) - mean(gene$start_position)
  diff <- neighborhood_length - gene_length
  neighborhood_start <- mean(gene$start_position) - diff/2
  neighborhood_end <- mean(gene$end_position) + diff/2
  neighborhood_table <- subset(
    gene_table,
    chromosome_name == unique(gene$chromosome_name) &
    exon_chrom_start > neighborhood_start &
    exon_chrom_end < neighborhood_end)
  # calculate relative genomic coordinates
  neighborhood_table$relative_start_position <- neighborhood_table$start_position - min(neighborhood_table$start_position)
  neighborhood_table$relative_end_position <- neighborhood_table$end_position - min(neighborhood_table$start_position)
  neighborhood_table$relative_exon_chrom_start <- neighborhood_table$exon_chrom_start - min(neighborhood_table$start_position)
  neighborhood_table$relative_exon_chrom_end <- neighborhood_table$exon_chrom_end - min(neighborhood_table$start_position)
  neighborhood_table$species <- deparse(substitute(gene_table)) # use name of gene_table argument as species name
  return(neighborhood_table)
}