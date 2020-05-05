# the function get_neighbor_distances retrieves the average distance
# between a gene and its nearest neighbor genes on the same chromosome.

get_neighbor_distances <- function(gene_table){
  gene_dist_table <- gene_table[,1:4]
  gene_dist_table <- unique(gene_dist_table)
  gene_dist_table <- gene_dist_table[order(
    gene_dist_table$chromosome_name,
    gene_dist_table$start_position),]
  gene_dist_table$length <- gene_dist_table$end_position - gene_dist_table$start_position
  gene_dist_table$space <- NA
  for (chromosome in unique(gene_dist_table$chromosome_name)){
    chrom <- subset(gene_dist_table, chromosome_name == chromosome)
    for (i in 1:nrow(chrom)){
      if (chrom[i,"start_position"] == min(chrom$start_position)){
        space <- chrom[i + 1, "start_position"] - chrom[i, "end_position"]
        gene_dist_table[gene_dist_table$ensembl_gene_id == chrom$ensembl_gene_id[i],"space"] <- space
      } else if (chrom[i, "end_position"] == max(chrom$end_position)){
        space <- chrom[i, "start_position"] - chrom[i - 1, "end_position"]
        gene_dist_table[gene_dist_table$ensembl_gene_id == chrom$ensembl_gene_id[i],"space"] <- space
      } else {
        lower <- chrom[i, "start_position"] - chrom[i - 1, "end_position"]
        upper <- chrom[i + 1, "start_position"] - chrom[i, "end_position"]
        space <- (lower + upper)/2
        gene_dist_table[gene_dist_table$ensembl_gene_id == chrom$ensembl_gene_id[i],"space"] <- space
      }
    }
  }
  gene_dist_table$space[gene_dist_table$space <= 0] <- 0
  return(gene_dist_table)
}