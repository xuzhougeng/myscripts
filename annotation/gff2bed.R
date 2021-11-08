
gff <- rtracklayer::import.gff("GCF_000092025.1_ASM9202v1_genomic.gff")
gff <- gff[gff$type == "gene", ]
geneID <- substring(unlist(gff$Dbxref),8)
df <- data.frame(seqname = GenomeInfoDb::seqnames(gff), 
				 start = BiocGenerics::start(gff) - 1, 
				 end = BiocGenerics::end(gff), 
				 name = paste0(gff$Name,"_", geneID), 
				 score= ".", strand = BiocGenerics::strand(gff))
write.table(df, file="C58_gene.bed", sep="\t", quote=FALSE, col.names=FALSE, row.names = FALSE)

df$start <- df$start - 2000
df$start <- ifelse(df$start < 0, 0, df$start)
df$end <- df$end + 2000
write.table(df, file="C58_2k_gene.bed", sep="\t", quote=FALSE, col.names=FALSE, row.names = FALSE)

df$start <- df$start - 3000
df$start <- ifelse(df$start < 0, 0, df$start)
df$end <- df$end + 3000
write.table(df, file="C58_3k_gene.bed", sep="\t", quote=FALSE, col.names=FALSE, row.names = FALSE)


#gene <- mygene::getGenes(geneID)
#df2 <- data.frame(symbol = gene$symbol, descrption = gene$name, ncbi_id = gene$query)
#write.csv(df2, "C58_description.csv", row.names=F)
