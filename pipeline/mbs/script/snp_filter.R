# R script for filter and plot

df <- read.csv(snakemake@input[[1]], header = TRUE,check.names = TRUE )
wildtype <- gsub("-",".",snakemake@params[['wildtype']], fixed=TRUE)
mutation <- gsub("-",".",snakemake@params[['mutation']], fixed=TRUE)

# Calculate the allele frequency of alternative base
# AD: reference
# AD.1: alternative
df$wildtype_AF <- df[paste0(wildtype, '.AD.1')] / (df[paste0(wildtype, '.AD')] + df[paste0(wildtype, '.AD.1')])
df$wildtype_AF <- df$wildtype_AF[[1]]
df$mutation_AF <- df[paste0(mutation, '.AD.1')] / (df[paste0(mutation, '.AD')] + df[paste0(mutation, '.AD.1')])
df$mutation_AF <- df$mutation_AF[[1]]

# Calculate the delta snp index
df$Index <- df$mutation_AF - df$wildtype_AF

# filter the snp below the threshold
threshold <- snakemake@params[['threshold']]
filtered <- dplyr::filter(df, abs(Index)> threshold  & ! stringr::str_detect(EFF, "MODIFIER"))

# annotation the remain gene with gene description downloaded from TAIR10
anno_df <- read.delim("script/gene_description_20131231.tsv", sep="\t", header=FALSE, quote="")
filtered$gene_name <- stringr::str_extract(filtered$EFF,"AT\\dG\\d{5}\\.\\d")
filtered <- dplyr::left_join(filtered, anno_df, by=c("gene_name"="V1"))

# write out the candidate snp site
write.csv(filtered, snakemake@output[['csv']])

# plot
library(ggplot2)
library(magrittr)
trend_plot <- function(raw_df, flt_df, chrom){
    p1 <- raw_df %>% dplyr::filter(CHROM == chrom) %>%
        ggplot(aes(x=POS,y=Index)) +
        geom_point(size=2) +
        scale_x_continuous(breaks = seq(0,3e7,1000000)) +
        scale_y_continuous(breaks = seq(-1,1,0.1), limits = c(-1,1))
    d2 <- flt_df %>% dplyr::filter(CHROM == chrom)
    p2 <- p1 + geom_point(data=d2, aes(x=POS, y =Index),color='red',size=2) +
        geom_text(data=d2, aes(x=POS, y=Index,label=gene_name),angle=40,hjust=-0.1) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
    return(p2)

}

pdf(snakemake@output[["pdf"]],width=20,height=10)
# snp index plot separately 
for (chrom in unique(df$CHROM)){
    p <- trend_plot(df, filtered, chrom)
    print(p)
}

# snp index plot together
chr_len <- read.delim("script/Athaliana_chr_len.tsv", sep="\t", header = T)
snp_flt <- dplyr::filter(df, CHROM != "ChrCh", CHROM != "ChrM")
plot_pos <- apply(snp_flt, 1, function(df1,df2){
  as.numeric(df1[2]) + df2$start[which(df2$chrom==df1[1])]-1
}, df2=chr_len)

snp_flt$plot_pos <- unlist(plot_pos)

ggplot(snp_flt, aes(x=plot_pos,y=Index)) +
    geom_line(size=0.75, colour='black') +   
    theme_minimal() +
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
    geom_segment(data=chr_len,aes(x=start,xend=end,y=0,yend=0,colour=chrom),size=2) +
    ylim(-1,1)

dev.off()
