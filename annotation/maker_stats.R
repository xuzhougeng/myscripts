#!/usr/bin/env Rscript

# ----------------------------------------------------------#
# This scripts will generate statistical plot based on      #
# AED and QI socre                                          #
# ----------------------------------------------------------#

# first running the following command to get the input data
# prefix=
# seqkit seq -n ${prefix}.all.maker.proteins.fasta | \
#  perl -n -e 'chomp;($AED,$QI)=(split/ /,$_)[2,4]; $AED=~s/AED://; $QI=~s/QI://; $QI=~s/\|/\t/g; ;print "$AED\t$QI\n";' > stat.txt

# column of stat.txt
#1. AED
#2. Length of the 5' UTR
#3. Fraction of splice sites confirmed by an EST alignment (-c)
#4. Fraction of exons that overlap an EST alignmetn(-e)
#5. Fraction of exons that overlap EST or Protein alignments(-o)
#6. Fraction of splice site confrimed by a ab-initio prediction(-a)
#7. Fraction of exons that overlap a ab-initio prediction(-t)
#8. Number of exons in the mRNA
#9. length of the 3' UTR
#10. Length of the protein sequence produced by the mRNA (-l)

suppressMessages(library(ggplot2))
suppressMessages(library(scales))

spec = matrix(c(
       'help', 'h', 0, "logical",
	   "file", 'f', 1, 'character'
	    ), byrow=TRUE, ncol=4)
opt = getopt::getopt(spec)

# if help was asked for print a friendly message 
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
	  cat(getopt::getopt(spec, usage=TRUE))
  q(status=1)
}


df <- read.table(opt$file, head=FALSE, sep="\t")

df <- within(df, {
  group <- NA
  group[V1 == 0] <- "0"
  group[V1 > 0 & V1 <=  0.25] <- "0~0.25"
  group[V1 > 0.25 & V1 <= 0.75] <- "0.25~0.75"
  group[V1 > 0.75 & V1 < 1] <- "0.75~1"
  group[V1 == 1] <- "1"
})

df$group <- factor(df$group, levels = c("0","0~0.25","0.25~0.75","0.75~1","1" ))

# AED count of each group
pdf("AED_distribution.pdf")

pie_df <- as.data.frame(table(df$group))
ggplot(pie_df, aes( x="", y=Freq, fill=Var1)) + 
  geom_bar(stat = "identity") +  
  geom_text(aes(y=Freq/length(Freq)+
                  c(0,cumsum(Freq)[-length(Freq)]),
                label=percent(Freq/sum(Freq))),size=5)+
  coord_polar("y",start=1) + 
  scale_fill_manual(values=c("#323695", "#8fc3dd","#fffdbf","#f88d51","#a51626")) +
  labs(x="",y="") + 
  theme_minimal() +
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
		panel.grid = element_blank()
  )
dev.off()

# exon number of each group
pdf("exon_distribution_of_each_group.pdf", width = 10, height=5)
exon_df <- as.data.frame(table(df$group,df$V8))

ggplot(exon_df, aes(x=Var2, y=Freq, fill=Var1)) + geom_bar(stat="identity") +
  scale_fill_manual(values=c("#323695", "#8fc3dd","#fffdbf","#f88d51","#a51626")) +
  theme_bw() +
  labs(x = "exon number", y= "count", fill="group")
dev.off()


# Fraction of EST overlap
pdf("est_overlap_fraction_of_each_group.pdf", width=5, height=10)
est_df <- as.data.frame(table(df$group,df$V3))
ggplot(est_df, aes(x=Var2, y=Freq, fill=Var1)) + geom_bar(stat="identity") +
  scale_fill_manual(values=c("#323695", "#8fc3dd","#fffdbf","#f88d51","#a51626")) +
  labs(x = "Fraction of splice sites confirmed by an EST alignment", y= "count", fill="group") + coord_flip() +
  scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) + theme_bw()

dev.off()
