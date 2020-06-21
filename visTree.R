#!/usr/bin/env Rscript

suppressMessages(library(ggtree))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

tree <- read.tree(args[1])

p <- ggplot(tree, aes(x, y)) + 
  geom_tree() + 
  theme_tree() +
  xlab("") + 
  ylab("") +
  geom_text(aes(label=label), hjust=1.0, vjust=1.5) +
  coord_cartesian(clip="off") +
  theme_tree2(plot.margin=margin(6,120,6,6))


fn <- gsub("\\..*","\\.png",basename(args[1]))
ggsave(fn,p)
