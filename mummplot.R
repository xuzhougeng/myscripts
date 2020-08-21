
# modified from https://jmonlong.github.io/Hippocamplus/2017/09/19/mummerplots-with-ggplot2/
# use base R function to replace the dplyr/tidyr


# read delta file from MUMmer 
read.delta <- function(deltafile){
 
   # read all lines into memory
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  
  # remove first lines
  lines = lines[-1]
  
  # get length of each row
  line_list = strsplit(lines, ' ')
  lines_len = as.numeric(lapply(line_list, length))
  
  # group row by row length
  line_list = line_list[lines_len != 1]
  lines_len = lines_len[lines_len != 1]
  head_pos = which(lines_len == 4)
  
  head_id = rep(head_pos, c(head_pos[-1], length(line_list)+1)-head_pos)
  mat = matrix(as.numeric(unlist(line_list[lines_len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  
  res$qid = unlist(lapply(line_list[head_id[lines_len==7]], '[', 2))
  res$rid = gsub('^>', '', unlist(lapply(line_list[head_id[lines_len==7]], '[', 1)) )
  
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

# filter the MUMmer
# only contigs with at least one aligned segment larger than a minimum size
# Smaller alignment in these contigs are kept if in the same range as the large aligned segments
filter.mum <- function(df, min_len=1000, flanks=1e4){
  
  # filter by length
  tmp <- df[ abs(df$re - df$rs) > min_len, ]
  
  # get the min, max, median of each combination of qid and rid
  fac <- as.factor(paste(tmp$qid, tmp$rid,sep = "-"))
  min_max_list <- lapply(split(tmp, fac ), function(x) { 
    c(min(x$qs)-flanks, max(x$qe)+flanks, median(x$rs) ) 
    } )
  coord <- do.call(rbind, min_max_list)
  coord <- as.data.frame.matrix(coord, stringsAsFactors = FALSE)
  colnames(coord) <- c("qsL", "qeL", "rs")
  coord <- coord[ order(coord$rs, decreasing = TRUE), ]
  coord$qid <- gsub("-.*","", row.names(coord))
  coord$rid <- gsub(".*-","", row.names(coord))
  coord$qid <- factor(coord$qid, levels = unique(coord$qid))
  coord <- coord[, -which(names(coord) %in% c("rs"))]
  
  df <- merge(df, coord, by = c("qid", "rid"))
  
  df <- df[ (df$qs > df$qsL) & (df$qe < df$qeL),    ]
  df$qid <- factor(df$qid, levels = unique(df$qid))
  df <- df[ , -which(names(df) %in% c("qsL", "qeL")) ]
  df
}

#
diag.mum <- function(df) {
  # find best qid order
  fac <- as.factor(paste(df$qid, df$rid,sep = "-"))
  rid <- lapply(split(df, fac ), function(x) { 
    c( sum(abs(x$qe - x$qs) ), weighted.mean(x$rs, abs(x$qe- x$qs)) ) 
  } )
  rid <- do.call(rbind, rid)
  rid <- as.data.frame.matrix(rid, stringsAsFactors = FALSE)
  colnames(rid) <- c("base", "rs")
  rid <- rid[ order(rid$base),  ]

  # find best qid  
  
}



p <- ggplot(mumgp.filt, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
  geom_segment() +
  geom_point(alpha=1, size=.5) + 
  facet_grid(qid~rid, scales='free', space='free', switch='y') +
  theme_bw() 
  
p

p <- p + theme( panel.spacing = unit(0, "lines") )
p <- p + theme( strip.background=element_blank()) 

p <- p +   theme( axis.text.y=element_blank(), axis.ticks.y=element_blank() ) +
  theme( axis.text.x=element_blank(), axis.ticks.x=element_blank() ) 
  

p + theme( panel.grid.minor  = element_blank()) 

p <- p + theme( panel.border = element_blank())   

p + theme(axis.line.y.right = element_line(colour="black"))


  
p +   theme( axis.line.y = element_line())
  
xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1')


ggsave("output.pdf", p)



element_grob.element_custom <- function(element, ...)  {
  # line (1,0)->(1,1)
  grid::segmentsGrob(x0=1, y0=0, x1=1,y1=1, gp=gpar(lwd=2))
}

border_custom <- function(...){
  structure(
    list(...), # this ... information is not used, btw
    class = c("element_custom","element_blank", "element") # inheritance test workaround
  ) 
  
}

p + theme(panel.border=border_custom())

