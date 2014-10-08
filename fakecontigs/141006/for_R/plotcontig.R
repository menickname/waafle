plotcontig <- function(contigNAME, combtable){
  info <- subset(combtable, combtable$V1 == contigNAME)
  infogenus <- subset(info, info$V4 == 'genus')
  infogenus_sorted <- infogenus[order(-infogenus$V5),]
  if (nrow(infogenus_sorted) > 10){
    cut_sorted <- infogenus_sorted[1:10,]
  } else{
    cut_sorted <- infogenus_sorted
  }
  basicgrid <- ggplot(info) + geom_point(aes(x=V6,y=V5), size=3) + geom_segment(aes(x=V6,y=V5,xend=V7,yend=V5), size=1.2) + facet_grid(V4 ~ .)
  coloredgrid <- basicgrid + geom_point(data=cut_sorted, aes(x=V6,y=V5,colour=V3), size=3) + geom_segment(data=cut_sorted, aes(x=V6,y=V5,xend=V7,yend=V5,colour=V3), size=1.2) + geom_line(data=cut_sorted, aes(x=V6,y=V5, colour=V3))
  coloredgrid + labs(x="Contig Coordinates (bp)", y="Scores", title=contigNAME) + scale_colour_discrete(name='genus', breaks=cut_sorted$V3)
}
