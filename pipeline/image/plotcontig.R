#!/usr/bin/R

library(ggplot2)

# This function plots the contig alone
plotcontig <- function(CONTIGNAME, GENETABLE){
  # Get the contig we care about (or want to look at)
  contig <- subset(GENETABLE, GENETABLE$contig == CONTIGNAME)
  
  # Sort taxa by scores, so legend only outputs top 10
  contig_sorted <- contig[order(-contig$score), ]
  
  # Get the length of the contig
  contiglen = contig$contiglen[1] + 1
  
  # Plot the groups
  plot1 <- ggplot(contig, aes(x=groupstart, y=0)) + geom_point() + geom_segment(aes(x=groupstart, y=0, xend=groupend, yend=0)) + xlim(0, contiglen)
  
  # Plot the BLAST results
  plot2 <- plot1 + geom_point(data=contig, aes(x=genestart, y=score, colour=taxa)) + geom_point() + geom_segment(aes(x=genestart, y=score, xend=geneend, yend=score, colour=taxa), size=1.2) + geom_line(data=contig, aes(x=genestart, y=score, colour=taxa))
  
  # Add labels
  plot2 + labs(x='Contig Coordinates', y='Score', title=CONTIGNAME, colour='Taxa') + scale_colour_discrete(breaks=contig_sorted$taxa[1:10])
}

#This function plots the synthetic contig answers
plotanswercontig <- function(CONTIGNAME, ANSWERKEY){
  # Get the contig we care about (or want to look at)
  answer <- subset(ANSWERKEY, ANSWERKEY$contig == CONTIGNAME)
  
  # Get the length of the contig
  contiglen = answer$contiglen[1] + 1
  
  # Plot coordinates
  plot1 <- ggplot(answer, aes(x=genestart, y=score, colour=taxa)) + geom_point() + geom_segment(aes(x=genestart, y=score, xend=geneend, yend=score, colour=taxa), size=1.2)
  
  # Label axes
  plot1 + labs(x='Contig Coordinates', y='Score', title=CONTIGNAME, colour='Taxa')
}


#This function plots the contig along with the answer
plotall <- function(CONTIGNAME, GENETABLE, ANSWERKEY){
  # Get the contig we care about (or want to look at)
  answer <- subset(ANSWERKEY, ANSWERKEY$contig == CONTIGNAME)
  contig <- subset(GENETABLE, GENETABLE$contig == CONTIGNAME)
  total <- rbind(answer, contig)
  
  # Sort taxa by scores, so legend only outputs top 10
  contig_sorted <- contig[order(-contig$score), ]
  
  # Get the length of the contig
  answerlen <- total$contiglen[1] + 1
  
  # Plot the groups
  plot1 <- ggplot() + geom_point(data=contig, aes(x=groupstart, y=0)) + geom_segment(data=contig, aes(x=groupstart, y=0, xend=groupend, yend=0)) + xlim(0, answerlen)
  
  # Plot the BLAST results
  plot2 <- plot1 + geom_point(data=total, aes(x=genestart, y=score, colour=taxa)) + geom_segment(data=total, aes(x=genestart, y=score, xend=geneend, yend=score, colour=taxa), size=1.2) + geom_line(data=total, aes(x=genestart, y=score, colour=taxa))
  
  # Plot the answers 
  #plot3 <- plot2 + geom_point(data=answer, aes(x=genestart, y=score, colour=taxa)) + geom_segment(data=answer, aes(x=genestart, y=score, xend=geneend, yend=score, colour=taxa), size=1.2)
  
  # Add labels
  plot4 <- plot2 + labs(x='Contig Coordinates', y='Score', title=CONTIGNAME, colour='Taxa') + scale_colour_discrete(breaks=contig_sorted$taxa[1:15])
  
  # Facet grid
  plot4 + facet_grid(status ~ .)
}

#This function prints the contig based on what you need printed
saveplot <- function(CONTIG, GENETABLE, ANSWERKEY, NAME){
  myplot <- plotall(CONTIG, GENETABLE, ANSWERKEY)
  filename <- paste(NAME, '.png', sep='')
  #pdf(filename, width=11, height=8.5, useDingbats = FALSE) 
  png(filename, width=1100, height=850, res=120)
  print(myplot)
  dev.off()
}
