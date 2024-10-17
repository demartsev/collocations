setwd("C:/Users/vdemartsev/My Cloud/HyraxPostDoc/Call_sequences")
library(igraph)


txt_files <- list.files(paste(getwd(), "/sep_females/", sep = ""), full.names = T)
pdf(file= "sample.pdf" , onefile = TRUE)

for (female in txt_files) {

input.matrix <- read.table(female, header=TRUE, sep="\t", quote="", comment.char="")

#input.matrix <- read.table(file.choose(), header=TRUE, sep="\t", quote="", comment.char="")
names(input.matrix) <- c("W_C", "Coll_Word")
input.matrix <- table(input.matrix$Coll_Word, input.matrix$W_C)

# computation
pearson.residuals <- as.data.frame.matrix(chisq.test(input.matrix, correct=FALSE)$residuals)
output.table <- data.frame(COLLOCATE=rownames(pearson.residuals), pearson.residuals,
                           SUMABSDEV=apply(pearson.residuals, 1, \(af) sum(abs(af))),
                           LARGESTPREF=colnames(pearson.residuals)[apply(pearson.residuals, 1, \(af) which.max(af))])
#output.table <- output.table[order(output.table$LARGESTPREF, -output.table$SUMABSDEV),]
    #write.csv(output.table, "all_females_results.csv")
    #write.csv(input.matrix, "all_females_raw_counts.csv")
    


    library(reshape2)
    library(ggplot2)
    
    longData<-melt(output.table[ , -(ncol(output.table) - 1)])
    longData<-longData[longData$value!=0,]
    
    
 longData$r_val <- round(longData$value, 3)
 
 ind <- stringr::str_sub(female, 70,71) 
    
 plot(ggplot(longData, aes(x = variable, y = COLLOCATE)) + 
      geom_raster(aes(fill= value)) + 
      scale_fill_gradient(low="grey90", high="red") +
      labs(x="Call_2", y="Call_1", title="All_females") +
      geom_text(  aes(x=variable, y=COLLOCATE, label = r_val), color="black", size=5) + 
      theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                         axis.text.y=element_text(size=9),
                         plot.title=element_text(size=11)))
    
}
dev.off()