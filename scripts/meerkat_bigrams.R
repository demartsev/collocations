
library(gtools)
library(stringr)
library("dplyr")
library(tidyr)

setwd("C:/Users/vdemartsev/Dropbox/meerkats_shared/data")
all_call_data <- read.delim("2021-05-17_ALL_conflicts_resolved.csv")

all_call_data <- all_call_data[which(all_call_data$pred_focalType == "F" & all_call_data$isCall == T) ,]


table(all_call_data$callType)

both_years <- all_call_data

#here set the vectors for unifying modifier and call names
hybrid <- c("hybrid|hyb")
sequence <- c("seq|sq")
move <- c("move|mov")
agression <- c("aggress|agress|chat|growl|grunt|bark")
alarm <- c("alarm|alrm|alert|ala")
lost <- c("lost|loc|lc")

both_years$callType <- str_replace_all(both_years$callType, hybrid, "hyb")
both_years$callType <- str_replace_all(both_years$callType, sequence, "sq")
both_years$callType <- str_replace_all(both_years$callType, move, "mo")
both_years$callType <- str_replace_all(both_years$callType, c("lead"), "ld")
both_years$callType <- str_replace_all(both_years$callType, c("social"), "soc")
both_years$callType <- str_replace_all(both_years$callType, c("ukn"), "unk")
both_years$callType <- str_replace_all(both_years$callType, agression, "agg")
both_years$callType <- str_replace_all(both_years$callType, alarm, "al")
both_years$callType <- str_replace_all(both_years$callType, lost, "lc")
both_years[which(both_years$callType == "s"), "callType"] <- "sn"  
both_years[which(both_years$callType == "c"), "callType"] <- "cc"  


both_years$callType <- str_replace(both_years$callType, fixed("s+"), "sn ") #renaming s as the first element
both_years$callType <- str_replace(both_years$callType, "\\+s$", " sn ")   #renaming s as the second element

### get hyb calls ###

both_years$hyb <- NA
both_years[grepl("fu|hyb", both_years$callType), "hyb"] <- "fu"
both_years[grepl("sq", both_years$callType), "hyb"] <- "sq"



### get call type elements ####
call_types <- c("agg|al|cc|ld|mo|sn|soc|lc|unk") #call types that we need
both_years <- cbind(both_years , str_extract_all(both_years$callType, call_types, simplify = TRUE)) #getting the call type elements 
both_years$`1` <- as.character(both_years$`1`) #removing factors
both_years$`2` <- as.character(both_years$`2`) #removing factors

#collecting the call type elements in an alphabetic order
both_years$final <- ifelse(both_years$`1` < both_years$`2`, paste(both_years$`1`, both_years$`2`), paste(both_years$`2`, both_years$`1`)) 
# keeping the original order for sequential calls
both_years[which(both_years$hyb == "sq") , "final"] <- paste(both_years[which(both_years$hyb == "sq") , "1"], both_years[which(both_years$hyb == "sq") , "2"])


#looking at the frequencies of the main call types as a self check
call_types <- c("cc", "soc", "al", "agg", "sn", "ld", "mo", "lc", "unk")


freq <- data.frame()
for (i in 1:length (call_types))
{
  freq[i,1] <- sum(str_count(rbind(both_years$`1`, both_years$`2`) , pattern = call_types[i]), na.rm = T)
  freq[i,2] <- call_types[i]
}
freq<- freq[order(freq$V1, decreasing = F),] ### if we decide to recode by the rarest category decreasing is to be set to TRUE


#get sample size plot for sanity check
par(mfrow=c(1,2))
bp <- barplot(freq$V1, names.arg = freq$V2, main = "Sample sizes per call type") 
text(bp, 0, freq$V1 ,cex=1,pos=3) 
#recode the call-types into main call categories. The order of recoding is frequency based 
# hybrid call_types are collapsed to the more frequent type. 

both_years$type_group <- NA  

for ( i in 1:nrow(freq))
{ type <- freq$V2[i]

both_years[which(grepl(type, both_years$final)), "type_group"] <- type 
}


#sample sizes per call category

freq <- data.frame()
for (i in 1:length (call_types))
{
  freq[i,1] <- sum(str_count(both_years$type_group , pattern = call_types[i]), na.rm = T)
  freq[i,2] <- call_types[i]
}
freq<- freq[order(freq$V1, decreasing = F),] 

bp <- barplot(freq$V1, names.arg = freq$V2, main = "Sample sizes per call category") #get sample size plot for sanity check
text(bp, 0, freq$V1 ,cex=1,pos=3) 


#cleaning
both_years <- subset(both_years, select=-c(`1`, `2`, `3`))
names(both_years)[names(both_years) == "final"] <- "stn_call_type"


## the columns containing call type info:
## hyb - indicator of fused or sequential calls
## stn_call_type - details the call type elements in each call
## type_group - only main calls, hybrids were collapsed into the parent type acording to overal call frequency 


both_years$year <- str_sub(both_years$date, 1, 4)
all_pairs <- data.frame()
for (year in unique(both_years$year)) {
  year_select <- both_years[which(both_years$year == year), ]
  
  for (ind in unique(year_select$ind)) {
    ind_select <- year_select[which(year_select$ind == ind) ,] 
    
    for (day in unique(ind_select$date)) {
      day_select <- ind_select[which(ind_select$date == day) ,]
      day_select$brake <- as.numeric(as.POSIXct(day_select$t0GPS_UTC ) - lag(as.POSIXct(day_select$tEndGPS_UTC)))
      pairs_per_day <-
        cbind(
          paste(day_select$type_group[-length(day_select$type_group)], "_",   day_select$type_group[-1], sep = ""),
          day_select$brake[-1],
          day_select$ind[1],
          day_select$wavFileName[1])
          pairs_per_day <- pairs_per_day[which(pairs_per_day[,2] < 2), ]
          
          all_pairs <- rbind(all_pairs, pairs_per_day)
    }
  }
}

colnames(all_pairs) <- c("pair", "brake", "ID", "file")

all_pairs <- all_pairs[-which(grepl("unk", all_pairs$pair)), ]

pair_counts <- as.data.frame(table(all_pairs$pair))

frequent_call_pairs <- pair_counts[which(pair_counts$Freq > 50),]

frequent_pairs <- all_pairs[which(all_pairs$pair %in% frequent_call_pairs$Var1) ,]

pairs_only <- separate(frequent_pairs, pair, c("call_1", "call_2"), sep = "_") [ , 1:2]

input.matrix <- pairs_only

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



plot(ggplot(longData, aes(x = variable, y = COLLOCATE)) + 
       geom_raster(aes(fill= value)) + 
       scale_fill_gradient(low="grey90", high="red") +
       labs(x="Call_2", y="Call_1", title="All_meerkats") +
       geom_text(  aes(x=variable, y=COLLOCATE, label = r_val), color="black", size=5) + 
       theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                          axis.text.y=element_text(size=9),
                          plot.title=element_text(size=11)))
