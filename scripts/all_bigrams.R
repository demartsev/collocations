#This is a script for running an MSCA and MICA analysis on annotated call sequences in
# Meerkats, Hyraxes, Coaties, Hyenas and Sifaka


library(gtools)
library(stringr)
library(stringi) 
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(chron)
library(lubridate) 
library(sjPlot)
library(svMisc)
library(cowplot)
library(gridExtra)
library(grid)
library(scales)
setwd("/mnt/EAS_ind/vdemartsev/analysis/collocations/data/labels")

#species list 
species <- c("hyrax", "meerkat", "coati", "sifaka", "hyena")

#max inter-call interval to consider for combinations
ICI <- 0.5

#minimal number of events for a call pair to be included in the analysis
pair_cutoff <- 49

specie <- "coati"

for (specie in species) {
  
  if(specie == "hyena") { 
    ICI <- 4
    hyena_csvs <- list.files(paste(getwd(), "/", specie,  sep = ""), full.names = T, recursive = T)
    all_hyena_calls <- data.frame()
    for (csv in hyena_csvs) {
      curr <- read.delim(csv, sep = "\t", header = F)
      curr$file <- basename(csv)
      all_hyena_calls <- rbind(all_hyena_calls , curr)
    }
   #make all lower case
   all_hyena_calls$V3 <- tolower(all_hyena_calls$V3) 
   #remove cubs and non_foc and other weird enteries
   all_hyena_calls <-  all_hyena_calls[!grepl("cub|soa|eoa|non|skp|checkpoint", all_hyena_calls$V3), ]
   #remove modifiers
   all_hyena_calls$V3 <- str_replace_all(all_hyena_calls$V3, "foc|unf|fcoc|no|\\?", "")
   #trim wihite spaces 
   all_hyena_calls$V3 <- trimws(all_hyena_calls$V3, "both")
   #write.csv(table(all_hyena_calls$V3), "hyena_calls.csv")
   
   #load call conversion table
   hyena_calls_conv <- read.csv("hyena_conversion_table.csv")
   
   for (i in 1:nrow(all_hyena_calls)) {
     all_hyena_calls$new_name[i] <- hyena_calls_conv[which(hyena_calls_conv$old == all_hyena_calls$V3[i]), "new"]
     all_hyena_calls$isCall[i] <- hyena_calls_conv[which(hyena_calls_conv$old == all_hyena_calls$V3[i]), "call"] }
   
   hyena_calls_only <- all_hyena_calls[which(all_hyena_calls$isCall == "yes") ,]
   
   options(digits.secs=3)  
   options(scipen = 999) 
   colnames(hyena_calls_only) <- c("StartS", "DurationS","call_type", "file", "new_name", "isCall")
   hyena_calls_only$EndS <- hyena_calls_only$StartS + hyena_calls_only$DurationS
   
   #reformatting sequential calls into two row format
   
   seq <- unique(hyena_calls_conv[hyena_calls_conv$seq == "yes", "new"])
   for (row in which(hyena_calls_only$new_name %in% seq)) {
     
     #print(row)
     ext_row <- hyena_calls_only[row, ]
     hyena_calls_only$new_name[row] <- hyena_calls_conv[which(hyena_calls_conv$new == hyena_calls_only$new_name[row]) , "comp1"][1]
     hyena_calls_only$EndS[row] <- hyena_calls_only$StartS[row] + hyena_calls_only$DurationS[row]/2 - 0.02
     ext_row$new_name <- hyena_calls_conv[which(hyena_calls_conv$new == ext_row$new_name ) , "comp2"][1]
     ext_row$StartS <- hyena_calls_only$EndS[row]+0.04
     hyena_calls_only <- rbind(hyena_calls_only,  ext_row)
   }
   
   hyena_calls_fin <- data.frame() 
   for (file in unique(hyena_calls_only$file)){
     curr <- hyena_calls_only[which(hyena_calls_only$file == file) ,]
     curr <- curr[order(curr$StartS),]
     curr$lag <- curr$EndS - lag(curr$StartS)
     hyena_calls_fin <- rbind(hyena_calls_fin, curr)}
   
   pairs_only <- data.frame()
   
   for (i in which(!is.na(hyena_calls_fin$lag))) {
     if (hyena_calls_fin$lag[i] < ICI) {pair <- c(hyena_calls_fin$new_name[i-1], 
                                                  hyena_calls_fin$new_name[i], 
                                                  hyena_calls_fin$lag[i])
     pairs_only <- rbind(pairs_only, pair)}
   }
   
   #input.matrix <- read.table(file.choose(), header=TRUE, sep="\t", quote="", comment.char="")
   names(pairs_only) <- c("W_C", "Coll_Word")
   
   pairs_only$pairs <- paste(pairs_only$W_C, pairs_only$Coll_Word, sep = "-")
   
   pair_counts <- as.data.frame(table(pairs_only$pairs))
   
   frequent_call_pairs <- pair_counts[which(pair_counts$Freq > pair_cutoff),]
   
   frequent_pairs <- pairs_only[which(pairs_only$pairs %in% frequent_call_pairs$Var1) ,]
   
   pairs_only <- frequent_pairs[, -3]
   input.matrix <- table(pairs_only$Coll_Word, pairs_only$W_C)
   
   # computation
   pearson.residuals <- as.data.frame.matrix(chisq.test(input.matrix, correct=FALSE)$residuals)
   output.table <- data.frame(COLLOCATE=rownames(pearson.residuals), pearson.residuals,
                              SUMABSDEV=apply(pearson.residuals, 1, \(af) sum(abs(af))),
                              LARGESTPREF=colnames(pearson.residuals)[apply(pearson.residuals, 1, \(af) which.max(af))])
   
   longData<-melt(output.table[ , -(ncol(output.table) - 1)])
   longData<-longData[longData$value!=0,]
   longData$r_val <- round(longData$value, 2)
   longData$pair <- paste(longData$variable, "-", longData$COLLOCATE, sep = "")
   longData_1 <- longData[which(longData$pair %in% frequent_call_pairs$Var1), ]
   
   MDCA <- plot(ggplot(longData_1, aes(x = variable, y = COLLOCATE)) + 
                  geom_raster(aes(fill= value)) + 
                  scale_fill_gradient2(low="blue", high="red", mid = "grey90" ) +
                  labs(x="Call_2", y="Call_1", title=element_blank()) +
                  geom_text(  aes(x=variable, y=COLLOCATE, label = r_val), color="black", size=5) + 
                  theme_bw() + theme(axis.text.x=element_text(size=12, angle=0, vjust=0.3),
                                     axis.text.y=element_text(size=12),
                                     plot.title=element_text(size=11),
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),
                                     panel.border = element_blank(),
                                     legend.title = element_blank()))
   
   
   ###calculating MICA#### 
   #for frequent call pairs only using the formula
   # PMI(x, y) = log2(P(x, y))/P(x)*P(y)
   PMI_summary <- data_frame()
   #get one pair
   for (pair in frequent_call_pairs$Var1){
     x_is <- strsplit(pair, "-")[[1]][1]
     y_is <- strsplit(pair, "-")[[1]][2]
     P_x_y <- frequent_call_pairs[which(frequent_call_pairs$Var1 == pair), "Freq"] / nrow(hyena_calls_fin)
     P_x <- length(which(hyena_calls_fin$new_name == x_is)) / nrow(hyena_calls_fin)
     P_y <- length(which(hyena_calls_fin$new_name == y_is)) / nrow(hyena_calls_fin)
     PMI <- log2(P_x_y) / P_x * P_y
     combination <- c(x_is, y_is, PMI)
     PMI_summary <- rbind(PMI_summary, combination)}
   
   colnames(PMI_summary) <- c("call_1", "call_2", "PMI")
   PMI_summary$PMI <- as.numeric(PMI_summary$PMI)
   
   MICA <- ggplot(PMI_summary, aes( call_1, call_2, fill= PMI)) + 
     geom_raster() + 
     scale_fill_gradient(low="blue", high="grey90")+
     labs(x="Call_2", y="Call_1", title="MICA") +
     geom_text(aes(x=call_1, y=call_2, label = round(PMI, 2)), color="black", size=5) + 
     theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                        axis.text.y=element_text(size=9),
                        plot.title=element_text(size=11),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.title = element_blank())
   grid.arrange(MDCA, MICA, nrow = 1, top = "Hyena")
   
   #### Calculating Articulation Rate ####
   hyena_calls_fin <- subset(hyena_calls_fin, lag >= 0)  
       
   
   pairs_list <- c()
   for (i in nrow(pearson.residuals):1) {
     for (j in 1:ncol(pearson.residuals)) {
       pair_order <-
         paste(rownames(pearson.residuals)[i], colnames(pearson.residuals)[j], sep = "-")
       pairs_list <- c(pairs_list,pair_order)
     }
     
   }
   
   colnames(frequent_pairs) <- c("W_C", "Coll_Word", "brake", "pairs")
   frequent_pairs$brake <- as.numeric(frequent_pairs$brake)
   
   ggplot(data = frequent_pairs, aes(x = brake)) + geom_density() +
     facet_wrap(~pairs)
   
   
   plots_list <- list()
   for (one_pair in pairs_list) {
     
     plots_list[[one_pair]] <- ggplot(data = subset(frequent_pairs, pairs == one_pair), aes(x = brake)) + 
       geom_density() + #ylim(0,30) + 
       theme(axis.title.x = element_blank(), 
             axis.title.y = element_blank(),
             axis.text.y=element_blank(), 
             axis.ticks.y=element_blank())+ 
       scale_x_continuous(
         labels = label_number(accuracy = 0.1))
   }
   x.grob <- textGrob("gig  grn rum sql whp \n Call_2", 
                      gp=gpar(fontsize=10))
   y.grob <- textGrob(paste("gig" , "grn", "rum", "sql", "whp", 
                            sep = "            "), 
                      gp=gpar(fontsize=10), rot=90)
   grid <- grid.arrange(grobs =  plots_list, ncol = 5, bottom = x.grob, left = y.grob, 
                        top = grid::textGrob("Time_lag", x = 0, hjust = 0))
   #replace non_frequent call pairs will NULL
   names_to_plot <- which(pairs_list %in% frequent_call_pairs$Var1)
   
   
   grid <- grid.arrange(NULL, plots_list[[2]], NULL, NULL, plots_list[[5]],  
                        NULL, NULL, NULL, plots_list[[9]],  NULL, 
                        NULL, NULL,  plots_list[[13]],  NULL , NULL, 
                        NULL,  plots_list[[17]], NULL,NULL, plots_list[[20]], 
                        plots_list[[21]], NULL,   NULL,  NULL, NULL,
                        ncol = 5, bottom = x.grob, left = y.grob, 
                        top = grid::textGrob("Time_lag", x = 0, hjust = 0)) 
   
   grid.arrange(MDCA, grid, nrow = 1, top = "Hyena")  
   
  }
  
  if (specie == "meerkat") {
call_data <- read.csv(paste(getwd(), "/", specie, "/all_calls_sync_resolved_2023-09-10_cc_sn_filt_with_oor.csv", sep = ""))

 
call_data$stn_call_type <- str_replace(call_data$stn_call_type, "lc", "ld")
call_data$type_group <- str_replace(call_data$type_group, "lc", "ld")


## the columns containing call type info:
## hyb - indicator of fused or sequential calls
## stn_call_type - details the call type elements in each call
## type_group - only main calls, hybrids were collapsed into the parent type acording to overal call frequency 
options(digits.secs = 3)

call_data$year <- str_sub(call_data$date, 1, 4)
call_data <- call_data[which(call_data$isCall == 1),]

call_data$t0GPS_UTC <- as.POSIXct(call_data$t0GPS_UTC, tz = "UTC")
call_data$tMidGPS_UTC <- as.POSIXct(call_data$tMidGPS_UTC, tz = "UTC")
call_data$tendGPS_UTC <- as.POSIXct(call_data$tendGPS_UTC, tz = "UTC")

call_data$t0.numeric[1]
 
#reformatting sequential calls into two row format
for (row in which(call_data$hyb == "sq")) {
  #print(row)
  ext_row <- call_data[row, ]
  call_data$type_group[row] <- stri_extract_first(call_data[row, "stn_call_type"], regex="\\w+")
  call_data$tendGPS_UTC[row] <- as.POSIXct(call_data$tMidGPS_UTC[row], tz = "UTC")-0.05
  ext_row$type_group <- stri_extract_last(call_data[row, "stn_call_type"], regex="\\w+")
  ext_row$t0GPS_UTC <- as.POSIXct(ext_row$tMidGPS_UTC,  tz = "UTC")+0.05
  call_data <- rbind(call_data,  ext_row)
  }


table(call_data$hyb)

all_pairs <- data.frame()
for (year in unique(call_data$year)) {
  year_select <- call_data[which(call_data$year == year), ]
  
  for (ind in unique(year_select$ind)) {
    ind_select <- year_select[which(year_select$ind == ind) ,] 
    
    for (day in unique(ind_select$date)) {
      day_select <- ind_select[which(ind_select$date == day) ,]
      day_select <-  day_select[order(day_select$t0GPS_UTC),]
      day_select$brake <- as.POSIXct(day_select$t0GPS_UTC , tz = "UTC") - 
        lag(as.POSIXct(day_select$tendGPS_UTC,  tz = "UTC"))
      pairs_per_day <-
        cbind(
          paste(day_select$type_group[-length(day_select$type_group)], "_",   day_select$type_group[-1], sep = ""),
          paste(day_select$hyb[-length(day_select$hyb)], "_",   day_select$hyb[-1], sep = ""),
          day_select$brake[-1],
          day_select$ind[1],
          day_select$wavFileName[1])
      
          pairs_per_day <- as.data.frame(pairs_per_day)
          pairs_per_day$V3 <- as.numeric(pairs_per_day$V3)
          pairs_per_day <- pairs_per_day[which(pairs_per_day$V3 < ICI), ] #here we filter by minimum lag
          
          all_pairs <- rbind(all_pairs, pairs_per_day)
    }
  }
}


colnames(all_pairs) <- c("pair", "hyb", "brake", "ID", "file")




all_pairs <- all_pairs[-which(grepl("unk", all_pairs$pair)), ]
all_pairs <- all_pairs[-which(grepl("NA", all_pairs$pair)), ]

pair_counts <- as.data.frame(table(all_pairs$pair))

frequent_call_pairs <- pair_counts[which(pair_counts$Freq > pair_cutoff),]

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

 



longData<-melt(output.table[ , -(ncol(output.table) - 1)])
longData<-longData[longData$value!=0,]


longData$r_val <- round(longData$value, 1)
longData$pair <- paste(longData$variable, "_", longData$COLLOCATE, sep = "")
longData_1 <- longData[which(longData$pair %in% frequent_call_pairs$Var1), ]


MDCA <- plot(ggplot(longData_1, aes(x = variable, y = COLLOCATE)) + 
       geom_raster(aes(fill= value)) + 
       scale_fill_gradient2(low="blue", high="red", mid = "grey90" ) +
       labs(x="Call_2", y="Call_1", title=element_blank()) +
       geom_text(  aes(x=variable, y=COLLOCATE, label = r_val), color="black", size=5) + 
       theme_bw() + theme(axis.text.x=element_text(size=12, angle=0, vjust=0.3),
                          axis.text.y=element_text(size=12),
                          plot.title=element_text(size=11),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          legend.title = element_blank()))


      ###calculating MICA#### 
#for frequent call pairs only using the formula
# PMI(x, y) = log2(P(x, y))/P(x)*P(y)
PMI_summary <- data_frame()
#get one pair
for (pair in frequent_call_pairs$Var1){
  x_is <- strsplit(pair, "_")[[1]][1]
  y_is <- strsplit(pair, "_")[[1]][2]
  P_x_y <- frequent_call_pairs[which(frequent_call_pairs$Var1 == pair), "Freq"] / nrow(call_data)
  P_x <- length(which(call_data$type_group == x_is)) / nrow(call_data)
  P_y <- length(which(call_data$type_group == y_is)) / nrow(call_data)
  PMI <- log2(P_x_y) / P_x * P_y
  combination <- c(x_is, y_is, PMI)
  PMI_summary <- rbind(PMI_summary, combination)}

  colnames(PMI_summary) <- c("call_1", "call_2", "PMI")
  PMI_summary$PMI <- as.numeric(PMI_summary$PMI)
  
  MICA <- ggplot(PMI_summary, aes( call_1, call_2, fill= PMI)) + 
    geom_raster() + 
    scale_fill_gradient(low="blue", high="grey90")+
    labs(x="Call_2", y="Call_1", title="MICA") +
    geom_text(aes(x=call_1, y=call_2, label = round(PMI, 2)), color="black", size=5) + 
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       legend.title = element_blank())
  grid.arrange(MDCA, MICA, nrow = 1, top = "Meerkats")
  
  
  
  ggplot(data = all_pairs, aes(x = brake)) + 
    geom_density() + ylim(0,30) + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) + facet_wrap(~ pair)
  

  #### Calculating Articulation Rate ####
  all_pairs <- subset(all_pairs, brake >= 0)        
  
  pairs_list <- c()
  for (i in nrow(pearson.residuals):1) {
    for (j in 1:ncol(pearson.residuals)) {
      pair_order <-
        paste(rownames(pearson.residuals)[i], colnames(pearson.residuals)[j], sep = "_")
      pairs_list <- c(pairs_list,pair_order)
    }
    
  }
  


plots_list <- list()
for (one_pair in pairs_list) {
  
  plots_list[[one_pair]] <- ggplot(data = subset(all_pairs, pair == one_pair), aes(x = brake)) + 
    geom_density() + ylim(0,30) + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank())+ 
    scale_x_continuous(
      labels = label_number(accuracy = 0.1))
}
x.grob <- textGrob(paste("agg al cc mo sn soc\n Call_2", 
                         sep = "                           "), 
                   gp=gpar(fontsize=10))
y.grob <- textGrob(paste("agg", "al", "cc", "ld", "mo", "sn", "soc",
                         sep = "                     "), 
                   gp=gpar(fontsize=10), rot=90)
grid <- grid.arrange(grobs =  plots_list, ncol = 6, bottom = x.grob, left = y.grob, 
                     top = grid::textGrob("Time_lag", x = 0, hjust = 0))

#replace non_frequent call pairs will NULL
names_to_plot <- which(pairs_list %in% frequent_call_pairs$Var1)


grid <- grid.arrange(NULL, NULL, plots_list[[3]], NULL, plots_list[[5]], plots_list[[6]],
                     NULL, NULL, plots_list[[9]], plots_list[[10]],  plots_list[[11]], plots_list[[12]],
                     NULL, NULL, NULL, plots_list[[16]],  plots_list[[17]], NULL,
                     NULL, NULL, NULL, NULL, plots_list[[23]], NULL,
                     plots_list[[25]], NULL, plots_list[[27]], NULL, plots_list[[29]], plots_list[[30]], 
                     NULL, plots_list[[32]], NULL, NULL, NULL, NULL,
                     plots_list[[37]], NULL, plots_list[[39]], NULL, NULL, plots_list[[42]],
                     
                     ncol = 6, bottom = x.grob, left = y.grob, 
                     top = grid::textGrob("Time_lag", x = 0, hjust = 0))   
 
grid.arrange(MDCA, grid, nrow = 1, top = "Meerkats")

}

  if(specie == "hyrax") { 
    AVISOFT_extraction <- read.delim (paste(getwd(), "/", specie, "/logfile.txt", sep = ""), header = F)
    
    AVISOFT_extraction <-
      separate(
        data = AVISOFT_extraction,
        V7,
        into = c(
          'drive',
          'Users',
          "vdemartsev",
          "desktop",
          "vocals",
          "ID",
          "File"
        ),
        sep = "\\\\"
      )
    
    AVISOFT_extraction <- AVISOFT_extraction[,-c(7:11)]
    
    colnames(AVISOFT_extraction) <-
      c("#",
        "call_type",
        "brake_sec",
        "t_start",
        "t_end",
        "F0",
        "ID",
        "File")
    
    AVISOFT_extraction[which(AVISOFT_extraction$brake_sec == "_") ,  "brake_sec"] <-
      NA
    
    AVISOFT_extraction$brake_sec <-
      as.numeric(AVISOFT_extraction$brake_sec)
    
    multiples <-  extract_numeric(AVISOFT_extraction$call_type)
    
    mult_calls <- data.frame()
    
    for (x in 1:length(multiples)) {
      if (!is.na(multiples[x])) { mult_calls <- rbind( mult_calls,  AVISOFT_extraction[rep(x, multiples[x]-1) , ]) }
    }
    
    mult_calls$brake_sec <- 0.5
    AVISOFT_extraction <- rbind(AVISOFT_extraction, mult_calls)
    AVISOFT_extraction$call_type_old <- AVISOFT_extraction$call_type
    
    AVISOFT_extraction$call_type <-
      gsub('[0-9]+', '', AVISOFT_extraction$call_type)
    AVISOFT_extraction$call_type <-
      gsub(" ", NA, AVISOFT_extraction$call_type)
    
    #all_elements <- as.data.frame(table(AVISOFT_extraction$call_type))
    #write.csv(all_elements, "conversion_table.csv")
    conv_table <- read.csv(paste(getwd(), "/", specie, "/conversion_table.csv", sep = ""))
    
    for (x in 1:nrow(AVISOFT_extraction)) {
      if (is.na(AVISOFT_extraction$call_type[x])) {
        AVISOFT_extraction$call_type_upd[x] <- NA
      } else{
        AVISOFT_extraction$call_type_upd[x] <-
          conv_table[which(conv_table$old == AVISOFT_extraction$call_type[x]) , "new"]
      }
    }
    as.data.frame(table(AVISOFT_extraction$call_type_upd)) 
    
    
    
    
    files <- unique(AVISOFT_extraction$File)
    
    all_pairs <- data.frame()
    
    for (file in files) {
      file_select <-
        AVISOFT_extraction[which(AVISOFT_extraction$File == file) , ]
      
      file_select <-   file_select[order(file_select$t_start),]
      
      if (nrow(file_select) > 1) {
        pairs_per_file <-
          cbind(
            paste(file_select$call_type_upd[-length(file_select$call_type_upd)], "_",   file_select$call_type_upd[-1], sep = ""),
            file_select$brake_sec[-1],
            file_select$ID[1],
            file_select$File[1]
          )
        
        all_pairs <- rbind(all_pairs, pairs_per_file)
      }
      
      
    }
    
    
    colnames(all_pairs) <- c("pair", "brake", "ID", "file")
    #write.csv(all_pairs, "all_call_pairs_full_meta.csv")
    
    all_pairs$brake <- as.numeric(all_pairs$brake)
    
    all_pairs <- all_pairs[-which(grepl("NA", all_pairs$pair)), ]
    
    all_pairs <- all_pairs[which(all_pairs$brake < ICI), ]
    
    pair_counts <- as.data.frame(table(all_pairs$pair))
    
    frequent_call_pairs <- pair_counts[which(pair_counts$Freq > pair_cutoff),]
    
    frequent_pairs <- all_pairs[which(all_pairs$pair %in% frequent_call_pairs$Var1) ,]
    
    #ggplot(data = frequent_pairs , aes(x = brake)) + xlim(0, 20)  + geom_histogram(bins = 100) + 
    # facet_wrap(~ pair, scales = "free")
    
    
    pairs_only <- separate(frequent_pairs, pair, c("call_1", "call_2"), sep = "_") [ , 1:2]
    
    #input.matrix <- read.table(file.choose(), header=TRUE, sep="\t", quote="", comment.char="")
    names(pairs_only) <- c("W_C", "Coll_Word")
    input.matrix <- table(pairs_only$Coll_Word, pairs_only$W_C)
    
    # computation
    pearson.residuals <- as.data.frame.matrix(chisq.test(input.matrix, correct=FALSE)$residuals)
    output.table <- data.frame(COLLOCATE=rownames(pearson.residuals), pearson.residuals,
                               SUMABSDEV=apply(pearson.residuals, 1, \(af) sum(abs(af))),
                               LARGESTPREF=colnames(pearson.residuals)[apply(pearson.residuals, 1, \(af) which.max(af))])
    
    longData<-melt(output.table[ , -(ncol(output.table) - 1)])
    longData<-longData[longData$value!=0,]
    longData$r_val <- round(longData$value, 2)
    longData$pair <- paste(longData$variable, "_", longData$COLLOCATE, sep = "")
    longData_1 <- longData[which(longData$pair %in% frequent_call_pairs$Var1), ]
    #ind <- stringr::str_sub(female, 70,71) 
    
 MDCA <-  plot(ggplot(longData_1, aes(x = variable, y = COLLOCATE)) + 
           geom_raster(aes(fill= value)) + 
           scale_fill_gradient2(low="blue", high="red", mid = "grey90" ) +
           labs(x="Call_2", y="Call_1", title=element_blank()) +
           geom_text(  aes(x=variable, y=COLLOCATE, label = r_val), color="black", size=5) + 
           theme_bw() + theme(axis.text.x=element_text(size=12, angle=0, vjust=0.3),
                              axis.text.y=element_text(size=12),
                              plot.title=element_text(size=11),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              legend.title = element_blank()))
    
    ###calculating MICA#### 
    #for frequent call pairs only using the formula
    # PMI(x, y) = log2(P(x, y))/P(x)*P(y)
    PMI_summary <- data_frame()
    #get one pair
    for (pair in frequent_call_pairs$Var1){
      x_is <- strsplit(pair, "_")[[1]][1]
      y_is <- strsplit(pair, "_")[[1]][2]
      P_x_y <- frequent_call_pairs[which(frequent_call_pairs$Var1 == pair), "Freq"] / nrow(AVISOFT_extraction)
      P_x <- length(which(AVISOFT_extraction$call_type_upd == x_is)) / nrow(AVISOFT_extraction)
      P_y <- length(which(AVISOFT_extraction$call_type_upd == y_is)) / nrow(AVISOFT_extraction)
      PMI <- log2(P_x_y) / P_x * P_y
      combination <- c(x_is, y_is, PMI)
      PMI_summary <- rbind(PMI_summary, combination)}
 
      colnames(PMI_summary) <- c("call_1", "call_2", "PMI")
      PMI_summary$PMI <- as.numeric(PMI_summary$PMI)
      
MICA <- ggplot(PMI_summary, aes( call_1, call_2, fill= PMI)) + 
        geom_raster() + 
        scale_fill_gradient(low="blue", high="grey90")+
        labs(x="Call_2", y="Call_1", title="MICA") +
        geom_text(aes(x=call_1, y=call_2, label = round(PMI, 2)), color="black", size=5) + 
        theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                           axis.text.y=element_text(size=9),
                           plot.title=element_text(size=11),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           legend.title = element_blank())
      grid.arrange(MDCA, MICA, nrow = 1, top = "Hyraxes")
      
      
      #### Calculating Articulation Rate ####
      all_pairs <- subset(all_pairs, brake >= 0)        
      
      pairs_list <- c()
      for (i in nrow(pearson.residuals):1) {
        for (j in 1:ncol(pearson.residuals)) {
          pair_order <-
            paste(rownames(pearson.residuals)[i], colnames(pearson.residuals)[j], sep = "_")
          pairs_list <- c(pairs_list,pair_order)
        }
        
      }
      
      
      
      plots_list <- list()
      for (one_pair in pairs_list) {
        
        plots_list[[one_pair]] <- ggplot(data = subset(all_pairs, pair == one_pair), aes(x = brake)) + 
          geom_density() + ylim(0,8) + 
          theme(axis.title.x = element_blank(), 
                axis.title.y = element_blank(),
                axis.text.y=element_blank(), 
                axis.ticks.y=element_blank())+ 
          scale_x_continuous(
            labels = label_number(accuracy = 0.1))
      }
      x.grob <- textGrob(paste("Bark Click Growl Howl Twitter\n Call_2", 
                               sep = "                           "), 
                         gp=gpar(fontsize=10))
      y.grob <- textGrob(paste("Bark", "Click", "Growl",  "Howl", "Twitter",
                               sep = "                     "), 
                         gp=gpar(fontsize=10), rot=90)
      grid <- grid.arrange(grobs =  plots_list, ncol = 5, bottom = x.grob, left = y.grob, 
                           top = grid::textGrob("Time_lag", x = 0, hjust = 0))
      
      #replace non_frequent call pairs will NULL
      names_to_plot <- which(pairs_list %in% frequent_call_pairs$Var1)
      
      
      grid <- grid.arrange(NULL, plots_list[[2]], plots_list[[3]], NULL, plots_list[[5]], 
                           NULL, plots_list[[7]], plots_list[[8]], plots_list[[9]],  NULL, 
                           NULL, plots_list[[12]], plots_list[[13]], plots_list[[14]],  plots_list[[15]], 
                          
                           NULL, plots_list[[17]], plots_list[[18]], plots_list[[19]], plots_list[[20]], 
                           plots_list[[21]], NULL, NULL, NULL, NULL,

                           ncol = 5, bottom = x.grob, left = y.grob, 
                           top = grid::textGrob("Time_lag", x = 0, hjust = 0))    
      
      grid.arrange(MDCA, grid, nrow = 1, top = "Hyraxes")
      
      
    
    }

  if(specie == "coati") { 
    coati_csvs <- list.files(paste(getwd(), "/", specie, sep = ""), full.names = T)
    all_coati_calls <- data.frame()
    for (csv in coati_csvs) {
      curr <- read.csv(csv, sep = "\t")
      curr$file <- basename(csv)
      all_coati_calls <- rbind(all_coati_calls , curr)
      }
    
  all_coati_calls$Name <- trimws(all_coati_calls$Name, which =  "both")  
  
  #write.csv(table(all_coati_calls$Name), "coati_calls.csv")
  coati_calls_conv <- read.csv("coati_conversion_table.csv")
  
  #write.csv(setdiff(unique(all_coati_calls$Name), unique(coati_calls_conv$old)), "coati_diff_calls.csv")
  
   for (i in 1:nrow(all_coati_calls)) {
     all_coati_calls$new_name[i] <- coati_calls_conv[which(coati_calls_conv$old == all_coati_calls$Name[i]), "new"][1] 
     all_coati_calls$isCall[i] <- coati_calls_conv[which(coati_calls_conv$old == all_coati_calls$Name[i]), "call"] }
  
  coati_calls_only <- all_coati_calls[which(all_coati_calls$isCall == "yes") ,]
options(digits.secs=3)  
options(scipen = 999) 
  
  coati_calls_only$StartS <- difftime(parse_date_time(coati_calls_only$Start, orders = c("%H %M %OS", "%M %OS")), 
                                     parse_date_time("0", orders = "%0S"))
  coati_calls_only$DurationS <- difftime(parse_date_time(coati_calls_only$Duration, orders = c("%H %M %OS", "%M %OS")), 
                                      parse_date_time("0", orders = "%0S"))
  
  coati_calls_only$EndS <- coati_calls_only$StartS + coati_calls_only$DurationS
  
  #get below 1 sec calls only
  coati_calls_only <- coati_calls_only[which(coati_calls_only$DurationS < 1) ,]
  
  seq <- unique(coati_calls_conv[coati_calls_conv$seq == "yes", "new"])
  
  #reformatting sequential calls into two row format
  for (row in which(coati_calls_only$new_name %in% seq)) {
   
    #print(row)
    ext_row <- coati_calls_only[row, ]
    coati_calls_only$new_name[row] <- coati_calls_conv[which(coati_calls_conv$new == coati_calls_only$new_name[row]) , "comp1"]
    coati_calls_only$EndS[row] <- coati_calls_only$StartS[row] + coati_calls_only$DurationS[row]/2 - 0.02
    ext_row$new_name <- coati_calls_conv[which(coati_calls_conv$new == ext_row$new_name ) , "comp2"][1]
    ext_row$StartS <- coati_calls_only$EndS[row]+0.04
    coati_calls_only <- rbind(coati_calls_only,  ext_row)
  }
  
  
  
  
 coati_calls_fin <- data.frame() 
  for (file in unique(coati_calls_only$file)){
    curr <- coati_calls_only[which(coati_calls_only$file == file) ,]
    curr <- curr[order(curr$StartS),]
    curr$lag <- curr$EndS - lag(curr$StartS)
 coati_calls_fin <- rbind(coati_calls_fin, curr)}
 
 
pairs_only <- data.frame()
 
 for (i in which(!is.na(coati_calls_fin$lag))) {
   if (coati_calls_fin$lag[i] < ICI) {pair <- c(coati_calls_fin$new_name[i-1], 
                                                coati_calls_fin$new_name[i], 
                                                coati_calls_fin$lag[i])
   pairs_only <- rbind(pairs_only, pair)}
   }
 #input.matrix <- read.table(file.choose(), header=TRUE, sep="\t", quote="", comment.char="")
 names(pairs_only) <- c("W_C", "Coll_Word")
 
 pairs_only$pairs <- paste(pairs_only$W_C, pairs_only$Coll_Word, sep = "-")
 
 pair_counts <- as.data.frame(table(pairs_only$pairs))
 
 frequent_call_pairs <- pair_counts[which(pair_counts$Freq > pair_cutoff),]
 
 frequent_pairs <- pairs_only[which(pairs_only$pairs %in% frequent_call_pairs$Var1) ,]
 
 pairs_only <- frequent_pairs[, -3]
 input.matrix <- table(pairs_only$Coll_Word, pairs_only$W_C)
 
 # computation
 pearson.residuals <- as.data.frame.matrix(chisq.test(input.matrix, correct=FALSE)$residuals)
 output.table <- data.frame(COLLOCATE=rownames(pearson.residuals), pearson.residuals,
                            SUMABSDEV=apply(pearson.residuals, 1, \(af) sum(abs(af))),
                            LARGESTPREF=colnames(pearson.residuals)[apply(pearson.residuals, 1, \(af) which.max(af))])
 
 longData<-melt(output.table[ , -(ncol(output.table) - 1)])
 longData<-longData[longData$value!=0,]
 longData$r_val <- round(longData$value, 1)
 longData$pair <- paste(longData$variable, "-", longData$COLLOCATE, sep = "")
 longData_1 <- longData[which(longData$pair %in% frequent_call_pairs$Var1), ]
 longData_1$value <- round(longData_1$value, 1)
 longData_1$r_val <- round(longData_1$r_val, 1)
 
 MDCA <- plot(ggplot(longData_1, aes(x = variable, y = COLLOCATE)) + 
        geom_raster(aes(fill= value)) + 
          scale_fill_gradient2(low="blue", high="red", mid = "grey90" ) +
          labs(x="Call_2", y="Call_1", title=element_blank()) +
          geom_text(  aes(x=variable, y=COLLOCATE, label = r_val), color="black", size=5) + 
          theme_bw() + theme(axis.text.x=element_text(size=12, angle=0, vjust=0.3),
                             axis.text.y=element_text(size=12),
                             plot.title=element_text(size=11),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             legend.title = element_blank()))
 
 
 
 ###calculating MICA#### 
 #for frequent call pairs only using the formula
 # PMI(x, y) = log2(P(x, y))/P(x)*P(y)
 PMI_summary <- data_frame()
 #get one pair
 for (pair in frequent_call_pairs$Var1){
   x_is <- strsplit(pair, "-")[[1]][1]
   y_is <- strsplit(pair, "-")[[1]][2]
   P_x_y <- frequent_call_pairs[which(frequent_call_pairs$Var1 == pair), "Freq"] / nrow(coati_calls_fin)
   P_x <- length(which(coati_calls_fin$new_name == x_is)) / nrow(coati_calls_fin)
   P_y <- length(which(coati_calls_fin$new_name == y_is)) / nrow(coati_calls_fin)
   PMI <- log2(P_x_y) / P_x * P_y
   combination <- c(x_is, y_is, PMI)
   PMI_summary <- rbind(PMI_summary, combination)}
 
 colnames(PMI_summary) <- c("call_1", "call_2", "PMI")
 PMI_summary$PMI <- as.numeric(PMI_summary$PMI)
 
 MICA <- ggplot(PMI_summary, aes( call_1, call_2, fill= PMI)) + 
   geom_raster() + 
   scale_fill_gradient(low="blue", high="grey90")+
   labs(x="Call_2", y="Call_1", title="MICA") +
   geom_text(aes(x=call_1, y=call_2, label = round(PMI, 2)), color="black", size=5) + 
   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                      axis.text.y=element_text(size=9),
                      plot.title=element_text(size=11),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.title = element_blank())
 grid.arrange(MDCA, MICA, nrow = 1, top = "Coaties")
 
 
 
 #### Calculating Articulation Rate ####
 coati_calls_fin <- subset(coati_calls_fin, lag >= 0)  
 pairs_list <- c()
 for (i in nrow(pearson.residuals):1) {
   for (j in 1:ncol(pearson.residuals)) {
     pair_order <-
       paste(rownames(pearson.residuals)[i], colnames(pearson.residuals)[j], sep = "-")
     pairs_list <- c(pairs_list,pair_order)
   }
   
 }
 
 colnames(frequent_pairs) <- c("W_C", "Coll_Word", "brake", "pairs")
 frequent_pairs$brake <- as.numeric(frequent_pairs$brake)
 plots_list <- list()
 for (one_pair in pairs_list) {
   
   plots_list[[one_pair]] <- ggplot(data = subset(frequent_pairs, pairs == one_pair), aes(x = brake)) + 
     geom_density() + ylim(0,30) + 
     theme(axis.title.x = element_blank(), 
           axis.title.y = element_blank(),
           axis.text.y=element_blank(), 
           axis.ticks.y=element_blank())+ 
     scale_x_continuous(
       labels = label_number(accuracy = 0.1))
 }
 x.grob <- textGrob(paste("bark bop chirp chitter click dc grunt peep 
                          squeel\n Call_2", 
                          sep = "               "), 
                    gp=gpar(fontsize=10))
 y.grob <- textGrob(paste("bark", "bop", "chirp",  "chitter", "click", "dc", "grunt", "peep", 
                          "squeel",
                          sep = "            "), 
                    gp=gpar(fontsize=10), rot=90)
 grid <- grid.arrange(grobs =  plots_list, ncol = 9, bottom = x.grob, left = y.grob, 
                      top = grid::textGrob("Time_lag", x = 0, hjust = 0))
 #replace non_frequent call pairs will NULL
 names_to_plot <- which(pairs_list %in% frequent_call_pairs$Var1)
 
 
 grid <- grid.arrange(NULL, NULL, NULL, plots_list[[4]],  NULL, NULL, plots_list[[7]], NULL, plots_list[[9]],
                      NULL, NULL, NULL, plots_list[[13]],  NULL, NULL, NULL, plots_list[[17]], NULL,
                      NULL, NULL, plots_list[[21]], plots_list[[22]],  plots_list[[23]], NULL ,plots_list[[25]], NULL ,plots_list[[27]],
                      plots_list[[28]], NULL, plots_list[[30]], plots_list[[31]],  plots_list[[32]], plots_list[[33]],plots_list[[34]], NULL,NULL,
                      NULL, NULL, plots_list[[39]], NULL, plots_list[[41]],  NULL, plots_list[[43]], NULL,NULL, 
                      NULL,NULL, plots_list[[48]], plots_list[[49]], plots_list[[50]], plots_list[[51]],  plots_list[[52]], plots_list[[53]],plots_list[[54]],
                      NULL,NULL, plots_list[[57]], plots_list[[58]], plots_list[[59]], plots_list[[60]],  plots_list[[61]], NULL,plots_list[[63]], 
                      NULL, plots_list[[65]],NULL, NULL,  NULL, NULL,NULL, NULL,NULL,
                      plots_list[[73]], NULL, NULL, NULL,  NULL,  plots_list[[78]], plots_list[[79]], NULL,NULL,
                      ncol = 9, bottom = x.grob, left = y.grob, 
                      top = grid::textGrob("Time_lag", x = 0, hjust = 0)) 
 
 
 
 
 
 grid.arrange(MDCA, grid, nrow = 1, top = "Coaties")
 
 }
  
  if(specie == "sifaka") { 
    sifaka_csvs <- list.files(paste(getwd(), "/", specie, sep = ""), full.names = T)
    all_sifaka_calls <- data.frame()
    for (csv in sifaka_csvs) {
      curr <- read.csv(csv, sep = "\t")
      curr <- curr[,1:5 ]
      curr$file <- basename(csv)
      all_sifaka_calls <- rbind(all_sifaka_calls , curr)
    }
    
    all_sifaka_calls$Name <- trimws(all_sifaka_calls$Name, which =  "both")  
    
    #write.csv(table(all_sifaka_calls$Name), "sifaka_calls_1.csv")
    sifaka_calls_conv <- read.csv("sifaka_conversion_table.csv")
    
    
    #write.csv( setdiff(unique(all_sifaka_calls$Name), unique(sifaka_calls_conv$old)), "diff_calls.csv")
      
    for (i in 1:nrow(all_sifaka_calls)) {
      all_sifaka_calls$new_name[i] <- sifaka_calls_conv[which(sifaka_calls_conv$old == all_sifaka_calls$Name[i]), "new"] 
      all_sifaka_calls$isCall[i] <- sifaka_calls_conv[which(sifaka_calls_conv$old == all_sifaka_calls$Name[i]), "call"] }
    
    sifaka_calls_only <- all_sifaka_calls[which(all_sifaka_calls$isCall == "yes") ,]
    options(digits.secs=3)  
    options(scipen = 999) 
    
     
    
    sifaka_calls_only$StartS <- difftime(parse_date_time(sifaka_calls_only$Start, orders = c("%H %M %OS", "%M %OS")), 
                                        parse_date_time("0", orders = "%0S"), units="sec")
    sifaka_calls_only$DurationS <- difftime(parse_date_time(sifaka_calls_only$Duration, orders = c("%H %M %OS", "%M %OS")), 
                                           parse_date_time("0", orders = "%0S"))
    
    sifaka_calls_only$EndS <- sifaka_calls_only$StartS + sifaka_calls_only$DurationS
    
    
    
    #get below 10 sec calls only
    sifaka_calls_only <- sifaka_calls_only[which(sifaka_calls_only$DurationS < 10) ,]
    
    seq <- unique(sifaka_calls_conv[sifaka_calls_conv$seq == "yes", "new"])
    
    #reformatting sequential calls into two row format
    for (row in which(sifaka_calls_only$new_name %in% seq)) {
      #row <- 30
      #print(row)
      ext_row <- sifaka_calls_only[row, ]
      sifaka_calls_only$new_name[row] <- sifaka_calls_conv[which(sifaka_calls_conv$new == sifaka_calls_only$new_name[row]) , "comp1"][1]
      sifaka_calls_only$EndS[row] <- sifaka_calls_only$StartS[row] + sifaka_calls_only$DurationS[row]/2 - 0.02
      ext_row$new_name <- sifaka_calls_conv[which(sifaka_calls_conv$new == ext_row$new_name ) , "comp2"][1]
      ext_row$StartS <- sifaka_calls_only$EndS[row]+0.04
      sifaka_calls_only <- rbind(sifaka_calls_only,  ext_row)
    }
    
    
    
    
    sifaka_calls_fin <- data.frame() 
    for (file in unique(sifaka_calls_only$file)){
      curr <- sifaka_calls_only[which(sifaka_calls_only$file == file) ,]
      curr <- curr[order(curr$StartS),]
      curr$lag <- curr$EndS - lag(curr$StartS)
      sifaka_calls_fin <- rbind(sifaka_calls_fin, curr)}
    
    
    pairs_only <- data.frame()
    
    for (i in which(!is.na(sifaka_calls_fin$lag))) {
      if (sifaka_calls_fin$lag[i] < ICI) {pair <- c(sifaka_calls_fin$new_name[i-1], sifaka_calls_fin$new_name[i], sifaka_calls_fin$lag[i])
      pairs_only <- rbind(pairs_only, pair)}
    }
    #input.matrix <- read.table(file.choose(), header=TRUE, sep="\t", quote="", comment.char="")
    names(pairs_only) <- c("W_C", "Coll_Word", "brake")
    
    pairs_only$pairs <- paste(pairs_only$W_C, pairs_only$Coll_Word, sep = "-")
    
    pair_counts <- as.data.frame(table(pairs_only$pairs))
    
    frequent_call_pairs <- pair_counts[which(pair_counts$Freq > pair_cutoff),]
    
    frequent_pairs <- pairs_only[which(pairs_only$pairs %in% frequent_call_pairs$Var1) ,]
    
    pairs_only <- frequent_pairs[, -3]
    input.matrix <- table(pairs_only$Coll_Word, pairs_only$W_C)
    
    # computation
    pearson.residuals <- as.data.frame.matrix(chisq.test(input.matrix, correct=FALSE)$residuals)
    output.table <- data.frame(COLLOCATE=rownames(pearson.residuals), pearson.residuals,
                               SUMABSDEV=apply(pearson.residuals, 1, \(af) sum(abs(af))),
                               LARGESTPREF=colnames(pearson.residuals)[apply(pearson.residuals, 1, \(af) which.max(af))])
    
    longData<-melt(output.table[ , -(ncol(output.table) - 1)])
    longData<-longData[longData$value!=0,]
    longData$r_val <- round(longData$value, 2)
    longData$pair <- paste(longData$variable, "-", longData$COLLOCATE, sep = "")
    longData_1 <- longData[which(longData$pair %in% frequent_call_pairs$Var1), ]
    #longData$COLLOCATE <- as.factor(longData$COLLOCATE)
    
  MDCA <-  plot(ggplot(longData_1, aes(x = variable, y = COLLOCATE)) + 
           geom_raster(aes(fill= value)) +  
           scale_fill_gradient2(low="blue", high="red", mid = "grey90" ) +
           labs(x="Call_2", y="Call_1", title=element_blank()) +
           geom_text(  aes(x=variable, y=COLLOCATE, label = r_val), color="black", size=5) + 
           theme_bw() + theme(axis.text.x=element_text(size=12, angle=0, vjust=0.3),
                              axis.text.y=element_text(size=12),
                              plot.title=element_text(size=11),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              axis.line = element_blank(), 
                              legend.title = element_blank()))
  
  ###calculating MICA#### 
  #for frequent call pairs only using the formula
  # PMI(x, y) = log2(P(x, y))/P(x)*P(y)
  PMI_summary <- data_frame()
  #get one pair
  for (pair in frequent_call_pairs$Var1){
    x_is <- strsplit(pair, "-")[[1]][1]
    y_is <- strsplit(pair, "-")[[1]][2]
    P_x_y <- frequent_call_pairs[which(frequent_call_pairs$Var1 == pair), "Freq"] / nrow(sifaka_calls_fin)
    P_x <- length(which(sifaka_calls_fin$new_name == x_is)) / nrow(sifaka_calls_fin)
    P_y <- length(which(sifaka_calls_fin$new_name == y_is)) / nrow(sifaka_calls_fin)
    PMI <- log2(P_x_y) / P_x * P_y
    combination <- c(x_is, y_is, PMI)
    PMI_summary <- rbind(PMI_summary, combination)}
  
  colnames(PMI_summary) <- c("call_1", "call_2", "PMI")
  PMI_summary$PMI <- as.numeric(PMI_summary$PMI)
  
  MICA <- ggplot(PMI_summary, aes( call_1, call_2, fill= PMI)) + 
    geom_raster() + 
    scale_fill_gradient(low="blue", high="grey90")+
    labs(x="Call_2", y="Call_1", title="MICA") +
    geom_text(aes(x=call_1, y=call_2, label = round(PMI, 2)), color="black", size=5) + 
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       legend.title = element_blank())
  grid.arrange(MDCA, MICA, nrow = 1, top = "Sifakas")
  
  #### Calculating Articulation Rate ####
  sifaka_calls_fin <- subset(sifaka_calls_fin, lag >= 0)  
  pairs_list <- c()
  for (i in nrow(pearson.residuals):1) {
    for (j in 1:ncol(pearson.residuals)) {
      pair_order <-
        paste(rownames(pearson.residuals)[i], colnames(pearson.residuals)[j], sep = "-")
      pairs_list <- c(pairs_list,pair_order)
    }
    
  }
  
  
  
  
  
  
  #colnames(frequent_pairs) <- c("W_C", "Coll_Word", "brake", "pairs")
  frequent_pairs$brake <- as.numeric(frequent_pairs$brake)
  plots_list <- list()
  for (one_pair in pairs_list) {
    
    plots_list[[one_pair]] <- ggplot(data = subset(frequent_pairs, pairs == one_pair), aes(x = brake)) + 
      geom_density() + ylim(0,10) + 
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            axis.text.y=element_blank(), 
            axis.ticks.y=element_blank()) + 
      scale_x_continuous(
        labels = label_number(accuracy = 0.1))
  }
 
  
  x.grob <- textGrob("brk              chg               ctr               gof               grn               gwl\n Call_2", 
                     gp=gpar(fontsize=8))
  
  y.grob <- textGrob(paste("brk", "chg", "ctr",  "gof", "grn", "gwl",
                           sep = "            "), 
                     gp=gpar(fontsize=10), rot=90)
  
  

  
  #replace non_frequent call pairs will NULL
 names_to_plot <- which(pairs_list %in% frequent_call_pairs$Var1)
 
 
 grid <- grid.arrange(NULL, NULL, NULL, NULL,  NULL, plots_list[[6]], 
                      NULL, plots_list[[8]], NULL, NULL,  plots_list[[11]], NULL,
                      NULL, NULL, NULL, plots_list[[16]],  NULL, NULL,
                      NULL, NULL, plots_list[[21]], NULL,  NULL, NULL,
                      NULL, plots_list[[26]], NULL, NULL,  plots_list[[29]], NULL,
                      plots_list[[31]], NULL, NULL, NULL,  NULL, NULL,
                      ncol = 6, bottom = x.grob, #left = y.grob, 
                      top = grid::textGrob("Time_lag", x = 0, hjust = 0)) 
 

 
  grid.arrange(MDCA, grid, nrow = 1, top = "Sifaka", widths = c(1, 1))
  
  }
}
