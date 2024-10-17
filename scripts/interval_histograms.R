setwd("C:/Users/vdemartsev/My Cloud/HyraxPostDoc/Call_sequences")
library(tidyr)
library(ggplot2)

AVISOFT_extraction <- read.delim ("logfile.txt", header = F)

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
conv_table <- read.csv("conversion_table.csv")

for (x in 1:nrow(AVISOFT_extraction)) {
  if (is.na(AVISOFT_extraction$call_type[x])) {
    AVISOFT_extraction$call_type_upd[x] <- NA
  } else{
    AVISOFT_extraction$call_type_upd[x] <-
      conv_table[which(conv_table$old == AVISOFT_extraction$call_type[x]) , "new"]
  }
}
as.data.frame(table(AVISOFT_extraction$call_type_upd))

write.csv(AVISOFT_extraction, "all_calls_renamed.csv")

ggplot(data = AVISOFT_extraction , aes(x = brake_sec)) + xlim(0, 20)  + geom_histogram(bins = 100)


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

all_pairs <- all_pairs[which(all_pairs$brake <2), ]

pair_counts <- as.data.frame(table(all_pairs$pair))

frequent_call_pairs <- pair_counts[which(pair_counts$Freq > 199),]

frequent_pairs <- all_pairs[which(all_pairs$pair %in% frequent_call_pairs$Var1) ,]

ggplot(data = frequent_pairs , aes(x = brake)) + xlim(0, 20)  + geom_histogram(bins = 100) + 
  facet_wrap(~ pair, scales = "free")


pairs_only <- separate(frequent_pairs, pair, c("call_1", "call_2"), sep = "_") [ , 1:2]

write.table(pairs_only, "hyrax_call_pairs_only.txt", row.names=FALSE,sep="\t", quote = FALSE)

IDS <- unique(frequent_pairs$ID)

for (ID in IDS) {
  ID_select <- frequent_pairs[which(frequent_pairs$ID == ID) , ] 
  ID_pairs_only <- separate(ID_select, pair, c("call_1", "call_2"), sep = "_") [ , 1:2]
  write.table(ID_pairs_only, paste(getwd(), "/sep_females/",  ID, "_call_pairs_only.txt", sep = ""), row.names=FALSE,sep="\t", quote = FALSE)
  }

