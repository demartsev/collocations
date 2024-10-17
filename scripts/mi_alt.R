install.packages(c("quanteda", "quanteda.textmodels", "quanteda.textstats"))

library(quanteda)
library(quanteda.textmodels)
library(quanteda.textstats)
library(tidyverse)


df <- read.delim("logfile.txt")
colnames(df)[7] <- "file"

all_bigrams <- data.frame()

for (one_file in unique(df$file)) {
  curr_file <- subset(df, file == one_file)
  curr_file$lag <- curr_file$X0.08 - lag(curr_file$X2.62)
  curr_file$text <- NA
  if (nrow(curr_file) < 2) {next}
  for (i in 2:nrow(curr_file)) {
    if (curr_file$lag[i] >= 0.5) {curr_file$text[i] <- paste(curr_file$GW8[i-1])
    }else{curr_file$text[i] <- paste(curr_file$GW8[i-1], curr_file$GW8[i])}
  }
  
  all_bigrams <- rbind(all_bigrams, curr_file)
}

# Assuming you have a data frame 'df' with a column 'text' containing your text data

# Assuming you have a data frame 'df' with a column 'text' containing your text data
corpus <- corpus(all_bigrams$text)

# Tokenization and preprocessing
tokens <- tokens(corpus, remove_punct = TRUE, remove_numbers = TRUE)
dfm <- dfm(tokens, dfm_remove = stopwords("english"))

# Custom function to calculate mutual information
calculate_mi <- function(dfm) {
  mi_matrix <- dfm %*% t(dfm)
  diag(mi_matrix) <- 0

  # Check for empty matrices
  if (nrow(mi_matrix) == 0 || ncol(mi_matrix) == 0) {
    stop("Empty matrix encountered")
  }

  # Ensure both matrices have the same dimensions
  if (nrow(mi_matrix) != ncol(mi_matrix)) {
    stop("Matrix dimensions are not compatible")
  }

  # Calculate row and column sums
  row_sums <- rowSums(dfm)
  col_sums <- colSums(dfm)

  # Handle zero row or column sums
  row_sums[row_sums == 0] <- 1
  col_sums[col_sums == 0] <- 1
  
  # Check for zero row or column sums
  if (any(row_sums == 0) || any(col_sums == 0)) {
    stop("Zero row or column sum encountered")
  }



  mi_matrix <- mi_matrix / (row_sums %o% col_sums)
  mi_matrix <- log2(mi_matrix)
  return(mi_matrix)
}

# Calculate mutual information
mi <- calculate_mi(dfm)



# Extract collocations
collocations <- textstat_collocations(dfm, mi = mi)

# Visualize collocations
ggplot(collocations, aes(x = MI, y = reorder(collocation, MI))) +
  geom_bar(stat = "identity") +
  labs(x = "Mutual Information", y = "Collocation") +
  ggtitle("Mutual Information Collocations")
