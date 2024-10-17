#this script reads the exifdata of a wave file and extracts annotations generated in Adobe Audition
#it requires a wav file that was annotation with section labels in Audition. For a different type of 
#annotation changes might be required

#as is, the script reads one, named file. 
#for better productivity this can be made into a loop to read multiple files from a folder and save all the data in one output table

library(exiftoolr)
library(strex)
 

#read exif information of one, named WAV file located in the working directory
#this step might take some time for larger files!!!
test <- exif_read(
  "file_12_92_(2009_06_04-13_41_18)_BSWMUW29297.wav",
  tags = NULL,
  recursive = FALSE,
  args = NULL,
  quiet = TRUE,
  pipeline = "csv"
)

 
#the "test" objects includes all exif data from the file
#now we are selecting only the needed information from it

#get all call start times (in frames)
St <- str_extract_numbers(
  test$TracksMarkersStartTime,
  decimals = FALSE,
  negs = FALSE,
  sci = FALSE,
  commas = T,
  leave_as_string = FALSE
)

#get all call duration (in frames)
dur <- str_extract_numbers(
  test$TracksMarkersDuration,
  decimals = FALSE,
  negs = FALSE,
  sci = FALSE,
  commas = T,
  leave_as_string = FALSE
)


#get all call names
names <- str_replace_all(test$TracksMarkersName, " ", "")


#get frame rate for calculating the time in seconds
fr <- str_extract_numbers(
  test$TracksFrameRate,
  decimals = FALSE,
  negs = FALSE,
  sci = FALSE,
  commas = T,
  leave_as_string = FALSE
)[[1]][1]


#combine all info together 
file_lables_form <- as.data.frame(cbind(test$SourceFile, #file names
                                        St[[1]]/fr,      #start time
                                        dur[[1]]/fr,     #duration
                                        str_split(names, ",", n = Inf, simplify = F)[[1]])) #call names
#name the columns of the output table
colnames(file_lables_form) <- c("file", "S_time", "duration", "call_type")  

##this can be saved as a separate csv##