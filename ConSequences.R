---
  title: "ConSequences"
  author: "Beniamin Abramczyk"
  date: "29 of April, 2023"
  output: fasta_and_csv_files
---

##installing required packages##
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("Biostrings")
library(Biostrings)

### Name changing function ###
##input: list of DNAStringSet objects, data.frame with sequence metadata##
##output: fasta files, list of DNAStringSet objects##
##This function changes names of sequences in all your molecular markers##
##It returns a list of DNAStringSets with names changed according to uploaded table##
##it also checks if your sequences have the right accession numbers according to the uploaded table and prints a report##
##Regardless of order check, the function changes names of your sequences##
##The function also checks if strain IDs are in your genbank record name##
##stab = summary table with sequence data, marnames = molecular markers names##
##seqnames = sequences names column name in your stab, strains = strain name collumn in your stab ##
##slist = list of DNAStringSets##
change.seq.names <- function(slist, stab, marnames, seqnames, strains){
  changed.names <- list() #defining the list for fasta objects#
  stab <- replace(stab, stab=="", "no_accession") ##replaces lack of data with text "no_accession"##
  #function which changes names of sequences of all molecular markers#
  for(i in 1:length(marnames)){ 
    stab1 <- subset(stab, stab[[marnames[i]]] != "no_accession") ##subsets your stab only to rename sequences with accession numbers##
    stab1 <- stab1[order(stab1[[marnames[i]]]),] ##sorts alphabetically by accession##
    names.of.seq <- stab1[[seqnames[1]]] ##creates a vector of names by extracting a column from your stab##
    accession_numbers <- stab1[[marnames[i]]] ##creates a vector of accession numbers in the correct order##
    strains_names <- stab1[[strains[1]]]
    to_change <- slist[[marnames[i]]] ##extracts a single DNAStringSet from your list##
    to_change <- to_change[order(names(to_change), decreasing = F),]
    ##correction check of your sequences##
    ##defines a vectors, which function below will fill with FALSE/TRUE values for each sequence##
    accession_correction <- c() ##defines a vector for accession correction check##
    name_correction <- c() ##defines a vector for strain name correction check##
    for(g in 1:length(names(to_change))){
      accession_correction[g] <- grepl(accession_numbers[g], names(to_change)[g]) ##fills the vector checking each sequence name separately for the approperiate accession number abundance##
      name_correction[g] <- grepl(strains_names[g], names(to_change)[g])
    }
    accession_correction <- as.character(accession_correction) ##changes the format of accession correction logical values for text##
    name_correction <- as.character(name_correction) ##changes the format of name correction logical values for text##
    mistakes_accession <- grep("FALSE", accession_correction) ##searches for FALSE values position in vector "accession_correction"##
    mistakes_name <- grep("FALSE", name_correction) ##searches for FALSE values position in vector "name_correction"##
    if (length(mistakes_accession) == 0){ ##prints information about possible mistakes and their placement in your fasta file if accessions in metadata don't match those in names of your sequences in a DNAStringSet##
      print(paste("Sequences accession numbers in", marnames[i], "correct!", sep = " "))
    } else {
      print(paste("WARNING: Incorrect accession id of sequence number", mistakes_accession, "in", marnames[i]))
    }
    if (length(mistakes_name) == 0){ ##prints information about possible mistakes and their placement in your fasta file if accessions in metadata don't match those in names of your sequences in a DNAStringSet##
      print(paste("Strain IDs in", marnames[i], "correct!", sep = " "))
    } else {
      print(paste("WARNING: Incorrect strain id of sequence number", mistakes_name, "in", marnames[i]))
    }
    print(paste("Changing names of sequences in", marnames[i], sep = " "))
    names(to_change) <- names.of.seq ##changes the names##
    changed.names[[i]] <- to_change ##writes a DNAStringSet into a list##
  }
  names(changed.names) <- marnames ##names DNAStringSets in a list appropriately##
  ##saving fasta files of your sequences with changed names in your working directory##
  for(j in 1:length(marnames)){ 
    file_to_save <- changed.names[[marnames[j]]] ##defines a file from a list to save##
    writeXStringSet(file_to_save,  paste(marnames[[j]], "names_changed" , "fasta", sep = ".")) ##saves the file##
    print(paste("Fasta file for", marnames[j], "saved", sep = " ")) ##informs you about which file it saved##
  }
  return(changed.names) ##returns a list of DNAStringSets##
}

###missing data filling function###
##input: list of DNAStringSet objects, data.frame with sequence metadata##
##output: fasta files, list of DNAStringSet objects##
###function that adds "+" signs for molecular markers sequences###
###that do not have accession numbers in the uploaded table (empty cells not NAs)###
###the pluses should later (after concatenation) be replaced within any text editor by find and replace with question marks##
##stab = summary table with sequence data, marnames = molecular markers names##
##seqnames = sequences names column name in your stab##
##slist = list of DNAStringSets##
###Explanation: it is impossible to add question marks to a DNAStringSet###

fill.with.pluses <- function(slist, stab, marnames, seqnames){
  ##adding pluses instead of missing sequences##
  added_pluses <- list() #defining the list for fasta objects#
  stab <- replace(stab, stab=="", "no_accession") ##replaces lack of data with text "no_accession"##
  for(i in 1:length(marnames)){
    stab1 <- subset(stab, stab[[marnames[i]]] != "no_accession") ##extracts only rows with existing accession number for a certain molecular marker
    stab1 <- stab1[order(stab1[[seqnames[1]]]),] ##sorts stab1 alphabetically##
    present <- stab1[[seqnames[1]]] ##creates a vector of names for sequences present within certain molecular marker sequences##
    to_add <- slist[[marnames[i]]] ##extracts a DNAStringSet from a list##
    to_add <- to_add[order(names(to_add), decreasing = F),] ##sorts stab1 DNAStringSet alphabetically##
    stab1 <- subset(stab, stab[[marnames[i]]] == "no_accession")
    if (length(stab1[[seqnames[1]]])==0 ){ ##if all sequences are present it just adds your DNAStringSet to your DNAStringSet list##
      added_pluses[[i]] <- to_add ##adds a DNAStringSets to a list of DNAStringSets##
      print(paste("Sequences for all strains present in", marnames[i], sep = " ")) ##prints an information, that no sequence should be added to this stringset##
    } else {
    not_present <- stab1[[seqnames[1]]] ##names for sequences not present within certain molecular marker sequences##
    adding_length <- length(to_add)[1] ##creates a vector with a number of sequences##
    plus <- c() ##defines a vector for plus signs##
    ##creating an object with an acurate number of "+" signs##
    for(h in 1:width(to_add)[1]){
      plus[1] <- paste(plus, "+", sep="") ##adds a plus to one "cell" in a vector##
    }
    ##adding plus signs to DNAStringSet##
    for(g in 1:length(not_present)){
      to_add[[adding_length + g]] <- plus ##adds "+" signs vector to DNAStringSet ##
    }
    names(to_add)[(length(present)+1) : length(c(present,not_present))] <- c(not_present) ##naming sequences 
    added_pluses[[i]] <- to_add ##adds a DNAStringSets to a list of DNAStringSets##
    print(paste("Added", length(not_present), "sequences in", marnames[i], sep=" ")) ##prints an information about how many sequences it added##
    }
  }
  names(added_pluses) <- marnames ##names each DNAStringSet in your list##
  ##saving your sequence files##
  for(j in 1:length(marnames)){ 
    file_to_save <- added_pluses[[marnames[j]]] ##defines a file from a list to save##
    writeXStringSet(file_to_save,  paste(marnames[[j]], "pluses.added" , "fasta", sep = ".")) ##saves the file##
    print(paste("Fasta file for", marnames[j], "with pluses added saved", sep = " ")) ##informs you about which file it saved##
  }
  return(added_pluses) ##returns a list of DNAStringSets##
}

##concatenating function##
##input: list of DNAStringSet objects, data.frame with sequence metadata##
##output: fasta files, list of DNAStringSet objects##
##This function concatenates sequences of different molecular markers##
##It returns a DNAStringSet of concatenated sequences and saves it in your working directory##
##It checks if the names of your sequences are correct (refering to a metadata table), if they're incorrect it prints a report##
##containing numbers of sequences for each molecular marker separately, for which it found differences##
##in such situation the function returns sorted DNAStringSet##
##stab = summary table with sequence metadata, marnames = molecular markers names##
##seqnames = sequences names column name in your stab##
##slist = list of DNAStringSets##
##It also saves a partition file needed for further analysis##

concatenate.sequences <- function(slist, stab, marnames, seqnames){
  sorted_sequences <- list() ##creates a list for sorted DNAStringSets##
  ##sorting all your DNAStringSets alphabetically##
  for (k in 1:length(marnames)){
    to_sort <- slist[[marnames[k]]] ##extracts a DNAStringSet of certain molecular marker##
    sorted <- to_sort[order(names(to_sort), decreasing = F),] #creates a sorted DNAStringSet#
    sorted_sequences[[k]] <- sorted ##saves sorted DNAStringSet in a list##
  }
  names(sorted_sequences) <- marnames ##adds names for each of your DNAStringSets##
  stab <- stab[order(stab[[seqnames[1]]]),] ##sorts your metadata table alphabetically by names of sequences##
  ##checking if all of your DNAStringSets have the same number of sequences##
  combinations <- combn(names(slist), 2) ##creates a matrix of combinations of your markers names##
  combinations <- as.data.frame(combinations) ##transforms the matrix into a data frame##
  errors_of_lengths <- c() ##creates a vector for errors of numbers of sequences##
  errors_of_lengths_report <- c() ##creates a vector for names of molecular markers of unmatching numbers of sequences##
  ##checking if numbers of sequences in files to concatenate are the same##
  for(w in 1:length(combinations)){
    to_check <- combinations[[w]] ##creates a vector of two names of markers to check in a round##
    to_check <- as.character(to_check) ##changes the format from as.factor##
    length_check_1 <- length(sorted_sequences[[to_check[1]]]) ##creates a vector of number of sequences in marker 1##
    length_check_2 <- length(sorted_sequences[[to_check[2]]]) ##creates a vector of number of sequences in marker 2##
    if(sum(length_check_1)==sum(length_check_2)){ 
      errors_of_lengths[w] <- c("TRUE") ##saves a value "TRUE" in a vector errors of length if number of sequences is the same##
    }else{
      errors_of_lengths[w] <- c("FALSE") ##saves a value "FALSE" in a vector errors of length if number of sequences is not the same##
      errors_of_lengths_report[2*w-1] <- to_check[1] ##saves names of molecular markers of unmatching sequences numbers##
      errors_of_lengths_report[2*w] <- to_check[2] ##saves names of molecular markers of unmatching sequences numbers##
    }
  }
  mistakes_lengths <- grep("FALSE", errors_of_lengths) ##searches for FALSES in a vector errors_of_lengths##
  if(length(mistakes_lengths)==0){
  ##checking if order of your sequences is correct##
  error_of_names_correction_report <- list() ##defines a list for report of names correction##
  errors_of_names <- list() ##defines a list for errors of names##
  for(n in 1:length(marnames)){
    names_to_check <- names(sorted_sequences[[marnames[n]]]) #creates a vector of names to check for certain molecular marker##
    reference <- stab[[seqnames[1]]] #creates a vector for names reference from your stab#
    names_correction <- c() ##defines a vector of values for incorrect number of sequences names##
    for(g in 1:length(names_to_check)){
      names_correction[g] <- grepl(reference[g], names_to_check[g])
    }
    names_correction <- as.character(names_correction) ##changes the format of name correction logical values for text##
    mistakes_names <- grep("FALSE", names_correction) ##searches for FALSE values position in vector "name_correction"##
    errors_of_names[[n]] <- mistakes_names ##saves the search result in a list##
    error_of_names_correction_report[[n]] <- mistakes_names ##saves the report for each molecular marker in a list##
  }
  if (length(unlist(errors_of_names))==0){
    print("Names of Your sequences correct!")
    print("Concatenating")
    ##concatenating the sequences##
    values_for_partitions <- c() ##defines a vector for values of sequences lengths##
    concatenated_sequences <- sorted_sequences[[marnames[1]]] ##creates a DNAStringSet of sequences of molecular marker 1##
    for(x in 2:length(marnames)){
      if(x==2){
        values_for_partitions[x-1] <- width(sorted_sequences[[marnames[x-1]]])[1]  ##saves value for a partition##
        values_for_partitions[x] <- width(sorted_sequences[[marnames[x]]])[1] + values_for_partitions[x-1]##saves value for a partition##
      }else{
      values_for_partitions[x-1] <- width(sorted_sequences[[marnames[x-1]]])[1] + values_for_partitions[x-2]##saves value for a partition##
      values_for_partitions[x] <- width(sorted_sequences[[marnames[x]]])[1] + values_for_partitions[x-1] ##saves value for a partition##
      }
      concatenated_sequences <- xscat(concatenated_sequences, sorted_sequences[[marnames[x]]]) ##concatenates the sequences##
    }
    names(concatenated_sequences) <- stab[[seqnames[1]]] ##names your sequences in concatenated sequences##
    writeXStringSet(concatenated_sequences, "concatenated.fasta", format = "fasta") ##saves a file in your working directory##
    print("Concatenated sequences saved in >>concatenated.fasta<< in Your working directory") ##prints an information about where your sequences were saved##
    data_type <- rep(c("DNA"), times = length(marnames)) ##creates a vector with information about data_type for partitions file##
    ##creating a vector with data of partitions lengths'##
    partitions_data <- c(paste(marnames[1], " = 1-", values_for_partitions[1], sep="")) ##creating a vector with data of first partition placement##
    for(y in 2:length(values_for_partitions)){ ##creates a vector of data of partitions placement##
      partitions_data[y] <- c(paste(marnames[y], " = ", values_for_partitions[y-1]+1, "-", values_for_partitions[y], sep=""))
    }
    partitions <- as.data.frame(cbind(data_type, partitions_data)) ##creates a data frame with data about partitions' placement##
    names(partitions) <- c() ##removes names of columns in the data frame##
    write.csv(partitions, "partitions.csv")
    print("Partitions file for Your data was saved in >>partitions.csv<< in your working directory")
    return(concatenated_sequences) ##returns a DNAStringSet of concatenated sequences##
  } else {
    names(error_of_names_correction_report) <- marnames ##names objects in your reports list##
    print("FATAL ERROR: Names of Your sequences are incorrect. Here's your report:") ##prints the message##
    print(error_of_names_correction_report) ##prints the error of names report##
    return(sorted_sequences) ##returns sorted sequences##
  }
  }else{
    print("FATAL ERROR: Numbers of sequences in your files to concatenate are uneven:") ##prints the message## 
    for(q in 1:(length(errors_of_lengths_report)/2)){ ##prints the report of which sequences' numbers are unmatching##
      print(paste("Number of sequences in", errors_of_lengths_report[2*q-1], "is different than in", errors_of_lengths_report[2*q], sep=" "))
    }
    return(sorted_sequences) ##returns sorted sequences##
  }
}

##### EXAMPLE OF DATA CURATION USING FUNCTIONS PRESENTED ABOVE WITH UMBELOPSIS MOLECULAR DATA ####

##sequence files should be downloaded using batchentrez from NCBI GenBank database and exported as fasta file##
##download the data using txt files that are named by molecular markers##

###Uploading data###

act1 <- readDNAStringSet("act1.fasta", format = "fasta") ##reads your sequence file 1##
COI <- readDNAStringSet("COI.fasta", format = "fasta") ##reads your sequence file 2##
ITS <- readDNAStringSet("ITS.fasta", format = "fasta") ##reads your sequence file 3##
LSU <- readDNAStringSet("LSU.fasta", format = "fasta") ##reads your sequence file 4##
mcm7 <- readDNAStringSet("mcm7.fasta", format = "fasta") ##reads your sequence file 5##
SSU <- readDNAStringSet("SSU.fasta", format = "fasta") ##reads your sequence file 6##
all_seq <- list(act1, COI, ITS, LSU, mcm7, SSU) ##combines stringsets into one list##
names(all_seq) <- c("act1", "COI", "ITS", "LSU", "mcm7", "SSU") ##gives names to objects in a list##
##your column with names of sequences to be changed for in metadata table must not contain any "."##
table <- as.data.frame(read.table("Umbelopsis.csv", header=TRUE, sep=';')) ##uploads a table of metadata##
##this part of code should be adjusted to your need - maybe in your version of R the data frame is uploaded automatically as.characters##
##this part of code changes format of columns in your metadata table to text##
table$species <- as.character(table$species)
table$strain_ID <- as.character(table$strain_ID)
table$ITS <- as.character(table$ITS)
table$LSU <- as.character(table$LSU)
table$act1 <- as.character(table$act1)
table$SSU <- as.character(table$SSU)
table$mcm7 <- as.character(table$mcm7)
table$COI <- as.character(table$COI)
table$name <- as.character(table$name)

##FUNCTION change.seq.names USE##
after_name_change <- change.seq.names(slist = all_seq, stab = table, marnames = c("act1", "COI", "ITS", "LSU", "mcm7", "SSU"), seqnames = c("name"), strains = c("strain_ID"))

##Now take your files with changed names to any server to allign and trim##
##remember about deleting the sequence length information that trimal prints##
##Come back here##

###data upload###

act1_trim <- readDNAStringSet("act1_trim.fasta", format = "fasta") ##reads your trimmed sequence file 1##
COI_trim <- readDNAStringSet("COI_trim.fasta", format = "fasta") ##reads your trimmed sequence file 2##
ITS_trim <- readDNAStringSet("ITS_trim.fasta", format = "fasta") ##reads your trimmed sequence file 3##
LSU_trim <- readDNAStringSet("LSU_trim.fasta", format = "fasta") ##reads your trimmed sequence file 4##
mcm7_trim <- readDNAStringSet("mcm7_trim.fasta", format = "fasta") ##reads your trimmed sequence file 5##
SSU_trim <- readDNAStringSet("SSU_trim.fasta", format = "fasta") ##reads your trimmed sequence file 6##
all_seq_trim <- list(act1_trim, COI_trim, ITS_trim, LSU_trim, mcm7_trim, SSU_trim) ##combines stringsets into one list##
names(all_seq_trim) <- c("act1", "COI", "ITS", "LSU", "mcm7", "SSU") ##gives names to objects in a list##


##FUNCTION fill.with.pluses USE##
Umbelopsis.added <- fill.with.pluses(slist = all_seq_trim, stab = table, marnames = c("act1", "COI", "ITS", "LSU", "mcm7", "SSU"), seqnames = c("name"))

##FUNCTION concatenate.sequences USE##
Umbelopsis_conc <- concatenate.sequences(slist = Umbelopsis.added, stab = table, marnames = c("act1", "COI", "ITS", "LSU", "mcm7", "SSU"), seqnames = c("name"))

##Open your concatenation file in notepad, nano, textedit or any other simple text editor an replace "+" with "?"##
##Open your partitions file in a simple text editor, delete the first row (there should be nothing there) and delete all quotation marks##
##Upload all your files to a server - they are ready for analysis with modeltest-ng and raxml-ng##


