
TA <- function(nt) {
  if (nt == "A"){
    Result <- "T"
  } else if (nt == "T"){
    Result <- "A"
  } else if (nt == "G"){
    Result <- "C"
  } else if (nt == "C"){
    Result <- "G"
  }
  return(Result)
}
TG <- function(nt) {
  if (nt == "A"){
    Result <- "C"
  } else if (nt == "T"){
    Result <- "G"
  } else if (nt == "G"){
    Result <- "T"
  } else if (nt == "C"){
    Result <- "A"
  }
  return(Result)
}
OneBpToAll <- function(nt) {
  if (nt == "A"){
    Result <- c("C", "G", "T")
  } else if (nt == "T"){
    Result <- c("G", "C", "A")
  } else if (nt == "G"){
    Result <- c("T", "C", "A")
  } else if (nt == "C"){
    Result <- c("A","T","G")
  }
  return(Result)
}
SaveCSV <- function(dataframe, filename){
  write.table(dataframe, filename,
              na = "",
              row.names = FALSE,
              col.names = FALSE,
              append = FALSE,
              sep = ",")
}

