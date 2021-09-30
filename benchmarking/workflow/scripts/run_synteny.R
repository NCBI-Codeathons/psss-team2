###### -- build base level gold standard data ---------------------------------

suppressMessages(library(DECIPHER))

###### -- ad hoc function, convert synteny tables to blast like matches -------

ConvertResults <- function(SyntenyObject) {
  
  if (all(dim(SYN) != c(2, 2))) {
    stop ("Function was written to work on a single pairwise comparison.")
  }
  
  # subset out hit table and block table
  HT <- SyntenyObject[[1, 2]]
  BT <- SyntenyObject[[2, 1]]
  
  o1 <- order(BT[, 1],
              BT[, 2])
  BT <- BT[o1, ]
  i1 <- BT[, 1L]
  i2 <- BT[, 2L]
  
  # loop through all blocks, give them a unique index to collapse on
  Count <- Pi1 <- Pi2 <- 0L
  UI <- vector(mode = "integer",
               length = length(o1))
  for (m1 in seq_along(o1)) {
    # if current i1 and current i2 do not match previous
    # iterate count
    # if they do, don't iterate count
    
    if (Pi1 == i1[m1] &
        Pi2 == i2[m1]) {
      # do nothing
    } else {
      # iterate count
      Count <- Count + 1L
    }
    
    # reassign Pi1 and Pi2 for next iteration
    Pi1 <- i1[m1]
    Pi2 <- i2[m1]
    
    # assign current count
    UI[m1] <- Count
  }
  
  SoH <- mapply(FUN = function(x, y) {
    # total nucleotide matches
    sum(unname(HT[x:y, 4L, drop = TRUE]))
  },
  x = BT[, 9],
  y = BT[, 10])
  
  # this could likely be more succinct, but it works ...
  TNM <- split(x = SoH,
               f = UI)
  TNM <- sapply(TNM,
                function(x) sum(x))
  QStart <- split(x = BT[, 5L, drop = TRUE],
                  f = UI)
  QStart <- sapply(QStart,
                   function(x) min(x))
  SStart <- split(x = BT[, 6L, drop = TRUE],
                  f = UI)
  SStart <- sapply(SStart,
                   function(x) min(x))
  QEnd <- split(x = BT[, 7L, drop = TRUE],
                f = UI)
  QEnd <- sapply(QEnd,
                 function(x) max(x))
  SEnd <- split(x = BT[, 8L, drop = TRUE],
                f = UI)
  SEnd <- sapply(SEnd,
                 function(x) max(x))
  C1 <- split(x = BT[, 1L, drop = TRUE],
              f = UI)
  C1 <- sapply(C1,
               function(x) unique(x))
  C2 <- split(x = BT[, 2L, drop = TRUE],
              f = UI)
  C2 <- sapply(C2,
               function(x) unique(x))

  
  InitialTable <- cbind(C1,
                        C2,
                        QStart,
                        QEnd,
                        SStart,
                        SEnd,
                        TNM)
  MatchSpace <- unname(apply(X = InitialTable,
                             MARGIN = 1L,
                             FUN = function(x) {
                               min((x[4] - x[3]) + 1L, (x[6] - x[5]) + 1L)
                             }))
  
  Pident <- TNM / MatchSpace
  
  CN1 <- names(SyntenyObject[[1, 1]])[C1]
  CN2 <- names(SyntenyObject[[2, 2]])[C2]
  
  Res <- data.frame("qseqid" = CN1,
                    "sseqid" = CN2,
                    "pident" = Pident,
                    "length" = MatchSpace,
                    "qstart" = InitialTable[, 3L],
                    "qend" = InitialTable[, 4L],
                    "sstart" = InitialTable[, 5L],
                    "send" = InitialTable[, 6L],
                    stringsAsFactors = FALSE)
  
  # return(list(UI,
  #             BT,
  #             HT,
  #             SoH,
  #             TNM,
  #             InitialTable,
  #             MatchSpace,
  #             Pident))
  return(Res)
}

###### -- code body -----------------------------------------------------------

TIMESTART <- Sys.time()

# deal with flexible input data
# ARGS <- commandArgs(trailingOnly = TRUE)
TargetFolder <- snakemake@params[["target_folder"]]
FocalSeqs <- snakemake@input[["query_fna_filtered"]]

# local testing
# TargetFolder <- "~/PSSS_CodeathonData"
# FocalSeqs <- "~/PSSS_CodeathonData/nmdc_mga04781_contigs.fna"

FileVector01 <- list.files(path = TargetFolder,
                           full.names = TRUE)

# check for focal seqs in target folder, if present remove them

ch1 <- strsplit(x = FocalSeqs,
                split = "/",
                fixed = TRUE)
ch1 <- ch1[[1]][length(ch1[[1]])]

ch2 <- strsplit(x = FileVector01,
                split = "/",
                fixed = TRUE)
ch2 <- sapply(ch2,
              function(x) x[length(x)],
              simplify = TRUE)

if (ch1 %in% ch2) {
  FileVector01 <- FileVector01[-which(ch2 == ch1)]
  ch2 <- ch2[-which(ch2 == ch1)]
}

# build DB as tempfile - will discard upon session closure
tmp1 <- tempfile()

focalseqs <- readDNAStringSet(FocalSeqs)
focalseqs <- focalseqs[width(focalseqs) > 500]

Seqs2DB(seqs = focalseqs,
        type = "XStringSet",
        dbFile = tmp1,
        identifier = "1")

SearchIDs <- as.character(seq(length(FileVector01)) + 1L)

for (m1 in seq_along(SearchIDs)) {
  
  currentseqs <- readDNAStringSet(FileVector01[m1])
  currentseqs <- currentseqs[width(currentseqs) > 500]
  
  Seqs2DB(seqs = currentseqs,
          type = "XStringSet",
          dbFile = tmp1,
          identifier = SearchIDs[m1])
}

# find synteny with defaults

cat("\nBeginning synteny map construction...\n")

CONN <- dbConnect(SQLite(),
                  tmp1)

Res <- vector(mode = "list",
              length = length(SearchIDs))

for (m1 in seq_along(Res)) {
  # for ( m1 in 2:length(Res)) {
  ID <- SearchIDs[m1]
  FILE <- ch2[m1]
  
  SYN <- FindSynteny(dbFile = CONN,
                     identifier = c("1", ID),
                     verbose = TRUE)
  
  # save(SYN,
  #      ID,
  #      FILE,
  #      file = paste0("~/PSSS_Standard_SYN_",
  #                    formatC(x = as.integer(ID),
  #                            flag = 0,
  #                            width = 2,
  #                            format = "d"),
  #                    ".RData"),
  #      compress = "xz")
  
  Res[[m1]] <- SYN
}

dbDisconnect(conn = CONN)

TIMEEND <- Sys.time()

cat("\nScript completed in:\n")
print(TIMEEND - TIMESTART)
cat("\nSaving data ...\n")

# save(Res,
#      ch2,
#      FileVector01,
#      file = "~/PSSS_Standard_SYN_ALL_v1.RData",
#      compress = "xz")

ResTable <- vector(mode = "list",
                   length = length(Res))

for (m1 in seq_along(Res)) {
  x <- ConvertResults(SyntenyObject = Res[[m1]])
  ResTable[[m1]] <- x
}

ResTable <- do.call(rbind,
                    ResTable)

# save(Res,
#      ResTable,
#      ch2,
#      FileVector01,
#      file = "~/PSSS_Standard_SYN_ALL_v2.RData",
#      compress = "xz")

ResTable <- ResTable[ResTable$pident > 0.95, ]

write.table(x = ResTable,
            file = snakemake@output[["synteny_outfile"]],
            append = FALSE,
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

cat("\nDone\n")








