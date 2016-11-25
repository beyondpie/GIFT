### get protein information for gift.
### szu
### 2016-11-25



f1_nm <- "prot_uniprot.txt"
f2_nm <- "protein_explanation_uniprot"
out_nm <- "protinfo_gift.txt"

fromdata <- read.table(f1_nm,sep="\n",header = FALSE,stringsAsFactors = FALSE)
todata <- read.table(f2_nm, sep ="\n", stringsAsFactors = FALSE)
fileConn <- file(out_nm)

result <- c()
for (ele in fromdata$V1){
  pos <- grep(ele, todata$V1, value = TRUE)
  if(length(pos) > 1){
    break
  }
  if(pos > 0){
    result <- c(result, todata$V1[pos])
  } else{
    result <- c(result, "NA")
  }
}
write(result, fileConn,sep="\n")
close(fileConn)
