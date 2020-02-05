#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "duplicate.pdf"
}

data <- read.table(args[1])
rate1 <- data[,2]
rate2 <- data[,3]
rate3 <- data[,4]
pdf(args[2])
par(bg='white')
hist(rate1,col="blue",xlim=c(0,1),breaks=20,main="distribution of barcode rate in duplicate reads",xlab="barcode rate")
#box()
#hist(rate2,col="blue",xlim=c(0,1),breaks=20,main="distribution of dupbar_in_reads_rate",xlab="dupbar rate")
#box()
#hist(rate3,col="blue",xlim=c(0,1),breaks= 20,main="distribution of identity reads rate in dupbar",xlab="un_snps rate")
#box()
dev.off()
