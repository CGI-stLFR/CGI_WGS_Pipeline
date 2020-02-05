
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "fraglen.pdf"
}

fraglen_mincutoff_file=args[1]
fraglen_mincutoff_data=read.table(fraglen_mincutoff_file,header = T,colClasses = c(rep("NULL",5),"numeric","numeric"),sep="\t")

library(ggplot2)
library(scales)
pdf(args[2])
ggplot(data=fraglen_mincutoff_data)+
  geom_histogram(aes(x=length))+
  scale_y_continuous(labels = scales::comma)+
  scale_x_continuous(labels = scales::comma)+
  labs(title="fragment_length_distribution")+
  theme_bw()

ggplot(data=fraglen_mincutoff_data)+
  geom_histogram(aes(x=readcount))+
  scale_y_continuous(labels = scales::comma)+
  scale_x_continuous(labels = scales::comma)+
  labs(title="fragment_readcount_distribution")+
  theme_bw()


dev.off()



