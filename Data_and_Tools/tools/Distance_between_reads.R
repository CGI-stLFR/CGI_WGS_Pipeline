args = commandArgs(TRUE)
IN   = args[1]
ymax = args[2]
OUT  = args[3]

pdf(OUT,width=12,height=6)
par(bg='white')
y<-read.table(IN)
YMAX=as.numeric(ymax)
plot(y$V1,y$V3,type='l',col="blue",xlim=c(-100,100),ylim=c(0,YMAX),xlab="distance",ylab="percentage (%)",main="Distance between reads")
axis(1,at=seq(-100,100,10))
axis(2,at=seq(0,YMAX,1))

dev.off()

