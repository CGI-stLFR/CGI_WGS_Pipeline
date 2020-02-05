args = commandArgs(TRUE)
IN   = args[1]
xmax = args[2]
ymax = args[3]
OUT  = args[4]

pdf(OUT,width=12,height=6)
par(bg='white')
y<-read.table(IN)
XMAX=as.numeric(xmax)
YMAX=as.numeric(ymax)
plot(y$V1,y$V3,pch='.',col="blue",xlim=c(0,XMAX),ylim=c(0,YMAX),xlab="distance",ylab="percentage (%)",main="Distance between mate pairs reads")
#axis(1,at=seq(0,XMAX,10))
#axis(2,at=seq(0,YMAX,1))

dev.off()

