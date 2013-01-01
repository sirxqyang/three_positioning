require(graphics)
args <- commandArgs()
filename <- args[6]
Data <- read.table(filename, header=FALSE)
x <- Data[,1]
y <- Data[,2]

## will make the in pdf format
dir.create('results')
dir.create('results/translational_positioning')
#Path <- 'results/translational_positioning/'
#Name <- sub("_tss", "", filename)
#Suffix <- ".pdf"
#transp=paste (Path, Name, Suffix, sep = "", collapse = NULL) 
pdf('results/translational_positioning/translational_positioning.pdf', width=10, height=4)
plot(x, y, col='white', xlab='Distance to TSS (bp)', ylab='Average nucleosome density')
lines(spline(x, y), col=1)
dev.off()