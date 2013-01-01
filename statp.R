source("Statistical_positioning",local=TRUE)

x <- -1000:1000
y <- Same_strand[1:2001]
z <- Oppo_strand[1:2001]

dir.create('results/statisitcal_positioning')

pdf( "results/statisitcal_positioning/statistical_position1.pdf", width=16, height=6)
plot(x, y, type="l", col = "red",
xlab="Tag start to start distance(nt)", ylab="reads coincidence number", ylim=c(100000,500000),
main="statistical positioning for same strand")
abline(v=seq(-1000,1000,250), col = "gray60", lty="dotted")
dev.off()

pdf( "results/statisitcal_positioning/statistical_position2.pdf", width=16, height=6)
plot(x, z, type="l", col = "red",
xlab="Tag start to start distance(nt)", ylab="reads coincidence number", 
main="statistical positioning for oppo strand")
abline(v=seq(-1000,1000,250), col = "gray60", lty="dotted")
dev.off()