source("thing", local=TRUE)

x <- 11:150
y <- Code[11:150]

dir.create('results/rotational_positioning')

pdf( "results/rotational_positioning/rotational_positioning.pdf", width=10, height=4)
plot(x, y, type="l", col = "red",
xlab="Distance to 5'end(nt)", ylab="AA/TT/AT fraction", 
main="rotational positioning")
abline(v=(seq(0,150,10)), col="lightgray", lty="dotted")
dev.off()
