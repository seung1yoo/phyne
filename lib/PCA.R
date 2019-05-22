## "INPUT"
#DISTANCE        example_2       example_3       example_1       example_4
#example_2       0.000000        0.020977        1.060653        0.007050
#example_3       0.020977        0.000000        1.069947        0.022595
#example_1       1.060653        1.069947        0.000000        1.069183
#example_4       0.007050        0.022595        1.069183        0.000000

library(RColorBrewer)
data <- read.table([INPUT], sep='\t', header=T, row.names=1)
pc <- princomp(as.matrix(data), scores=TRUE)
png(file=[OUTPUT], height=1000, width=1000)
plot(pc$score[,1:2], type='n', xlab='PC 1', ylab='PC 2', main = 'PCA Plot')
text(pc$score[,1:2], labels=colnames(data), col=brewer.pal(length(colnames(data)),"Set1"))
dev.off()
