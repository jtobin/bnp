require(ggplot2)
require(reshape)
require(dplyr)

d = read.csv('tmp.dat', header = F, colClasses = 'numeric')

melted = melt(d)

mixturePriorPlot =
  ggplot(melted, aes(value, fill = variable, colour = variable)) +
    geom_density(alpha = 0.2)

