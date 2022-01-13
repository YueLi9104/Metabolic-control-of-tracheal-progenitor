

# load ggplot2
library(ggplot2)

# read the csv file, the first column is reserved for row names.
pathway <- read.csv("DAVID using FlybaseID with generatio.csv")

pathway$pathway <- factor(pathway$pathway, levels = rev(pathway$pathway))

p <- ggplot(pathway, aes(GeneRatio, pathway))

p <- p + geom_point()
# point size is adjusted to a proper size corresponding to its count.
p <- p + geom_point(aes(size = Count))


#  Display in a 3-dims way
pbubble <- p + geom_point(aes(size = Count, color = PValue))

#  Set gradient colors corresponding to PValue.
pr <- pbubble + scale_color_gradient(low = "blue", high = "red")

pr <- pr + labs(color = expression(PValue), size = "Count", x = "GeneRatio")

# Draw bubble chart

pr + theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
)
ggsave("out.pdf")