library(ggplot2)
install.packages("Cairo",repos = "http://cran.us.r-project.org")
library(Cairo)

plot <- ggplot(mtcars, aes(mpg, wt)) +
  geom_point()
#ggsave("mtcars.png")

Cairo(90, 150, file="mttest.png", type="png", bg="white", res = 300, units = "mm")
plot
dev.off()

write.csv(mtcars,"mtcars.csv")