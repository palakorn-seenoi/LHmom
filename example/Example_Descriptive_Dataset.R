library(LHmom)
data(package = "LHmom")
ls("package:LHmom")

#install.packages("moments")
library(moments)

data("bangkok1")
bangkok1$rainfall
summary(bangkok1$rainfall)
sd(bangkok1$rainfall)
moments::skewness(bangkok1$rainfall)

#install.packages("e1071")
#library(e1071)
#e1071::skewness(bangkok1$rainfall)
#e1071::skewness(bangkok1$rainfall, type=1)

data("khonkaen")
khonkaen
khonkaen$rainfall
summary(khonkaen$rainfall)
sd(khonkaen$rainfall)
moments::skewness(khonkaen$rainfall)

data("sarakham")
sarakham
sarakham$temperature
summary(sarakham$temperature)
sd(sarakham$temperature)
moments::skewness(sarakham$temperature)
