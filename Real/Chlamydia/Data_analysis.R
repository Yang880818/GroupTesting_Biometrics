## Trans xpt to csv
rm(list = ls())
library("foreign")
library("plyr")
## Chlamydia Data
Data_Chlmda_1 <- read.xport("./Data/1999_2000/LAB05.XPT")
Data_Chlmda_2 <- read.xport("./Data/2001_2002/L05_B.XPT")
Data_Chlmda_3 <- read.xport("./Data/2003_2004/L05_C.XPT")
Data_Chlmda_4 <- read.xport("./Data/2005_2006/CHLMDA_D.XPT")
Data_Chlmda_5 <- read.xport("./Data/2007_2008/CHLMDA_E.XPT")
Data_Chlmda_6 <- read.xport("./Data/2009_2010/CHLMDA_F.XPT")
Data_Chlmda_7 <- read.xport("./Data/2011_2012/CHLMDA_G.XPT")
Data_Chlmda_8 <- read.xport("./Data/2013_2014/CHLMDA_H.XPT")
Data_Chlmda_9 <- read.xport("./Data/2015_2016/CHLMDA_I.XPT")

write.csv(Data_Chlmda_1, file = "./Data/1999_2000/Data_Chlmda_1.csv")
write.csv(Data_Chlmda_2, file = "./Data/2001_2002/Data_Chlmda_2.csv")
write.csv(Data_Chlmda_3, file = "./Data/2003_2004/Data_Chlmda_3.csv")
write.csv(Data_Chlmda_4, file = "./Data/2005_2006/Data_Chlmda_4.csv")
write.csv(Data_Chlmda_5, file = "./Data/2007_2008/Data_Chlmda_5.csv")
write.csv(Data_Chlmda_6, file = "./Data/2009_2010/Data_Chlmda_6.csv")
write.csv(Data_Chlmda_7, file = "./Data/2011_2012/Data_Chlmda_7.csv")
write.csv(Data_Chlmda_8, file = "./Data/2013_2014/Data_Chlmda_8.csv")
write.csv(Data_Chlmda_9, file = "./Data/2015_2016/Data_Chlmda_9.csv")


Data_WholeBlood_1 <- read.xport("./Data/1999_2000/LAB25.XPT")
Data_WholeBlood_2 <- read.xport("./Data/2001_2002/L25_B.XPT")
Data_WholeBlood_3 <- read.xport("./Data/2003_2004/L25_C.XPT")
Data_WholeBlood_4 <- read.xport("./Data/2005_2006/CBC_D.XPT")
Data_WholeBlood_5 <- read.xport("./Data/2007_2008/CBC_E.XPT")
Data_WholeBlood_6 <- read.xport("./Data/2009_2010/CBC_F.XPT")
Data_WholeBlood_7 <- read.xport("./Data/2011_2012/CBC_G.XPT")
Data_WholeBlood_8 <- read.xport("./Data/2013_2014/CBC_H.XPT")
Data_WholeBlood_9 <- read.xport("./Data/2015_2016/CBC_I.XPT")

write.csv(Data_WholeBlood_1, file = "./Data/1999_2000/Data_WholeBlood_1.csv")
write.csv(Data_WholeBlood_2, file = "./Data/2001_2002/Data_WholeBlood_2.csv")
write.csv(Data_WholeBlood_3, file = "./Data/2003_2004/Data_WholeBlood_3.csv")
write.csv(Data_WholeBlood_4, file = "./Data/2005_2006/Data_WholeBlood_4.csv")
write.csv(Data_WholeBlood_5, file = "./Data/2007_2008/Data_WholeBlood_5.csv")
write.csv(Data_WholeBlood_6, file = "./Data/2009_2010/Data_WholeBlood_6.csv")
write.csv(Data_WholeBlood_7, file = "./Data/2011_2012/Data_WholeBlood_7.csv")
write.csv(Data_WholeBlood_8, file = "./Data/2013_2014/Data_WholeBlood_8.csv")
write.csv(Data_WholeBlood_9, file = "./Data/2015_2016/Data_WholeBlood_9.csv")


 

Data_Albumin_1 <- read.xport("./Data/1999_2000/LAB16.XPT")
Data_Albumin_2 <- read.xport("./Data/2001_2002/L16_B.XPT")
Data_Albumin_3 <- read.xport("./Data/2003_2004/L16_C.XPT")
Data_Albumin_4 <- read.xport("./Data/2005_2006/ALB_CR_D.XPT")
Data_Albumin_5 <- read.xport("./Data/2007_2008/ALB_CR_E.XPT")
Data_Albumin_6 <- read.xport("./Data/2009_2010/ALB_CR_F.XPT")
Data_Albumin_7 <- read.xport("./Data/2011_2012/ALB_CR_G.XPT")
Data_Albumin_8 <- read.xport("./Data/2013_2014/ALB_CR_H.XPT")
Data_Albumin_9 <- read.xport("./Data/2015_2016/ALB_CR_I.XPT")

write.csv(Data_Albumin_1, file = "./Data/1999_2000/Data_Albumin_1.csv")
write.csv(Data_Albumin_2, file = "./Data/2001_2002/Data_Albumin_2.csv")
write.csv(Data_Albumin_3, file = "./Data/2003_2004/Data_Albumin_3.csv")
write.csv(Data_Albumin_4, file = "./Data/2005_2006/Data_Albumin_4.csv")
write.csv(Data_Albumin_5, file = "./Data/2007_2008/Data_Albumin_5.csv")
write.csv(Data_Albumin_6, file = "./Data/2009_2010/Data_Albumin_6.csv")
write.csv(Data_Albumin_7, file = "./Data/2011_2012/Data_Albumin_7.csv")
write.csv(Data_Albumin_8, file = "./Data/2013_2014/Data_Albumin_8.csv")
write.csv(Data_Albumin_9, file = "./Data/2015_2016/Data_Albumin_9.csv")


## Demo Data
Data_Demo_1 <- read.xport("./Data/1999_2000/DEMO.XPT")
Data_Demo_2 <- read.xport("./Data/2001_2002/DEMO_B.XPT")
Data_Demo_3 <- read.xport("./Data/2003_2004/DEMO_C.XPT")
Data_Demo_4 <- read.xport("./Data/2005_2006/DEMO_D.XPT")
Data_Demo_5 <- read.xport("./Data/2007_2008/DEMO_E.XPT")
Data_Demo_6 <- read.xport("./Data/2009_2010/DEMO_F.XPT")
Data_Demo_7 <- read.xport("./Data/2011_2012/DEMO_G.XPT")
Data_Demo_8 <- read.xport("./Data/2013_2014/DEMO_H.XPT")
Data_Demo_9 <- read.xport("./Data/2015_2016/DEMO_I.XPT")

write.csv(Data_Demo_1, file = "./Data/1999_2000/Data_Demo_1.csv")
write.csv(Data_Demo_2, file = "./Data/2001_2002/Data_Demo_2.csv")
write.csv(Data_Demo_3, file = "./Data/2003_2004/Data_Demo_3.csv")
write.csv(Data_Demo_4, file = "./Data/2005_2006/Data_Demo_4.csv")
write.csv(Data_Demo_5, file = "./Data/2007_2008/Data_Demo_5.csv")
write.csv(Data_Demo_6, file = "./Data/2009_2010/Data_Demo_6.csv")
write.csv(Data_Demo_7, file = "./Data/2011_2012/Data_Demo_7.csv")
write.csv(Data_Demo_8, file = "./Data/2013_2014/Data_Demo_8.csv")
write.csv(Data_Demo_9, file = "./Data/2015_2016/Data_Demo_9.csv")



### Filter and Rearrange
## 1999_2000
Chlmda_1 <- read.csv("./Data/1999_2000/Data_Chlmda_1.csv", header = T)

## Biomarkers 
WholeBlood_1 <- read.csv("./Data/1999_2000/Data_WholeBlood_1.csv", header = T)
Albumin_1 <- read.csv("./Data/1999_2000/Data_Albumin_1.csv", header = T)
Chlmda_Bios_1 <- na.omit(merge(Chlmda_1[c("SEQN", "URXUCL")], WholeBlood_1[c("SEQN", "LBDMONO", "LBDNENO")]))
dim(Chlmda_Bios_1)
Chlmda_Bios_1 <- Chlmda_Bios_1[which(Chlmda_Bios_1$LBDNENO > 0 & Chlmda_Bios_1$LBDMONO > 0), ]
dim(Chlmda_Bios_1)
Chlmda_Bios_1 <- na.omit(merge(Chlmda_Bios_1, Albumin_1[c("SEQN", "URXUMA")]))
dim(Chlmda_Bios_1)
Chlmda_Bios_1$URXUCL[which(Chlmda_Bios_1$URXUCL == 2)] <- 0 ### Replace 2 by 0 in "URXUCL"
write.csv(Chlmda_Bios_1, file = "./Data/1999_2000/Chlmda_Bios_1.csv")

## Weight
Demo_1 <- read.csv("./Data/1999_2000/Data_Demo_1.csv", header = T)
Chlmda_Bios_Demo_1_final <- na.omit(merge(Chlmda_Bios_1, Demo_1[c("SEQN", "WTINT2YR")]))
write.csv(Chlmda_Bios_Demo_1_final, file = "./Data/1999_2000/Chlmda_Bios_Demo_1_final.csv")
dim(Chlmda_Bios_Demo_1_final)



## 2001_2002
Chlmda_2 <- read.csv("./Data/2001_2002/Data_Chlmda_2.csv", header = T)

## Biomarkers
WholeBlood_2 <- read.csv("./Data/2001_2002/Data_WholeBlood_2.csv", header = T)
Albumin_2 <- read.csv("./Data/2001_2002/Data_Albumin_2.csv", header = T)
Chlmda_Bios_2 <- na.omit(merge(Chlmda_2[c("SEQN", "URXUCL")], WholeBlood_2[c("SEQN", "LBDMONO", "LBDNENO")]))
dim(Chlmda_Bios_2)
Chlmda_Bios_2 <- Chlmda_Bios_2[which(Chlmda_Bios_2$LBDNENO > 0 & Chlmda_Bios_2$LBDMONO > 0), ]
dim(Chlmda_Bios_2)
Chlmda_Bios_2 <- na.omit(merge(Chlmda_Bios_2, Albumin_2[c("SEQN", "URXUMA")]))
dim(Chlmda_Bios_2)
Chlmda_Bios_2$URXUCL[which(Chlmda_Bios_2$URXUCL == 2)] <- 0 ### Replace 2 by 0 in "URXUCL"
write.csv(Chlmda_Bios_2, file = "./Data/2001_2002/Chlmda_Bios_2.csv")

## Weight
Demo_2 <- read.csv("./Data/2001_2002/Data_Demo_2.csv", header = T)
Chlmda_Bios_Demo_2_final <- na.omit(merge(Chlmda_Bios_2, Demo_2[c("SEQN", "WTINT2YR")]))
write.csv(Chlmda_Bios_Demo_2_final, file = "./Data/2001_2002/Chlmda_Bios_Demo_2_final.csv")
dim(Chlmda_Bios_Demo_2_final)



## 2003_2004
Chlmda_3 <- read.csv("./Data/2003_2004/Data_Chlmda_3.csv", header = T)

## Biomarkers
WholeBlood_3 <- read.csv("./Data/2003_2004/Data_WholeBlood_3.csv", header = T)
Albumin_3 <- read.csv("./Data/2003_2004/Data_Albumin_3.csv", header = T)
Chlmda_Bios_3 <- na.omit(merge(Chlmda_3[c("SEQN", "URXUCL")], WholeBlood_3[c("SEQN", "LBDMONO", "LBDNENO")]))
dim(Chlmda_Bios_3)
Chlmda_Bios_3 <- Chlmda_Bios_3[which(Chlmda_Bios_3$LBDNENO > 0 & Chlmda_Bios_3$LBDMONO > 0), ]
dim(Chlmda_Bios_3)
Chlmda_Bios_3 <- na.omit(merge(Chlmda_Bios_3, Albumin_3[c("SEQN", "URXUMA")]))
dim(Chlmda_Bios_3)
Chlmda_Bios_3$URXUCL[which(Chlmda_Bios_3$URXUCL == 2)] <- 0 ### Replace 2 by 0 in "URXUCL"
write.csv(Chlmda_Bios_3, file = "./Data/2003_2004/Chlmda_Bios_3.csv")

## Weight
Demo_3   <- read.csv("./Data/2003_2004/Data_Demo_3.csv", header = T)
Chlmda_Bios_Demo_3_final <- na.omit(merge(Chlmda_Bios_3, Demo_3[c("SEQN", "WTINT2YR")]))
write.csv(Chlmda_Bios_Demo_3_final, file = "./Data/2003_2004/Chlmda_Bios_Demo_3_final.csv")
dim(Chlmda_Bios_Demo_3_final)



## 2005_2006
Chlmda_4 <- read.csv("./Data/2005_2006/Data_Chlmda_4.csv", header = T)

## Biomarkers
WholeBlood_4 <- read.csv("./Data/2005_2006/Data_WholeBlood_4.csv", header = T)
Albumin_4 <- read.csv("./Data/2005_2006/Data_Albumin_4.csv", header = T)
Chlmda_Bios_4 <- na.omit(merge(Chlmda_4[c("SEQN", "URXUCL")], WholeBlood_4[c("SEQN", "LBDMONO", "LBDNENO")]))
dim(Chlmda_Bios_4)
Chlmda_Bios_4 <- Chlmda_Bios_4[which(Chlmda_Bios_4$LBDNENO > 0 & Chlmda_Bios_4$LBDMONO > 0), ]
dim(Chlmda_Bios_4)
Chlmda_Bios_4 <- na.omit(merge(Chlmda_Bios_4, Albumin_4[c("SEQN", "URXUMA")]))
dim(Chlmda_Bios_4)
Chlmda_Bios_4$URXUCL[which(Chlmda_Bios_4$URXUCL == 2)] <- 0 ### Replace 2 by 0 in "URXUCL"
write.csv(Chlmda_Bios_4, file = "./Data/2005_2006/Chlmda_Bios_4.csv")

## Weight
Demo_4   <- read.csv("./Data/2005_2006/Data_Demo_4.csv", header = T)
Chlmda_Bios_Demo_4_final <- na.omit(merge(Chlmda_Bios_4, Demo_4[c("SEQN", "WTINT2YR")]))
write.csv(Chlmda_Bios_Demo_4_final, file = "./Data/2005_2006/Chlmda_Bios_Demo_4_final.csv")
dim(Chlmda_Bios_Demo_4_final)



## 2007_2008
Chlmda_5 <- read.csv("./Data/2007_2008/Data_Chlmda_5.csv", header = T)

## Biomarkers
WholeBlood_5 <- read.csv("./Data/2007_2008/Data_WholeBlood_5.csv", header = T)
Albumin_5 <- read.csv("./Data/2007_2008/Data_Albumin_5.csv", header = T)
Chlmda_Bios_5 <- na.omit(merge(Chlmda_5[c("SEQN", "URXUCL")], WholeBlood_5[c("SEQN", "LBDMONO", "LBDNENO")]))
dim(Chlmda_Bios_5)
Chlmda_Bios_5 <- Chlmda_Bios_5[which(Chlmda_Bios_5$LBDNENO > 0 & Chlmda_Bios_5$LBDMONO > 0), ]
dim(Chlmda_Bios_5)
Chlmda_Bios_5 <- na.omit(merge(Chlmda_Bios_5, Albumin_5[c("SEQN", "URXUMA")]))
dim(Chlmda_Bios_5)
Chlmda_Bios_5$URXUCL[which(Chlmda_Bios_5$URXUCL == 2)] <- 0 ### Replace 2 by 0 in "URXUCL"
write.csv(Chlmda_Bios_5, file = "./Data/2007_2008/Chlmda_Bios_5.csv")

## Weight
Demo_5   <- read.csv("./Data/2007_2008/Data_Demo_5.csv", header = T)
Chlmda_Bios_Demo_5_final <- na.omit(merge(Chlmda_Bios_5, Demo_5[c("SEQN", "WTINT2YR")]))
write.csv(Chlmda_Bios_Demo_5_final, file = "./Data/2007_2008/Chlmda_Bios_Demo_5_final.csv")
dim(Chlmda_Bios_Demo_5_final)



## 2009_2010
Chlmda_6 <- read.csv("./Data/2009_2010/Data_Chlmda_6.csv", header = T)

## Biomarkers
WholeBlood_6 <- read.csv("./Data/2009_2010/Data_WholeBlood_6.csv", header = T)
Albumin_6 <- read.csv("./Data/2009_2010/Data_Albumin_6.csv", header = T)
Chlmda_Bios_6 <- na.omit(merge(Chlmda_6[c("SEQN", "URXUCL")], WholeBlood_6[c("SEQN", "LBDMONO", "LBDNENO")]))
dim(Chlmda_Bios_6)
Chlmda_Bios_6 <- Chlmda_Bios_6[which(Chlmda_Bios_6$LBDNENO > 0 & Chlmda_Bios_6$LBDMONO > 0), ]
dim(Chlmda_Bios_6)
Chlmda_Bios_6 <- na.omit(merge(Chlmda_Bios_6, Albumin_6[c("SEQN", "URXUMA")]))
dim(Chlmda_Bios_6)
Chlmda_Bios_6$URXUCL[which(Chlmda_Bios_6$URXUCL == 2)] <- 0 ### Replace 2 by 0 in "URXUCL"
write.csv(Chlmda_Bios_6, file = "./Data/2009_2010/Chlmda_Bios_6.csv")

## Weight
Demo_6   <- read.csv("./Data/2009_2010/Data_Demo_6.csv", header = T)
Chlmda_Bios_Demo_6_final <- na.omit(merge(Chlmda_Bios_6, Demo_6[c("SEQN", "WTINT2YR")]))
write.csv(Chlmda_Bios_Demo_6_final, file = "./Data/2009_2010/Chlmda_Bios_Demo_6_final.csv")
dim(Chlmda_Bios_Demo_6_final)






(dim(Chlmda_Bios_Demo_1_final)[1] + dim(Chlmda_Bios_Demo_2_final)[1] + dim(Chlmda_Bios_Demo_3_final)[1] 
  + dim(Chlmda_Bios_Demo_4_final)[1] + dim(Chlmda_Bios_Demo_5_final)[1] + dim(Chlmda_Bios_Demo_6_final)[1])


set.seed(1)
## 1999-2000
Y_1.weight <- Chlmda_Bios_Demo_1_final$WTINT2YR

### Calculate the true prevalence
prevalTrue_1 <- sum(Y_1.weight[Chlmda_Bios_Demo_1_final$URXUCL == 1])/sum(Y_1.weight)
print(paste("The prevalence in 1999-2000 data is", round(prevalTrue_1, 5)))

WeightStand_1 <- Y_1.weight/sum(Y_1.weight)
#WeightStand_1

Y_1.ID       <- Chlmda_Bios_Demo_1_final$SEQN
Y_1.disease  <- Chlmda_Bios_Demo_1_final$URXUCL
X_1.marker_1 <- Chlmda_Bios_Demo_1_final$LBDMONO
X_1.marker_2 <- Chlmda_Bios_Demo_1_final$LBDNENO
X_1.marker_3 <- Chlmda_Bios_Demo_1_final$URXUMA

Y_1.data <- cbind(Y_1.ID, Y_1.disease, X_1.marker_1, X_1.marker_2, X_1.marker_3, Y_1.weight)
colnames(Y_1.data) <- c("SEQN", "Disease_status", "Biomarker_MONO", "Biomarker_EOSI", "Biomarker_Albumin", "Weight")
numSam_1 <- nrow(Y_1.data)
resamLOC_1 <- sample(1:numSam_1, size = numSam_1, replace = T, prob = WeightStand_1)
Data_resample_1 <- Y_1.data[resamLOC_1,]
#Data_resample_1

## 2001-2002
Y_2.weight <- Chlmda_Bios_Demo_2_final$WTINT2YR

### Calculate the true prevalence
prevalTrue_2 <- sum(Y_2.weight[Chlmda_Bios_Demo_2_final$URXUCL == 1])/sum(Y_2.weight)
print(paste("The prevalence in 2001-2002 data is", round(prevalTrue_2, 5)))

WeightStand_2 <- Y_2.weight/sum(Y_2.weight)
#WeightStand_2

Y_2.ID       <- Chlmda_Bios_Demo_2_final$SEQN
Y_2.disease  <- Chlmda_Bios_Demo_2_final$URXUCL
X_2.marker_1 <- Chlmda_Bios_Demo_2_final$LBDMONO
X_2.marker_2 <- Chlmda_Bios_Demo_2_final$LBDNENO
X_2.marker_3 <- Chlmda_Bios_Demo_2_final$URXUMA

Y_2.data <- cbind(Y_2.ID, Y_2.disease, X_2.marker_1, X_2.marker_2, X_2.marker_3, Y_2.weight)
colnames(Y_2.data) <- c("SEQN", "Disease_status", "Biomarker_MONO", "Biomarker_EOSI", "Biomarker_Albumin", "Weight")
numSam_2 <- nrow(Y_2.data)
resamLOC_2 <- sample(1:numSam_2, size = numSam_2, replace = T, prob = WeightStand_2)
Data_resample_2 <- Y_2.data[resamLOC_2,]
#Data_resample_2


## 2003-2004
Y_3.weight <- Chlmda_Bios_Demo_3_final$WTINT2YR

### Calculate the true prevalence
prevalTrue_3 <- sum(Y_3.weight[Chlmda_Bios_Demo_3_final$URXUCL == 1])/sum(Y_3.weight)
print(paste("The prevalence in 2003-2004 data is", round(prevalTrue_3, 5)))

WeightStand_3 <- Y_3.weight/sum(Y_3.weight)
#WeightStand_3

Y_3.ID       <- Chlmda_Bios_Demo_3_final$SEQN
Y_3.disease  <- Chlmda_Bios_Demo_3_final$URXUCL
X_3.marker_1 <- Chlmda_Bios_Demo_3_final$LBDMONO
X_3.marker_2 <- Chlmda_Bios_Demo_3_final$LBDNENO
X_3.marker_3 <- Chlmda_Bios_Demo_3_final$URXUMA

Y_3.data <- cbind(Y_3.ID, Y_3.disease, X_3.marker_1, X_3.marker_2, X_3.marker_3, Y_3.weight)
colnames(Y_3.data) <- c("SEQN", "Disease_status", "Biomarker_MONO", "Biomarker_EOSI", "Biomarker_Albumin", "Weight")
numSam_3 <- nrow(Y_3.data)
resamLOC_3 <- sample(1:numSam_3, size = numSam_3, replace = T, prob = WeightStand_3)
Data_resample_3 <- Y_3.data[resamLOC_3,]
#Data_resample_3


## 2005-2006
Y_4.weight <- Chlmda_Bios_Demo_4_final$WTINT2YR

### Calculate the true prevalence
prevalTrue_4 <- sum(Y_4.weight[Chlmda_Bios_Demo_4_final$URXUCL == 1])/sum(Y_4.weight)
print(paste("The prevalence in 2005-2006 data is", round(prevalTrue_4, 5)))

WeightStand_4 <- Y_4.weight/sum(Y_4.weight)
#WeightStand_4

Y_4.ID       <- Chlmda_Bios_Demo_4_final$SEQN
Y_4.disease  <- Chlmda_Bios_Demo_4_final$URXUCL
X_4.marker_1 <- Chlmda_Bios_Demo_4_final$LBDMONO
X_4.marker_2 <- Chlmda_Bios_Demo_4_final$LBDNENO
X_4.marker_3 <- Chlmda_Bios_Demo_4_final$URXUMA

Y_4.data <- cbind(Y_4.ID, Y_4.disease, X_4.marker_1, X_4.marker_2, X_4.marker_3, Y_4.weight)
colnames(Y_4.data) <- c("SEQN", "Disease_status", "Biomarker_MONO", "Biomarker_EOSI", "Biomarker_Albumin", "Weight")
numSam_4 <- nrow(Y_4.data)
resamLOC_4 <- sample(1:numSam_4, size = numSam_4, replace = T, prob = WeightStand_4)
Data_resample_4 <- Y_4.data[resamLOC_4,]
#Data_resample_4


## 2007-2008
Y_5.weight <- Chlmda_Bios_Demo_5_final$WTINT2YR

### Calculate the true prevalence
prevalTrue_5 <- sum(Y_5.weight[Chlmda_Bios_Demo_5_final$URXUCL == 1])/sum(Y_5.weight)
print(paste("The prevalence in 2007-2008 data is", round(prevalTrue_5, 5)))

WeightStand_5 <- Y_5.weight/sum(Y_5.weight)
#WeightStand_5

Y_5.ID       <- Chlmda_Bios_Demo_5_final$SEQN
Y_5.disease  <- Chlmda_Bios_Demo_5_final$URXUCL
X_5.marker_1 <- Chlmda_Bios_Demo_5_final$LBDMONO
X_5.marker_2 <- Chlmda_Bios_Demo_5_final$LBDNENO
X_5.marker_3 <- Chlmda_Bios_Demo_5_final$URXUMA


Y_5.data <- cbind(Y_5.ID, Y_5.disease, X_5.marker_1, X_5.marker_2, X_5.marker_3, Y_5.weight)
colnames(Y_5.data) <- c("SEQN", "Disease_status", "Biomarker_MONO", "Biomarker_EOSI", "Biomarker_Albumin", "Weight")
numSam_5 <- nrow(Y_5.data)
resamLOC_5 <- sample(1:numSam_5, size = numSam_5, replace = T, prob = WeightStand_5)
Data_resample_5 <- Y_5.data[resamLOC_5,]
#Data_resample_5


## 2009-2010
Y_6.weight <- Chlmda_Bios_Demo_6_final$WTINT2YR

### Calculate the true prevalence
prevalTrue_6 <- sum(Y_6.weight[Chlmda_Bios_Demo_6_final$URXUCL == 1])/sum(Y_6.weight)
print(paste("The prevalence in 2009-2010 data is", round(prevalTrue_6, 5)))

WeightStand_6 <- Y_6.weight/sum(Y_6.weight)
#WeightStand_6

Y_6.ID       <- Chlmda_Bios_Demo_6_final$SEQN
Y_6.disease  <- Chlmda_Bios_Demo_6_final$URXUCL
X_6.marker_1 <- Chlmda_Bios_Demo_6_final$LBDMONO
X_6.marker_2 <- Chlmda_Bios_Demo_6_final$LBDNENO
X_6.marker_3 <- Chlmda_Bios_Demo_6_final$URXUMA

Y_6.data <- cbind(Y_6.ID, Y_6.disease, X_6.marker_1, X_6.marker_2, X_6.marker_3, Y_6.weight)
colnames(Y_6.data) <- c("SEQN", "Disease_status", "Biomarker_MONO", "Biomarker_EOSI", "Biomarker_Albumin", "Weight")
numSam_6 <- nrow(Y_6.data)
resamLOC_6 <- sample(1:numSam_6, size = numSam_6, replace = T, prob = WeightStand_6)
Data_resample_6 <- Y_6.data[resamLOC_6,]
#Data_resample_6








### Final Combination
Data_resample_final <- na.omit(rbind(Data_resample_1, Data_resample_2, Data_resample_3, Data_resample_4, Data_resample_5, Data_resample_6))
mean(Data_resample_final[ ,2] == 1) ## [1] 0.01805229
dim(Data_resample_final)[1]
write.csv(Data_resample_final, file = "./Data/Data_WholeBlood_final_test.csv")



