## Author: Evan Wiewiora, University of Michigan
## Based off sample code from Dr. Veronica Berrocal, University of Michigan

library(ggplot2)
# Reading in the patient information
patient.data <- read.csv("~/Biostatistics/GSRA2019/ClinicalTrial/Data/patient.csv",header=T)
dim(patient.data)

subj.id <- as.character(patient.data$SSID)
trt.label <- as.character(patient.data$trt_label)

keep.subj <- which(as.character(trt.label)!="")

subj.id.keep <- subj.id[keep.subj]
un.subj.id <- unique(subj.id.keep)
trt.label.keep <- trt.label[keep.subj]
n.id.keep <- length(subj.id.keep)

id.renal <- c("01-0072","11-1506","31-4518","04-0453")

# Indices for the Placebo and Abatacept subjects and for subjects with
# renal crisis
index.pla <- which(as.character(trt.label.keep)=="Placebo")
index.aba <- which(as.character(trt.label.keep)=="Abatacept")
index.renal <- c(which(as.character(un.subj.id)==id.renal[1]),which(as.character(un.subj.id)==id.renal[2]),which(as.character(un.subj.id)==id.renal[3]),which(as.character(un.subj.id)==id.renal[4]))


# Reading the genetic subgroups in
intrinsic.set <- read.csv(file="~/Biostatistics/GSRA2019/ClinicalTrial/Data/Genetic_signature_baseline.csv",sep=",",header=T)
bsl.intra.subset <- as.character(intrinsic.set$Gene_signature)


## Indices for the subjects in the various subgroups
index.pla.p <- which(as.character(trt.label.keep)=="Placebo" & as.character(bsl.intra.subset)=="P")
index.aba.p <- which(as.character(trt.label.keep)=="Abatacept" & as.character(bsl.intra.subset)=="P")
index.pla.n <- which(as.character(trt.label.keep)=="Placebo" & as.character(bsl.intra.subset)=="N")
index.aba.n <- which(as.character(trt.label.keep)=="Abatacept" & as.character(bsl.intra.subset)=="N")
index.pla.i <- which(as.character(trt.label.keep)=="Placebo" & as.character(bsl.intra.subset)=="I")
index.aba.i <- which(as.character(trt.label.keep)=="Abatacept" & as.character(bsl.intra.subset)=="I")


## Read in parameter data

th2.mat <- read.table(file="~/Biostatistics/GSRA2019/ClinicalTrial/Data/Th2_data_work_dataset.csv",sep=",",header=T)

th2.visit.vec <- th2.mat$Visit
th2.trt.label <- th2.mat$Trt
th2.subj.id <- th2.mat$Subj.ID

il4p_il17p <- as.numeric(as.character(th2.mat$Th2_V4))
il17p <- as.numeric(as.character(th2.mat$Th2_V3))

treg.mat <- read.table(file="~/Biostatistics/GSRA2019/ClinicalTrial/Data/Treg_data_work_dataset.csv",sep=",",header=T)

treg.visit.vec <- treg.mat$Visit
treg.trt.label <- treg.mat$Trt
treg.subj.id <- treg.mat$Subj.ID

treg.cd25p <- as.numeric(as.character(treg.mat$Treg_V4))

tfh.mat <- read.table(file="~/Biostatistics/GSRA2019/ClinicalTrial/Data/Tfh_data_work_dataset.csv",sep=",",header=T)

pd1p.cxcr5n <- as.numeric(as.character(tfh.mat$Tfh_V14))

B.mat <- read.table(file="~/Biostatistics/GSRA2019/ClinicalTrial/Data/Bcell_data_work_dataset_Feb2019.csv",sep=",",header=T)

cd19pB <- as.numeric(as.character(B.mat$CD19.B))

un.visit <- unique(th2.visit.vec)

# parameters of interest

biom1.val <- cd19pB
biom2.val <- pd1p.cxcr5n

# biom1.time is a matrix with 88 rows and 4 columns (same for biom2.time)
# the first column is the value of the cell counts at baseline; the
# second one contains the change in the cell counts from baseline to
# 1 month; the third one contains the change in the cell counts from
# baseline to 3 months and so forth.

biom1.time <- matrix(biom1.val,nrow=length(un.subj.id),ncol=length(un.visit))
biom2.time <- matrix(biom2.val,nrow=length(un.subj.id),ncol=length(un.visit))

change.from.bsl.biom1 <- matrix(NA,nrow=length(un.subj.id),ncol=(length(un.visit)-1))
change.from.bsl.biom2 <- matrix(NA,nrow=length(un.subj.id),ncol=(length(un.visit)-1))

for(k in 1:(length(un.visit)-1)){
  change.from.bsl.biom1[,k] <- biom1.time[,(k+1)]-biom1.time[,1]
  change.from.bsl.biom2[,k] <- biom2.time[,(k+1)]-biom2.time[,1]
}

#Correlation tests

indices = list(index.aba.n, index.aba.i, index.aba.p, index.pla.n, index.pla.i,index.pla.p)

correlation_tests = function(t) {
  gen_sub = data.frame(correlation = rep(0,7), p = rep(0,7))
  for (j in 1:length(indices)){
    #test correlations of parameters for normal, inflammatory, and proliferative subgroups at time t
    temp = cor.test(biom1.time[indices[[j]],t],biom2.time[indices[[j]],t], use = "pairwise.complete.obs")
    gen_sub[j,] = c(as.numeric(temp$estimate),temp$p.value)
  }  
    
  #test correlations of parameters separated by treatment at time t 
  temp2 = cor.test(biom1.time[index.aba,t],biom2.time[index.aba,t], use = "pairwise.complete.obs")
  gen_sub[5,] = c(as.numeric(temp2$estimate),temp2$p.value)
  temp3 = cor.test(biom1.time[index.pla,t],biom2.time[index.pla,t], use = "pairwise.complete.obs")
  gen_sub[6,] = c(as.numeric(temp3$estimate),temp3$p.value)
  
  #test correlations of parameters in whole group at time t
  temp4 = cor.test(biom1.time[,t],biom2.time[,t], use = "pairwise.complete.obs")
  gen_sub[7,] = c(as.numeric(temp4$estimate),temp4$p.value)
  
  return(gen_sub)
}

#Get correlation matrices at 4 measurement times
cor_test1 = correlation_tests(1)
cor_test2 = correlation_tests(2)
cor_test3 = correlation_tests(3)
cor_test4 = correlation_tests(4)



#Plot correlations
# This is to put 4 plots in one page
##Individuals with adverse effects remain black
par(mfrow=c(2,2))
plot(biom1.time[,1],biom2.time[,1],xlab="IL4+IL17+ of CD4+",ylab="IL17+ of CD4",pch=19,col="black",
     main="IL4+IL17+ vs IL17+ \n Baseline - All subjects")
points(biom1.time[index.pla,1],biom2.time[index.pla,1],pch=19,col="coral")
points(biom1.time[index.aba,1],biom2.time[index.aba,1],pch=19,col="darkturquoise")
legend("topleft",col=c("coral","darkturquoise"),c("Placebo","Abatacept"),pch=rep(19,2))
points(biom1.time[index.renal,1],biom2.time[index.renal,1],pch=19,col="black")
abline(coef=summary(lm(biom2.time[index.pla,1]~biom1.time[index.pla,1]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba,1]~biom1.time[index.aba,1]))$coeff,col="darkturquoise",lwd=2,lty=1)


plot(biom1.time[,1],biom2.time[,1],xlab="IL4+IL17+ of CD4+",ylab="IL17+ of CD4",type="n",pch=19,col="black",
     main="IL4+IL17+ vs IL17+ \n Baseline - Normal-like")
points(biom1.time[index.pla.n,1],biom2.time[index.pla.n,1],pch=19,col="coral")
points(biom1.time[index.aba.n,1],biom2.time[index.aba.n,1],pch=19,col="darkturquoise")
abline(coef=summary(lm(biom2.time[index.pla.n,1]~biom1.time[index.pla.n,1]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba.n,1]~biom1.time[index.aba.n,1]))$coeff,col="darkturquoise",lwd=2,lty=1)

plot(biom1.time[,1],biom2.time[,1],xlab="IL4+IL17+ of CD4+",ylab="IL17+ of CD4",pch=19,type="n",col="black",
     main="IL4+IL17+ vs IL17+ \n Baseline - Inflammatory")
points(biom1.time[index.pla.i,1],biom2.time[index.pla.i,1],pch=19,col="coral")
points(biom1.time[index.aba.i,1],biom2.time[index.aba.i,1],pch=19,col="darkturquoise")
points(biom1.time[index.renal,1],biom2.time[index.renal,1],pch=19,col="black")
abline(coef=summary(lm(biom2.time[index.pla.i,1]~biom1.time[index.pla.i,1]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba.i,1]~biom1.time[index.aba.i,1]))$coeff,col="darkturquoise",lwd=2,lty=1)

plot(biom1.time[,1],biom2.time[,1],xlab="IL4+IL17+ of CD4+",ylab="IL17+ of CD4",pch=19,type="n",col="black",
     main="IL4+IL17+ vs IL17+ \n Baseline - Proliferative")
points(biom1.time[index.pla.p,1],biom2.time[index.pla.p,1],pch=19,col="coral")
points(biom1.time[index.aba.p,1],biom2.time[index.aba.p,1],pch=19,col="darkturquoise")
abline(coef=summary(lm(biom2.time[index.pla.p,1]~biom1.time[index.pla.p,1]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba.p,1]~biom1.time[index.aba.p,1]))$coeff,col="darkturquoise",lwd=2,lty=1)






par(mfrow=c(2,2))
plot(biom1.time[,2],biom2.time[,2],xlab="CD19+ B Cell",ylab="PD1+ CXCR5- of CD4",pch=19,col="black",
     main="CD19+ B Cell vs PD1+ CXCR5- \n Change Baseline to 1 Month - All subjects")
points(biom1.time[index.pla,2],biom2.time[index.pla,2],pch=19,col="coral")
points(biom1.time[index.aba,2],biom2.time[index.aba,2],pch=19,col="darkturquoise")
legend("topleft",col=c("coral","darkturquoise"),c("Placebo","Abatacept"),pch=rep(19,2))
points(biom1.time[index.renal,2],biom2.time[index.renal,2],pch=19,col="black")
abline(coef=summary(lm(biom2.time[index.pla,2]~biom1.time[index.pla,2]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba,2]~biom1.time[index.aba,2]))$coeff,col="darkturquoise",lwd=2,lty=1)


plot(biom1.time[,2],biom2.time[,2],xlab="CD19+ B Cell",ylab="PD1+ CXCR5- of CD4",type="n",pch=19,col="black",
     main="CD19+ B Cell vs PD1+ CXCR5- \n Change Baseline to 1 Month - Normal-like")
points(biom1.time[index.pla.n,2],biom2.time[index.pla.n,2],pch=19,col="coral")
points(biom1.time[index.aba.n,2],biom2.time[index.aba.n,2],pch=19,col="darkturquoise")
abline(coef=summary(lm(biom2.time[index.pla.n,2]~biom1.time[index.pla.n,2]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba.n,2]~biom1.time[index.aba.n,2]))$coeff,col="darkturquoise",lwd=2,lty=1)

plot(biom1.time[,2],biom2.time[,2],xlab="CD19+ B Cell",ylab="PD1+ CXCR5- of CD4",pch=19,type="n",col="black",
     main="CD19+ B Cell vs PD1+ CXCR5- \n Change Baseline to 1 Month - Inflammatory")
points(biom1.time[index.pla.i,2],biom2.time[index.pla.i,2],pch=19,col="coral")
points(biom1.time[index.aba.i,2],biom2.time[index.aba.i,2],pch=19,col="darkturquoise")
points(biom1.time[index.renal,2],biom2.time[index.renal,2],pch=19,col="black")
abline(coef=summary(lm(biom2.time[index.pla.i,2]~biom1.time[index.pla.i,2]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba.i,2]~biom1.time[index.aba.i,2]))$coeff,col="darkturquoise",lwd=2,lty=1)

plot(biom1.time[,2],biom2.time[,2],xlab="CD19+ B Cell",ylab="PD1+ CXCR5- of CD4",pch=19,type="n",col="black",
     main="CD19+ B Cell vs PD1+ CXCR5- \n Change Baseline to 1 Month - Proliferative")
points(biom1.time[index.pla.p,2],biom2.time[index.pla.p,2],pch=19,col="coral")
points(biom1.time[index.aba.p,2],biom2.time[index.aba.p,2],pch=19,col="darkturquoise")
abline(coef=summary(lm(biom2.time[index.pla.p,2]~biom1.time[index.pla.p,2]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba.p,2]~biom1.time[index.aba.p,2]))$coeff,col="darkturquoise",lwd=2,lty=1)





par(mfrow=c(2,2))
plot(biom1.time[,3],biom2.time[,3],xlab="CD19+ B Cell",ylab="PD1+ CXCR5- of CD4",pch=19,col="black",
     main="CD19+ B Cell vs PD1+ CXCR5- \n Change Baseline to 3 Month - All subjects")
points(biom1.time[index.pla,3],biom2.time[index.pla,3],pch=19,col="coral")
points(biom1.time[index.aba,3],biom2.time[index.aba,3],pch=19,col="darkturquoise")
legend("topleft",col=c("coral","darkturquoise"),c("Placebo","Abatacept"),pch=rep(19,2))
points(biom1.time[index.renal,3],biom2.time[index.renal,3],pch=19,col="black")
abline(coef=summary(lm(biom2.time[index.pla,3]~biom1.time[index.pla,3]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba,3]~biom1.time[index.aba,3]))$coeff,col="darkturquoise",lwd=2,lty=1)


plot(biom1.time[,3],biom2.time[,3],xlab="CD19+ B Cell",ylab="PD1+ CXCR5- of CD4",type="n",pch=19,col="black",
     main="CD19+ B Cell vs PD1+ CXCR5- \n Change Baseline to 3 Month - Normal-like")
points(biom1.time[index.pla.n,3],biom2.time[index.pla.n,3],pch=19,col="coral")
points(biom1.time[index.aba.n,3],biom2.time[index.aba.n,3],pch=19,col="darkturquoise")
abline(coef=summary(lm(biom2.time[index.pla.n,3]~biom1.time[index.pla.n,3]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba.n,3]~biom1.time[index.aba.n,3]))$coeff,col="darkturquoise",lwd=2,lty=1)

plot(biom1.time[,3],biom2.time[,3],xlab="CD19+ B Cell",ylab="PD1+ CXCR5- of CD4",pch=19,type="n",col="black",
     main="CD19+ B Cell vs PD1+ CXCR5- \n Change Baseline to 3 Month - Inflammatory")
points(biom1.time[index.pla.i,3],biom2.time[index.pla.i,3],pch=19,col="coral")
points(biom1.time[index.aba.i,3],biom2.time[index.aba.i,3],pch=19,col="darkturquoise")
points(biom1.time[index.renal,3],biom2.time[index.renal,3],pch=19,col="black")
abline(coef=summary(lm(biom2.time[index.pla.i,3]~biom1.time[index.pla.i,3]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba.i,3]~biom1.time[index.aba.i,3]))$coeff,col="darkturquoise",lwd=2,lty=1)

plot(biom1.time[,3],biom2.time[,3],xlab="CD19+ B Cell",ylab="PD1+ CXCR5- of CD4",pch=19,type="n",col="black",
     main="CD19+ B Cell vs PD1+ CXCR5- \n Change Baseline to 3 Month - Proliferative")
points(biom1.time[index.pla.p,3],biom2.time[index.pla.p,3],pch=19,col="coral")
points(biom1.time[index.aba.p,3],biom2.time[index.aba.p,3],pch=19,col="darkturquoise")
abline(coef=summary(lm(biom2.time[index.pla.p,3]~biom1.time[index.pla.p,3]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba.p,3]~biom1.time[index.aba.p,3]))$coeff,col="darkturquoise",lwd=2,lty=1)




par(mfrow=c(2,2))
plot(biom1.time[,4],biom2.time[,4],xlab="CD19+ B Cell",ylab="PD1+ CXCR5- of CD4",pch=19,col="black",
     main="CD19+ B Cell vs PD1+ CXCR5- \n Change Baseline to 6 Month - All subjects")
points(biom1.time[index.pla,4],biom2.time[index.pla,4],pch=19,col="coral")
points(biom1.time[index.aba,4],biom2.time[index.aba,4],pch=19,col="darkturquoise")
legend("topleft",col=c("coral","darkturquoise"),c("Placebo","Abatacept"),pch=rep(19,2))
points(biom1.time[index.renal,4],biom2.time[index.renal,4],pch=19,col="black")
abline(coef=summary(lm(biom2.time[index.pla,4]~biom1.time[index.pla,4]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba,4]~biom1.time[index.aba,4]))$coeff,col="darkturquoise",lwd=2,lty=1)


plot(biom1.time[,4],biom2.time[,4],xlab="CD19+ B Cell",ylab="PD1+ CXCR5- of CD4",type="n",pch=19,col="black",
     main="CD19+ B Cell vs PD1+ CXCR5- \n Change Baseline to 6 Month - Normal-like")
points(biom1.time[index.pla.n,4],biom2.time[index.pla.n,4],pch=19,col="coral")
points(biom1.time[index.aba.n,4],biom2.time[index.aba.n,4],pch=19,col="darkturquoise")
abline(coef=summary(lm(biom2.time[index.pla.n,4]~biom1.time[index.pla.n,4]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba.n,4]~biom1.time[index.aba.n,4]))$coeff,col="darkturquoise",lwd=2,lty=1)

plot(biom1.time[,4],biom2.time[,4],xlab="CD19+ B Cell",ylab="PD1+ CXCR5- of CD4",pch=19,type="n",col="black",
     main="CD19+ B Cell vs PD1+ CXCR5- \n Change Baseline to 6 Month - Inflammatory")
points(biom1.time[index.pla.i,4],biom2.time[index.pla.i,4],pch=19,col="coral")
points(biom1.time[index.aba.i,4],biom2.time[index.aba.i,4],pch=19,col="darkturquoise")
points(biom1.time[index.renal,4],biom2.time[index.renal,4],pch=19,col="black")
abline(coef=summary(lm(biom2.time[index.pla.i,4]~biom1.time[index.pla.i,4]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba.i,4]~biom1.time[index.aba.i,4]))$coeff,col="darkturquoise",lwd=2,lty=1)

plot(biom1.time[,4],biom2.time[,4],xlab="CD19+ B Cell",ylab="PD1+ CXCR5- of CD4",pch=19,type="n",col="black",
     main="CD19+ B Cell vs PD1+ CXCR5- \n Change Baseline to 6 Month - Proliferative")
points(biom1.time[index.pla.p,4],biom2.time[index.pla.p,4],pch=19,col="coral")
points(biom1.time[index.aba.p,4],biom2.time[index.aba.p,4],pch=19,col="darkturquoise")
abline(coef=summary(lm(biom2.time[index.pla.p,4]~biom1.time[index.pla.p,4]))$coeff,col="coral",lwd=2,lty=1)
abline(coef=summary(lm(biom2.time[index.aba.p,4]~biom1.time[index.aba.p,4]))$coeff,col="darkturquoise",lwd=2,lty=1)

