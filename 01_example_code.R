# Example code: 
# Below is the code for a 3-facet design using the indicator mean method with ULS estimator. 
# make sure to run the gtheory_lvaan function code first.


# G-theory analysis results (lavaan object)
Results <- gtheory_lavaan(data= mydata, facet= c(f1="I", f2="O", f3="S"),  # format of the variable names in this data: "O1S1I1", "O1S1I2"...
                          d_n= c(12, 1, 4), method="indicator_mean",   # For the D-study, n_I'= 12, n_O'=1, n_S'=4
                          estimator="ULS")


# Monte Carlo confidence intervals for G, D coefficients, and variance components
set.seed(123)
output <- monteCarloCI(Results, level=0.95) # 95% confidence intervals 
output <- output[which(rownames(output)=="Gcoef"):nrow(output),]  
output # print the results


# Cut score-specific D coefficient graph
attach(as.list(setNames(output[3:nrow(output), "est"],
         rownames(output)[3:nrow(output)])))
d_n = c(12, 1, 4) # number of conditions for each facet (D-study)
np = 511 # number of sample
       
total_vc <- VC_p + ((VC_pXf1 + VC_f1)/d_n[1]) + ((VC_pXf2 + VC_f2)/d_n[2]) + ((VC_pXf3 + VC_f3)/d_n[3]) +
  ((VC_pXf1Xf2 + VC_f1Xf2)/(d_n[1]*d_n[2])) + ((VC_pXf1Xf3 + VC_f1Xf3)/(d_n[1]*d_n[3])) + ((VC_pXf2Xf3 + VC_f2Xf3)/(d_n[2]*d_n[3])) +
  ((VC_pXf1Xf2Xf3 + VC_f1Xf2Xf3)/(d_n[1]*d_n[2]*d_n[3]))

ybar <- VC_p/np + ((VC_pXf1/np + VC_f1)/d_n[1]) + ((VC_pXf2/np + VC_f2)/d_n[2]) +
  ((VC_pXf3/np + VC_f3)/d_n[3]) + ((VC_pXf1Xf2/np + VC_f1Xf2)/(d_n[1]*d_n[2])) + 
  ((VC_pXf1Xf3/np + VC_f1Xf3)/(d_n[1]*d_n[3])) + ((VC_pXf2Xf3/np + VC_f2Xf3)/(d_n[2]*d_n[3])) +
  ((VC_pXf1Xf2Xf3/np + VC_f1Xf2Xf3)/(d_n[1]*d_n[2]*d_n[3]))

score.grid <- seq(1,8, length.out=100) # 8 category data
mean <- mean(rowMeans(data))

cutscore <- data.frame(score = score.grid)
for (i in 1:length(score.grid)) {
  cutscore$D_coef[i] <- (VC_p + (mean-(cutscore[i,1]))^2 - ybar)/
    ((VC_p + (mean - (cutscore[i,1]))^2 - ybar) + (total_vc - VC_p))
}

plot(cutscore, 
     type="l",
     xlab = "Cut Score", 
     ylab = "D coefficient",
     las=1)
