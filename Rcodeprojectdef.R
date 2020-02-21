install.packages("AICcmodavg")
install.packages ("pastecs")
library("AICcmodavg")
library("pastecs")

critical.r <- function( n, alpha = .05 ) {
  df <- n - 2
  critical.t <- qt( alpha/2, df, lower.tail = F )
  critical.r <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
  return( critical.r )
}

##########################################################################
##                            ANALYSIS ONE                              ##
##########################################################################
##                                Part 1                                ##
##########################################################################

##We have to find a way to split RACE into two seperate columns (dummy variables) because it has three levels

lowbwt$RACE.w<-ifelse(lowbwt$RACE==1,1,0)
lowbwt$RACE.b<-ifelse(lowbwt$RACE==2,1,0)

##Keeping only the variables of interest

lowbwt.na <- subset(lowbwt,select=c("LOW","AGE","LWT","FTV","RACE.w","RACE.b"))
head(lowbwt.na)

##Basic and descriptive statics

stat.desc(lowbwt.na,basic=TRUE, desc=TRUE)

cor(lowbwt.na, method="spearman")
critical.r(189)
##the only significant interactions are AGE*LWT and AGE*FTV
##the only significant main effect is LWT, but since we can have interactions, we have to put AGE/FTV in there too

##Build all biologically feasible models

lowbwt.list=list()

##We chose to start with all four variables and add single interactions

lowbwt.list[[1]]=glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV, data=lowbwt.na, family=binomial(link="logit"))
lowbwt.list[[2]]=glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+AGE*FTV, data=lowbwt.na, family=binomial(link="logit"))
lowbwt.list[[3]]=glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+AGE*LWT, data=lowbwt.na, family=binomial(link="logit"))
lowbwt.list[[4]]=glm(LOW ~ LWT+AGE+FTV+AGE*FTV, data=lowbwt.na, family=binomial(link="logit"))
lowbwt.list[[5]]=glm(LOW ~ AGE+FTV+AGE*FTV, data=lowbwt.na, family=binomial(link="logit"))


##Evaluate Akaike score

lowbwt.modnames <- c("LWT+RACE.w+RACE.b+AGE+FTV","LWT+RACE.w+RACE.b+AGE+FTV+AGE*FTV","LWT+RACE.w+RACE.b+AGE+FTV+AGE*LWT",
                     "LWT+AGE+FTV+AGE*FTV","AGE+FTV+AGE*FTV")
lowbwt.aictab <- aictab(cand.set=lowbwt.list, modnames=lowbwt.modnames)
lowbwt.aictab

##########################################################################
##                                Part 2                                ##
##########################################################################

##take the full model and summarize it

full.model1 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+LWT*RACE.w+LWT*RACE.b+LWT*AGE+LWT*FTV+RACE.b*AGE+RACE.w*AGE+RACE.w*FTV+RACE.b*FTV+AGE*FTV, data=lowbwt.na, family=binomial(link="logit"))
summary(full.model1)

##We can see that RACE (RACE.w and RACE.b) *FTV has no significant effect on the model. RACE.w*FTV has the hightest, but because it's "all or nothing"
##we leave RACE.b as well

redefined.model1.1 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+LWT*RACE.w+LWT*RACE.b+LWT*AGE+LWT*FTV+RACE.b*AGE+RACE.w*AGE+AGE*FTV, data=lowbwt.na, family=binomial(link="logit"))
summary(redefined.model1.1)

##We chose to keep the main effects, and drop out the interactions with the highest p-value, this time it's RACE and AGE

redefined.model1.2 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+LWT*RACE.w+LWT*RACE.b+LWT*AGE+LWT*FTV+AGE*FTV, data=lowbwt.na, family=binomial(link="logit"))
summary(redefined.model1.2)

##LWT*AGE has the highest p-value

redefined.model1.3 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+LWT*RACE.w+LWT*RACE.b+LWT*FTV+AGE*FTV, data=lowbwt.na, family=binomial(link="logit"))
summary(redefined.model1.3)

##LWT*RACE has the highest p-value

redefined.model1.4 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+LWT*FTV+AGE*FTV, data=lowbwt.na, family=binomial(link="logit"))
summary(redefined.model1.4)

##LWT*FTV has the highest p-value

redefined.model1.5 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+AGE*FTV, data=lowbwt.na, family=binomial(link="logit"))
summary(redefined.model1.5)

##Now we are down to only one interaction, so we look at the main effects.
##Here we see that RACE has the highest score

redefined.model1.6 <- glm(LOW ~ LWT+AGE+FTV+AGE*FTV, data=lowbwt.na, family=binomial(link="logit"))
summary(redefined.model1.6)

##We get to the same model as in the previous part. AGE is not significant but we cannot 
##leave it out because it is involved in an interaction

##We make a graph of the AICs over the number of steps

AIC <- list(full.model1$aic,redefined.model1.1$aic,redefined.model1.2$aic,redefined.model1.3$aic,redefined.model1.4$aic,redefined.model1.5$aic,redefined.model1.6$aic)
xnames <- c("step 1","step 2","step 3", "step 4", "step 5","step 6","step 7")
x <- 1:7
plot(x,AIC,xaxt="n",type='o', main="AICs backwards 1")
axis(1,x,labels=xnames)

##########################################################################
##                            ANALYSIS TWO                              ##
##########################################################################
##                                Part 1                                ##
##########################################################################

lowbwt.na2 <- subset(lowbwt,select=c("LOW","AGE","LWT","FTV","RACE.w","RACE.b","SMOKE","PTL","HT","UI"))

###Basic and descriptive statics

stat.desc(lowbwt.na2,basic=TRUE, desc=TRUE)

cor(lowbwt.na2, method="spearman")
critical.r(189)
##the only significant interactions are LWT*AGE, LWT*HT, LWT*UI, AGE*FTV, SMOKE*PTL and PTL*UI
##the only significant main effects are LWT, SMOKE, PTL, HT and UI
##LWT*UI can be left out too


##Build all biologically feasible models

lowbwt.list2=list()

##We chose to start with all eight variables and single interactions

lowbwt.list2[[1]]=glm(LOW ~ LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI, data=lowbwt.na2, family=binomial(link="logit"))
lowbwt.list2[[2]]=glm(LOW ~ LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI+LWT*HT, data=lowbwt.na2, family=binomial(link="logit"))
lowbwt.list2[[3]]=glm(LOW ~ LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI+AGE*FTV, data=lowbwt.na2, family=binomial(link="logit"))
lowbwt.list2[[4]]=glm(LOW ~ LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI+SMOKE*PTL, data=lowbwt.na2, family=binomial(link="logit"))
lowbwt.list2[[5]]=glm(LOW ~ LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
lowbwt.list2[[6]]=glm(LOW ~ LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI+LWT*AGE, data=lowbwt.na2, family=binomial(link="logit"))
lowbwt.list2[[7]]=glm(LOW ~ LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI+AGE*FTV+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
lowbwt.list2[[8]]=glm(LOW ~ LWT+AGE+FTV+SMOKE+PTL+HT+UI+AGE*FTV+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
lowbwt.list2[[9]]=glm(LOW ~ LWT+AGE+FTV+SMOKE+PTL+UI+AGE*FTV+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
lowbwt.list2[[10]]=glm(LOW ~ LWT+AGE+FTV+PTL+HT+UI+AGE*FTV+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
lowbwt.list2[[11]]=glm(LOW ~ AGE+FTV+SMOKE+PTL+HT+UI+AGE*FTV+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))

##Evaluate Akaike score

lowbwt.modnames2 <- c("LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI","LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI+LWT*HT","LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI+AGE*FTV",
                      "LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI+SMOKE*PTL","LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI+PTL*UI","LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI+LWT*AGE",
                      "LWT+AGE+RACE.w+RACE.b+FTV+SMOKE+PTL+HT+UI+AGE*FTV+PTL*UI","LWT+AGE+FTV+SMOKE+PTL+HT+UI+AGE*FTV+PTL*UI","LWT+AGE+FTV+SMOKE+PTL+UI+AGE*FTV+PTL*UI",
                      "LWT+AGE+FTV+PTL+HT+UI+AGE*FTV+PTL*UI","AGE+FTV+SMOKE+PTL+HT+UI+AGE*FTV+PTL*UI")
lowbwt.aictab2 <- aictab(cand.set=lowbwt.list2, modnames=lowbwt.modnames2)
lowbwt.aictab2


##########################################################################
##                                Part 2                                ##
##########################################################################

##We start with the full model and summarize it

##full.model2 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*RACE.w+LWT*RACE.b+LWT*AGE+LWT*FTV+LWT*SMOKE+LWT*PTL+LWT*HT+LWT*UI+RACE.w*AGE+RACE.b*AGE+RACE.w*FTV+RACE.b*FTV+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*PTL+RACE.b*PTL+RACE.w*HT+RACE.b*HT+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*HT+AGE*UI+FTV*SMOKE+FTV*PTL+FTV*HT+FTV*UI+SMOKE*PTL+SMOKE*HT+SMOKE*UI+PTL*HT+PTL*UI+HT*UI, data=lowbwt.na2, family=binomial(link="logit"))
##summary(full.model2)

##HT*UI can be deleted, because it's an interaction between two binary variables, which results in mostly zeros

full.model2 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+LWT*PTL+LWT*HT+LWT*UI+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*PTL+RACE.b*PTL+RACE.w*HT+RACE.b*HT+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*HT+AGE*UI+FTV*SMOKE+FTV*PTL+FTV*HT+FTV*UI+SMOKE*PTL+SMOKE*HT+SMOKE*UI+PTL*HT+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(full.model2)

##PTL*HT has the highest p-value

redefined.model2.1 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+LWT*PTL+LWT*HT+LWT*UI+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*PTL+RACE.b*PTL+RACE.w*HT+RACE.b*HT+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*HT+AGE*UI+FTV*SMOKE+FTV*PTL+FTV*HT+FTV*UI+SMOKE*PTL+SMOKE*HT+SMOKE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.1)

##SMOKE*HT has the highest p-value

redefined.model2.2 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+LWT*PTL+LWT*HT+LWT*UI+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*PTL+RACE.b*PTL+RACE.w*HT+RACE.b*HT+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*HT+AGE*UI+FTV*SMOKE+FTV*PTL+FTV*HT+FTV*UI+SMOKE*PTL+SMOKE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.2)

##SMOKE*PTL has the highest p-value

redefined.model2.3 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+LWT*PTL+LWT*HT+LWT*UI+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*PTL+RACE.b*PTL+RACE.w*HT+RACE.b*HT+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*HT+AGE*UI+FTV*SMOKE+FTV*PTL+FTV*HT+FTV*UI+SMOKE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.3)

##FTV*HT has the highest p-value

redefined.model2.4 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+LWT*PTL+LWT*HT+LWT*UI+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*PTL+RACE.b*PTL+RACE.w*HT+RACE.b*HT+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*HT+AGE*UI+FTV*SMOKE+FTV*PTL+FTV*UI+SMOKE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.4)

##RACE*PTL has the highest p-value

redefined.model2.5 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+LWT*PTL+LWT*HT+LWT*UI+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*HT+RACE.b*HT+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*HT+AGE*UI+FTV*SMOKE+FTV*PTL+FTV*UI+SMOKE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.5)

##RACE*HT has the highest p-value

redefined.model2.6 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+LWT*PTL+LWT*HT+LWT*UI+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*HT+AGE*UI+FTV*SMOKE+FTV*PTL+FTV*UI+SMOKE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.6)

##LWT*HT has the highest p-value

redefined.model2.7 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+LWT*PTL+LWT*UI+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*HT+AGE*UI+FTV*SMOKE+FTV*PTL+FTV*UI+SMOKE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.7)

##HT has the highest p-value, but again we leave in the main effects, and only look at interactions first
##FTV*UI has the highest p-value

redefined.model2.8 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+LWT*PTL+LWT*UI+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*HT+AGE*UI+FTV*SMOKE+FTV*PTL+SMOKE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.8)

##LWT*PTL has the highest p-value

redefined.model2.9 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+LWT*UI+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*HT+AGE*UI+FTV*SMOKE+FTV*PTL+SMOKE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.9)

##LWT*UI has the highest p-value

redefined.model2.10 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*HT+AGE*UI+FTV*SMOKE+FTV*PTL+SMOKE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.10)

##AGE*HT has the highest p-value

redefined.model2.11 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*UI+FTV*SMOKE+FTV*PTL+SMOKE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.11)

##SMOKE*UI has the highest p-value
##Interesting to note is that, while the main effect HT has been rather high (+/-0.90), we could not remove it as
##it was still involved in an interaction. But after dropping the last interaction (AGE*HT), it all of a sudden
##became significant (0.01)

redefined.model2.12 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*UI+FTV*SMOKE+FTV*PTL+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.12)

##FTV*SMOKE has the highest p-value

redefined.model2.13 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*UI+FTV*PTL+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.13)

##FTV*PTL has the highest p-value

redefined.model2.14 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+RACE.w*SMOKE+RACE.b*SMOKE+RACE.w*UI+RACE.b*UI+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.14)

##RACE*UI has the highest p-value

redefined.model2.15 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+LWT*SMOKE+RACE.w*SMOKE+RACE.b*SMOKE+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.15)

##LWT*SMOKE has the highest p-value

redefined.model2.16 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+RACE.w*SMOKE+RACE.b*SMOKE+AGE*FTV+AGE*SMOKE+AGE*PTL+AGE*UI+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.16)

##AGE*UI has the highest p-value

redefined.model2.17 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+RACE.w*SMOKE+RACE.b*SMOKE+AGE*FTV+AGE*SMOKE+AGE*PTL+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.17)

##RACE*SMOKE has the highest p-value

redefined.model2.18 <- glm(LOW ~ LWT+RACE.w+RACE.b+AGE+FTV+SMOKE+PTL+HT+UI+AGE*FTV+AGE*SMOKE+AGE*PTL+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.18)

##Dropping this interaction with RACE, we are now allowed to also drop the main effect of RACE, as it has a high p-value for a while

redefined.model2.19 <- glm(LOW ~ LWT+AGE+FTV+SMOKE+PTL+HT+UI+AGE*FTV+AGE*SMOKE+AGE*PTL+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(redefined.model2.19)

##We make a graph of the AICs over the number of steps

AIC2 <- list(full.model2$aic, redefined.model2.1$aic, redefined.model2.2$aic, redefined.model2.3$aic, redefined.model2.4$aic,
        redefined.model2.5$aic, redefined.model2.6$aic, redefined.model2.7$aic, redefined.model2.8$aic, redefined.model2.9$aic,
        redefined.model2.10$aic, redefined.model2.11$aic, redefined.model2.12$aic, redefined.model2.13$aic, redefined.model2.14$aic,
        redefined.model2.15$aic, redefined.model2.16$aic, redefined.model2.17$aic, redefined.model2.18$aic, redefined.model2.19$aic)
xnames <- c("step 1","step 2","step 3", "step 4", "step 5","step 6","step 7","step 8","step 9","step 10", "step 11", "step 12","step 13",
            "step 14","step 15","step 16","step 17", "step 18", "step 19","step 20")
x<-1:20
plot(x,AIC2, xaxt="n",type='o',main="AICs backwards 2")
axis(1,x,labels=xnames)

##########################################################################
##                             Conclusion                               ##
##########################################################################

##Conclusions of the last model building (8 variables)

conclusion.model <- glm(LOW ~ LWT+AGE+FTV+SMOKE+PTL+HT+UI+AGE*FTV+PTL*UI, data=lowbwt.na2, family=binomial(link="logit"))
summary(conclusion.model)
exp(conclusion.model$coefficients)
