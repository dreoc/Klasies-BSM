####################
### LDA ANALYSIS ###
#################### 16366 23070 25539 36957 40910 54916 58641 62453 76377 80598 81208 84060 84925 98601
library(MASS); options(scipen = 10)
#### SELECT P FROM DUMMY ANALYSIS
p <- 26
scores <- data.frame(dat[,1:p], log(gpa$Csize),gp=as.character(gp))
names(scores); dim(scores)
experi <- c(which(gp == ("BUTCHERY")),  which(gp == ("TRAMPLE")),  which(gp == ("HYENA")))

success<-numeric(100000)
seed<-character(10000)
for (i in 1:length(success)){
  set.seed(i)
  cat("iter no.",i,"\n")
  numb<-sample(1:length(experi),.7*length(experi))
  # Create training and testing datasets
  # summary(pc)
  train = scores[numb,]
  test = as.data.frame(scores[-numb,])
  dim(train);dim(test)
  lda.fit <- NULL
  
  # classification rate for model
  lda.fit = lda(train[,-dim(train)[2]], train$gp, prior=c(1,1,1)/3, CV=F)
  pred.gp <- NULL
  pred.gp = predict(lda.fit, test[,-dim(test)[2]])
  table(pred.gp$class, test$gp)
  #success[i]<-sum(pred.gp$class==test$gp)/length(test$gp)
  success[i]<-sum(test$gp[which(test$gp=="BUTCHERY")]==pred.gp$class[which(test$gp=="BUTCHERY")])/length(test$gp[which(test$gp=="BUTCHERY")])
  #if (success[i]>.87){stop()}
  #pred.gp
  
}
hist(success)
range(success)
quantile(success, c(c(.025,.975)))
mean(success)
which(success>.90)
