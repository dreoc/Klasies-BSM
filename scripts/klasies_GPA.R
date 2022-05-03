rm(list=ls())
library(rgl) # 3d visualization
library(vegan) # for the adonis function
library(geomorph) # geometrics for landmark and procrustes
library(rstudioapi)  

setwd(dirname(rstudioapi:::getActiveDocumentContext()$path))
getwd()

#-------TLK machine-------------#
# setwd("~/Box/Desiree_cut marks/Harris/Rcode/TLK Code/Experimental_Kmeans/Kmeans_experimentaldown")
# files <- list.files(path = "~/Box/Desiree_cut marks/Harris/Rcode/TLK Code/Experimental_Kmeans/Kmeans_experimentaldown")
# files <- files[grep(x = files, pattern = ".nts")]
# files
# #
# 
# obs <- tk_choose.files()

#-------EOC machine-------------#
#setwd("C:/Users/eotarola/Box/Box_LCA/Cutmarks_project/Desiree_cut marks/Harris/Rcode/TLK Code/Experimental_Kmeans/Kmeans_experimentaldown")
files <- list.files(path ="../data/experimental/" , pattern=".nts", full.names = F, recursive = T)
files
#
obs <- c("../data/paleo/downsampledJH_Xa1_1_1 - Finaldwnsmpcoords.nts", 
         "../data/paleo/downsampledJH_Xa1_1_2 - Finaldwnsmpcoords.nts",
         "../data/paleo/downsampledJH_Xa1_1_3 - Finaldwnsmpcoords.nts", 
         "../data/paleo/downsampledJH_Xa3_1_1 - Finaldwnsmpcoords.nts",
         "../data/paleo/downsampledJH_Xa3_1_2 - Finaldwnsmpcoords.nts", 
         "../data/paleo/downsampledJH_Xa6_1_1 - Finaldwnsmpcoords.nts",
         "../data/paleo/downsampledJH_Xa6_1_3 - Finaldwnsmpcoords.nts", 
         "../data/paleo/downsampledJH_Xa8_1_1 - Finaldwnsmpcoords.nts",
         "../data/paleo/downsampledJH_Xa8_1_2 - Finaldwnsmpcoords.nts", 
         "../data/paleo/downsampledJH_Xb3_1_1 - Finaldwnsmpcoords.nts",
         "../data/paleo/downsampledJH_Xb4_1 - Finaldwnsmpcoords.nts",   
         "../data/paleo/downsampledJH_Xb8_1 - Finaldwnsmpcoords.nts",  
         "../data/paleo/downsampledJH_Xb8_3 - Finaldwnsmpcoords.nts")


outliers <- c(#### Outliers from Meshlab downsample ######
              "Butchery_17.12.3 - ",
              "Butchery_17.32.1 -",
              "Crocodile_1.1.3 - ",
              "Crocodile_2.1.10 - ",
              "Crocodile_2.1.101 - ",
              "Trampling_2.200.1 - ",
              "Trampling_4.401.1 - ",
              "Trampling_4.401.11 - ",
              "Trampling_4.405.6 - ",
              
              ### DIKIKA OUTLIERS HYENA -4/9//21 ####
              "Woo_B5-dorsal_4mat",
              "Butchery_17.12.1 - ",
              "Butchery_17.12.4 - ",
              "Butchery_17.17.3 - ",
              "Butchery_17.34.1 - ",
              "Butchery_18.20.2 - ",
              "Crocodile_1.1.27 - Final",
              "Crocodile_2.1.219 -",
              "Crocodile_2.1.24 - ",
              "Trampling_4.404.5 - ",
              "Trampling_4.406.8 - ",
              "Wolf_FIO.3252019.7 - ",
              "Hyena_44.2.12.5 - ",
              "Hyena_51.30.44.1 - ",
              "Hyena_51.58.51.4 - ",
              "Butchery_17.34.1 - ",
              "Butchery_18.20.3 - ",
              "Butchery_18.20.2 - ",
              "Crocodile_1.1.27 - ",
              "Crocodile_2.1.219 - ",
              "Trampling_2.200.4 - ",
              "Trampling_4.404.5 - ",
              "Wolf_TIM.2182019.2 - ",
              "Hyena_44.2.27.1 - ",
              "Hyena_51.30.44.1 - ",
              "Hyena_51.58.51.4 - ",
              
              
              #### Outliers from GPA ###########
              # "Butchery_17.12.4 - ",
              # "Butchery_17.34.1 - ",
              # "Butchery_17.5.1a",
              # "Butchery_18.20.2 - ",
              # "Butchery_18.32.1 - ",
              # "Crocodile_1.1.4 - ",
              # "Crocodile_1.2.11 - ",
              # "Crocodile_2.1.103",
              # "Crocodile_2.1.205 - ",
              # "Crocodile_2.1.8 - ",
              # "Crocodile_2.1.88 - ",
              # "Crocodile_4.1.73 - ",
              # "Trampling_2.200.4 - ",
              # "Trampling_4.402.2 - ",
              # "Trampling_4.404.5 - ",
              # "Wolf_FIO.3252019.7 - ",
              # "Wolf_TIM.2182019.2 - ",
              # "Wolf_TIM.2182019.8 - "
              
              ###KLASIES
              "Butchery_17.1.1.4 - ",
              "Butchery_17.7.9 - ",
              "Butchery_17.7.1 - ",
              "Butchery_17.5.1b - ",
              "Butchery_17.17.2 - ",
              "Butchery_17.31.4 - "
              # "Butchery_17.12.1 - ",
              # "Butchery_17.17.3 - ",
              # "Butchery_17.31.2 - ",
              # "Butchery_17.31.5 - ",
              # "Butchery_17.34.1 - ",
              # "Butchery_18.20.1 - ",
              # "Butchery_18.20.2 - ",
              # "Butchery_18.20.3 - ",
              # "Trampling_2.200.4 - ",
              # "Trampling_4.402.2 - ",
              # "Trampling_4.404.5 - ",
              # "Trampling_4.404.7 - ",
              # "Wolf_FIO.3252019.7 - ",
              # "Wolf_TIM.2182019.37 - ",
              # "Wolf_TIM.2182019.8 - "
)

out.locs <- tapply(X = outliers, INDEX = 1:length(outliers), FUN = grep, x=files)
files <- files [-out.locs]
files

# remove classes not in use
# files <- files[-grep(pattern = "Butchery", x = files)]
# files <- files[-grep(pattern = "Trampling", x = files)]
files <- files[-grep(pattern = "Crocodile", x = files)]
# files <- files[-grep(pattern = "Woo", x = files)]
# files <- files[-grep(pattern = "Antoine", x = files)]
files <- files[-grep(pattern = "EJS", x = files)]
#files <- files[-grep(pattern = "Wolf", x = files)]
files


# Create a new matrix of files
b <- readmulti.nts(paste("../data/experimental/",files, sep=""))
dim(b)
dimnames(b)[[3]]

#### RENAMING FILES
sp <- NULL
for (i in 1:length(dimnames(b)[[3]])){
  sp <- c(sp, unlist(strsplit(dimnames(b)[[3]][i], split=" - Final"))[1])
}
sp

#### RENAMING FILES
spl <- NULL
for (i in 1:length(dimnames(b)[[3]])){
  spl <- c(spl, unlist(strsplit(sp[i], split="../data/experimental/"))[2])
}
spl

new.names <- sub(x=spl, pattern="downsampledJH_", replacement = "")
new.names <- sub(x=spl, pattern="downsampled", replacement = "")

cbind(dimnames(b)[[3]], new.names)
dimnames(b)[[3]] <- new.names

# Create a new matrix of OBS files
bb <-  readmulti.nts(obs)
dim(bb)
dimnames(bb)[[3]]

##-------EOC machine-------------##
spl <- NULL
for (i in 1:length(dimnames(bb)[[3]])){
  spl <- c(spl, unlist(strsplit(dimnames(bb)[[3]][i], split="paleo/"))[2])
}
spl

# #-------TLK machine-------------#
# spl <- NULL
# for (i in 1:length(dimnames(bb)[[3]])){
#   spl <- c(spl, unlist(strsplit(dimnames(bb)[[3]][i], split="/"))[11])
# }
# spl

### 

dimnames(bb)[[3]] <- spl
spl <- NULL
for (i in 1:length(dimnames(bb)[[3]])){
  spl <- c(spl, unlist(strsplit(dimnames(bb)[[3]][i], split=" - "))[1])
}
spl
spl <- sub(x=spl, pattern="downsampledJH_", replacement = "")

dimnames(bb)[[3]] <- spl

new.names2 <- sub(x=spl, pattern="LCA_", replacement = "")

cbind(dimnames(bb)[[3]], new.names2)
dimnames(bb)[[3]] <- new.names2


#### Create a combined matrix of Experimental and Archaeological Marks from obs
dim(bb)
dim(b)
bc <- array(NA, c(dim(b)[1], dim(b)[2], dim(b)[3]+dim(bb)[3]))
bc[,,1:dim(b)[3]] <- b
bc[,,(dim(b)[3]+1):(dim(b)[3] + dim(bb)[3])] <- bb
dim(bc)


##################
### CREATE GPA ###
##################

surf <- matrix(1:1152, ncol=1)

# Set up sliding curves
# curves <- read.csv("~/Box/Desiree_cut marks/Harris/Rcode/TLK Code/SupportCode/CurveSlidingMatrix.csv")
# curve.sliders <- unique(sort(unlist(curves)))
# surf <- surf[-curve.sliders]

# gpa <- gpagen(A = bc)

# Create the GPA
system.time(
  gpa <- gpagen(A=bc
                ,surfaces=surf, ProcD=TRUE, curves= NULL)
)



#### READ IN ALREADY CREATED GPA
# gpa <- readRDS(file.choose())
# dimnames(gpa$coords)[[3]]<-c(new.names, new.names2)

### SAVE GPA
# today <- paste(strsplit(date(), split = " ")[[1]][c(2,3,4,5)], collapse = "_")
# today <- paste(strsplit(today, split=":")[[1]][1:3], collapse = "_")
# filename <- paste("../data/",today,".gpa.",gpa$slide.method,".",gpa$nsliders,"Curves.",gpa$nsurf,"surfs",".rds", sep="")
# saveRDS(gpa, file = filename)



# Plot the GPA in 3D
clear3d()
i=1
plot3d(gpa$coords[,,1],aspect=FALSE, col="red", size=5)
for (i in 1:dim(bc)[3]){
  points3d(gpa$coords[,,i],aspect=FALSE)
  
}

# Use the below code to quickly run through each individual data point in a GPA
# plot3d(gpa$coords[,,1],aspect=FALSE, col="red", size=5)
# points3d(gpa$coords[,,i],aspect=FALSE)
# i 
# files[i]
# i <- i+1


gpa$consensus
clear3d()
plot3d(gpa$consensus, aspect=FALSE, col="red", size=5)


co <- two.d.array(gpa$coords)
row.names(co) <- c(new.names, new.names2)

dim(co)
pc <- prcomp(co)
summary(pc)

dat <- pc$x


rownames(dat)
col <- rep("blue",length(rownames(dat)))
col
col[grep(x = rownames(dat), pattern="Butchery")] <- "dark green"
col[grep(x = rownames(dat), pattern="Antoine")] <- "dark green"
col[grep(x = rownames(dat), pattern="Woo")] <- "dark green"
#col[grep(x = rownames(dat), pattern="Crocodile")] <- "green"
col[grep(x = rownames(dat), pattern="Trampling")] <- "orange"
col[grep(x = rownames(dat), pattern="Wolf")] <- "salmon"
col[grep(x = rownames(dat), pattern="Hyena")] <- "red"
col

plot(dat,  type="n")
text(dat, labels=rownames(dat),col=col)


clear3d()
plot3d(dat,  type="n", aspect = FALSE)
text3d(dat, texts = rownames(dat),col=col)



gp <- rep("PALEO",length(col))
gp
gp[grep(x = rownames(dat), pattern="Butchery")] <- "BUTCHERY"
gp[grep(x = rownames(dat), pattern="Antoine")] <- "BUTCHERY"
gp[grep(x = rownames(dat), pattern="Woo")] <- "BUTCHERY"
#gp[grep(x = rownames(dat), pattern="Crocodile")] <- "CROC"
gp[grep(x = rownames(dat), pattern="Trampling")] <- "TRAMPLE"
gp[grep(x = rownames(dat), pattern="Wolf")] <- "WOLF"
gp[grep(x = rownames(dat), pattern="Hyena")] <- "HYENA"
gp


plot3d(dat, aspect = FALSE, type="n")
points3d(pc$x, col = col)
text3d(pc$x, texts = gp,col=col)

plot(dat, type = "n")
points(x=dat, col = col)
text(x= dat, labels = gp, col=col)


############################
### PERMUTATIONAL ANOVAS ###
############################

newcases <- which(gp=="BUTCHERY"|gp=="TRAMPLE"|gp == "HYENA"|gp == "WOLF")

mod1 <- adonis(dat ~ gp, permutations = 9999, method = "euclidean")
mod1

mod2 <- adonis(dat[newcases,] ~ gp[newcases] + log(gpa$Csize[newcases]) +gp[newcases]:log(gpa$Csize[newcases]), 
               permutations = 9999, method = "euclidean")
mod2

##### PAIRWISE ANOVAS
newcases <- which(gp=="BUTCHERY"|gp=="HYENA")
mod3 <- adonis(dat[newcases,] ~ gp[newcases] + log(gpa$Csize[newcases]) +gp[newcases]:log(gpa$Csize[newcases]), 
               permutations = 9999, method = "euclidean")
mod3
#####
newcases <- which(gp=="BUTCHERY"|gp=="TRAMPLE")
mod4 <- adonis(dat[newcases,] ~ gp[newcases] + log(gpa$Csize[newcases]), 
               permutations = 9999, method = "euclidean")
mod4
#####
newcases <- which(gp=="BUTCHERY"|gp=="WOLF")
mod5 <- adonis(dat[newcases,] ~ gp[newcases] + log(gpa$Csize[newcases]) +gp[newcases]:log(gpa$Csize[newcases]), 
               permutations = 9999, method = "euclidean")
mod5

#####
newcases <- which(gp=="TRAMPLE"|gp=="HYENA")
mod6 <- adonis(dat[newcases,] ~ gp[newcases] + log(gpa$Csize[newcases]) , 
               permutations = 9999, method = "euclidean")
mod6

#####
newcases <- which(gp=="WOLF"|gp=="HYENA")
mod7 <- adonis(dat[newcases,] ~ gp[newcases] , 
               permutations = 9999, method = "euclidean")
mod7
#####
newcases <- which(gp=="TRAMPLE"|gp=="WOLF")
mod8 <- adonis(dat[newcases,] ~ gp[newcases] + log(gpa$Csize[newcases]) +gp[newcases]:log(gpa$Csize[newcases]), 
               permutations = 9999, method = "euclidean")
mod8


####################
### LDA ANALYSIS ###
####################
library(MASS); options(scipen = 10)
#### SELECT P FROM DUMMY ANALYSIS
p <- 26
scores <- data.frame(dat[,1:p], log(gpa$Csize),gp=as.character(gp))
names(scores); dim(scores)
experi <- c(which(gp == ("WOLF")), which(gp == ("BUTCHERY")),  which(gp == ("TRAMPLE")),  which(gp == ("HYENA")))

# Create training and testing datasets
# summary(pc)
train = scores[c(experi),]
test = as.data.frame(scores[-c(experi),])
dim(train);dim(test)
lda.fit <- NULL
lda.fit = lda(train[,-dim(train)[2]], train$gp, prior=c(1,1,1,1)/4, CV=T)
lda.fit
summary(lda.fit)
#write.csv(lda.fit$posterior, file = "post.csv")

# classification rate for model
table(train$gp,lda.fit$class)
sum(lda.fit$class==train$gp)/length(train$gp)
sum(train$gp[which(train$gp=="BUTCHERY")]==lda.fit$class[which(train$gp=="BUTCHERY")])/length(train$gp[which(train$gp=="BUTCHERY")])
sum(train$gp[which(train$gp=="HYENA")]==lda.fit$class[which(train$gp=="HYENA")])/length(train$gp[which(train$gp=="HYENA")])
sum(train$gp[which(train$gp=="TRAMPLE")]==lda.fit$class[which(train$gp=="TRAMPLE")])/length(train$gp[which(train$gp=="TRAMPLE")])
sum(train$gp[which(train$gp=="WOLF")]==lda.fit$class[which(train$gp=="WOLF")])/length(train$gp[which(train$gp=="WOLF")])

#
lda.fit = lda(train[,-dim(train)[2]], train$gp, prior=c(1,1,1,1)/4, CV=F)
pred.gp <- NULL
pred.gp = predict(lda.fit, test[,-dim(test)[2]])
table(pred.gp$class, test$gp)
pred.gp


####################
### LDA ANALYSIS ###
####################
set.seed(62453)
numb<-sample(1:length(experi),.7*length(experi))
# Create training and testing datasets
# summary(pc)
train = scores[numb,]
test = as.data.frame(scores[-numb,])
dim(train);dim(test)
lda.fit <- NULL

# classification rate for model
lda.fit = lda(train[,-dim(train)[2]], train$gp, prior=c(1,1,1,1)/4, CV=F)
pred.gp <- NULL
pred.gp = predict(lda.fit, test[,-dim(test)[2]])
table(pred.gp$class, test$gp)
sum(pred.gp$class==test$gp)/length(test$gp)
sum(test$gp[which(test$gp=="BUTCHERY")]==pred.gp$class[which(test$gp=="BUTCHERY")])/length(test$gp[which(test$gp=="BUTCHERY")])

pred.gp

test = as.data.frame(scores[-c(experi),])
lda.fit = lda(train[,-dim(train)[2]], train$gp, prior=c(1,1,1,1)/4, CV=F)
pred.gp <- NULL
pred.gp = predict(lda.fit, test[,-dim(test)[2]])
table(pred.gp$class, test$gp)