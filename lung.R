
#
library(Rtsne)
library(CancerSubtypes)
library("plot3D")
library("robustHD")
library("hyperSpec")
par(mar=c(1,1,1,1))
graphics.off()
set.seed(235)
#
# #Obtain the lung gene expression data as matrix. The data is in the form of samples as cols and genes(featues) as rows
lung_data <- as.matrix(read.delim("Lung/LUNG_Gene_Expression.txt"))
survivalData = read.table("Lung/LUNG_Survival.txt",sep="\t",header=TRUE)
lung_data[1:3,1:3]
dim(lung_data)
dim(survivalData)
d1 <- (dim(lung_data)[1])
d2 <- (dim(lung_data)[2])
#remove the first column because we dont need it
dummy1 <- lung_data[1:d1,2:d2]
index=which(is.na(dummy1))
dummy2 <- data.imputation(dummy1,fun="median")
dummy = data.normalization(dummy2,type="feature_Mean",log2=FALSE)
dim(dummy)
#transpose the data
lung_mat <- t(as.matrix(dummy))
dim(lung_mat)



lung_mat[1:3,1:3]

##########################################



####2. Normalise ranking and MAD
##### a.Normalise ranking
#GeneExp <- BRCA_gene #FxC
GeneExp<-t(lung_mat)
dim(GeneExp)#FxC
# dim(GeneExp)
# length(ranking)
# Gene_pageRank=ranking[1:nrow(GeneExp)]
# Gene_pageRank=Gene_pageRank/sum(Gene_pageRank)
# miRNA_pageRank=ranking[(nrow(GeneExp)+1):length(ranking)]
# miRNA_pageRank=miRNA_pageRank/sum(miRNA_pageRank)
##### b.Normalise MAD
Gene_MAD=as.vector(apply(GeneExp,1,mad))
Gene_MAD=Gene_MAD/sum(Gene_MAD)
# miRNAExp <- BRCA_miRNA
# miRNA_MAD=as.vector(apply(miRNAExp,1,mad))
# miRNA_MAD=miRNA_MAD/sum(miRNA_MAD)

##### c.combine ranking and MAD
# beta = 0.7
# Gene_combineweight = Gene_pageRank*beta+Gene_MAD*(1-beta)
Gene_combineweight = Gene_MAD
length(Gene_combineweight)
# miRNA_combineweight = miRNA_pageRank*beta+miRNA_MAD*(1-beta)
# length(miRNA_pageRank)
####3.Weighted similartiy fusion network
###weighted distance fuction
distanceWeighted<-function(X,weight)
  ##X is the expression Matrix(Row is featrue, column is sample)
{
  X_row = nrow(X)
  weight_diag<-diag(weight)
  X2<-(X^2)%*%weight_diag
  sumsqX = rowSums(X2)
  X1<-X%*%weight_diag
  XY = X1 %*% t(X)
  XX=matrix(rep(sumsqX, times = X_row), X_row, X_row)
  res=XX+t(XX)-2*XY
  res[res < 0] = 0
  diag(res)=0
  res<-sqrt(res)
  return(res)
}
dim(t(GeneExp))

tGeneExp<-t(GeneExp)
#Gw <- matrix(c(1,2,3,4,5,6),nrow=2,ncol=3)
#plot(log(Gene_combineweight),type='l')
length(Gene_combineweight)
Gw <- apply(tGeneExp,1, function(tGeneExp) tGeneExp*Gene_combineweight ) 
dim(Gw)
tGeneExp[1:3,1:3]
Gene_combineweight[1:3]
Gw[1:3,1:3]
tGw<-t(Gw)
dim(tGw) #RxF





par(mar=c(1,1,1,1))

set.seed(235)


  lung_tsne <- Rtsne(tGw, dim=3, 
                     perplexity =9,
                     theta = 0,
                     check_duplicates = FALSE,
                     pca = FALSE,
                     partial_pca = FALSE,
                     max_iter = 1000,
                     verbose = getOption("verbose", FALSE), 
                     is_distance = FALSE,
                     Y_init = NULL,
                     pca_center = FALSE, 
                     pca_scale = FALSE,
                     normalize = FALSE)
  #str(lung_tsne)
  fil <- as.data.frame(lung_tsne$Y)
  dim(fil)
  set.seed(235)
  
  library(otrimle)
 
    otri=otrimle(fil, 4, initial = NULL, logicd = NULL, npr.max =6/100, erc = 50,
                 iter.max = 100, tol = 1e-06, ncores = NULL, monitor = TRUE)
    
    
    otri$cluster
    ############%%%%%%%%%%%%%%%%%%%######################
    C<-otri$cluster
    
    
    
    ###############################
    set.seed(235)
    
    
    
    
    ########################################
    survivalData1 <- cbind(survivalData,C)
    colnames(survivalData1)<-c('PatientID', 'Survival', 'Death', 'C')
    vd<-length(which(C==0))
    if(vd<4 && vd > 0 ){
      survivalData1<-survivalData1[-which(survivalData1$C==0),]
      # i_j<-append(i_j,c(i,j))
      # Pvalue<-append(Pvalue,0)
      # Rled_value<-append(Rled_value,0)
      # next
    }
    if(length(which(C==1))<4){
      i_j<-append(i_j,c(i,j))
      Pvalue<-append(Pvalue,0)
      Rled_value<-append(Rled_value,0)
      next
    }
    if(length(which(C==2))<4){
      i_j<-append(i_j,c(i,j))
      Pvalue<-append(Pvalue,0)
      Rled_value<-append(Rled_value,0)
      next
    }
    
    
    if(length(which(C==3))<4){
      i_j<-append(i_j,c(i,j))
      Pvalue<-append(Pvalue,0)
      Rled_value<-append(Rled_value,0)
      next
    }
    if(length(which(C==4))<4){
      i_j<-append(i_j,c(i,j))
      Pvalue<-append(Pvalue,0)
      Rled_value<-append(Rled_value,0)
      next
    }
    
    ###########################################################################
    
    
    # #########################################
    
    # source("clust_sep.R")
    # xxxx<-clust_sep(c(survivalData1$Survival),c(survivalData1$C))
    # xxxx
    # Write CSV in R:write the survival dat alongwith a column of cluster indicies in t .csv file
    #write.csv(data_for_kme1, file = "1_Lung_Clustered_survival_TSNE_TDA.csv")
    write.csv(survivalData1, file = "1_Lung_Clustered_survival_TSNE_TDA.csv")
    
    # Loading
    library("survminer")
    # Fit survival curves
    require("survival")
    #Read the survivla data from the csv file
    survivalData1 <- read.csv("1_Lung_Clustered_survival_TSNE_TDA.csv")
    dim(survivalData)
    
    #LUNG2 <- lung2[1:111,2:4]
    #str(LUNG2)
    
    #str(lung2)
    #fit the survivla data into the function
    survivalData1<-as.data.frame(survivalData1)
    fitQ <- survfit(Surv(Survival,Death) ~ as.factor(C), data = survivalData1)
    
    #calculate Pvalue
    p <- surv_pvalue(fitQ)
    print("pvalue is")
    print(p[2])
    pval <- c(p[2])
    C
    # if(pval > 0.05) {
    #   print("pval not satisfied")
    #   array_null <- append(array_null,0)
    #   next
    # }
    # 
    cv<-otri$size
    # if(cv[1]==1) {
    #   print("Noise cluster Size is 1")
    #   i_j<-append(i_j,c(i,j))
    #   Rled_value<-append(Rled_value,0)
    #   Pvalue<-append(Pvalue,p[2])
    #   next
    # }
    # Drawing curves
    ggsurvplot(fitQ,
               pval = TRUE,
               palette = c("red","blue","green","black","orange"),
               conf.int = FALSE,
               font.x = c(28, "plain", "black"),
               font.family="Roman",
               font.y = c(28, "plain", "black"),
               pval.size=12,
               censor.size=4,
               font.title=c(30,"Roman", "black"),
               font.tickslab=c(23, "plain","black"),
               font.legend=c(30, "plain","black"),
               legend.title = "Fig.A",
               legend.labs = c("C=1", "C=2", "C=3", "C=4"),
               legend=c(0.8,0.75))
    #########################################################
    require(survival)
    require(survRM2)
    
    ## Computes L1 based measures of separation between survival curves
    ## note: higher ==> more separation
    ##
    ## Inputs:
    
    ##
    ##
    survival.curves.separation <- function(data, cluster, tau){
      
      ans <- list()
      tmp           <- table(cluster)
      ClusterLabel  <- as.numeric(names(tmp))
      ClusterSize   <- as.numeric(tmp)
      K             <- length(ClusterLabel)
      n             <- sum(ClusterSize)
      
      ## Set a time grid
      time_grid  <- min(data$Survival) : min(max(data$Survival), tau)
      time_scale <- 1/diff(range(time_grid))
      
      ## Create estimated survival curve matrix: [time x clusters]
      H <- matrix(NA, nrow=length(time_grid), ncol=K)
      colnames(H) <- paste('Cluster_', ClusterLabel, sep='')
      
      ## Estimate the KM curve clusterwise on a common support
      for(k in 1:K){
        ## Compute Kaplan-Meier estimator on the kth cluster
        km    <- survfit(Surv(Survival, Death)~1, data = data[cluster==ClusterLabel[k] , ] )
        ## Construct the KM estimator function
        KMfun <- stepfun(x=km$time[-1], y=km$surv)
        H[,k] <- KMfun(time_grid)
      }
      
      ## construct matrix of pairwise L1 distances
      D <- matrix(0, ncol=K, nrow=K)
      for (i in 1:K){
        for(j in 1:K)
          if(i!=j){
            D[i,j] <- D[j,i] <- sum( abs( H[ , i]  -  H[ , j] ))
          }
      }
      ## Some scaling is given so that these numbers are somewhow interpretable
      ## for the same number of clusters independently of the time interval
      D <- D * time_scale
      
      
      ## Metric 1: min pairwise L1 distance
      iut <- which(upper.tri(D, diag=FALSE))
      ans$L1min <- min(D[iut])
      
      ## Metric 2: compute the summed L1 distance of each
      ## cluster to the nearest one
      diag(D)   <- NA
      ans$L1sum <- sum( D[iut] ) / length(iut)
      
      
      return(ans)
    }
    
    
    ## Computes the discrepancy between survival curves in temrs of RMST
    ## (restricted mean survival time).
    
    ##
    ## Outputs: a list with multiple separation measures
    ##
    ##
    rmst.separation <- function(data, cluster, tau) {
      ans <- list()
      tmp           <- table(cluster)
      ClusterLabel  <- as.numeric(names(tmp))
      ClusterSize   <- as.numeric(tmp)
      K             <- length(ClusterLabel)
      n             <- sum(ClusterSize)
      
      
      ## Compute the minimum of the largest observed time in each of the two groups
      max.time <- rep(0, K)
      for (k in 1:K) {
        max.time[k] <- max(data$Survival[cluster == ClusterLabel[k]])
      }
      TAU  <- min(max.time, tau)
      
      
      ## Names
      ##    * RMST = restricted mean survival time
      ##    * LER  = life expectancy ratio
      ##
      ## LER: Life Expectancy Ratio  Matrix
      ##      LER[i,j] = LER[j,i] = max{RMST[i] / RMST[j], RMST[j] / RMST[i]}
      ##      note that here we don't have a baseline group so we define the ratio
      ##      always using in the denominator the group that have smaller RMST
      ##
      ## LED: Life Expectancy Difference
      ##    LED[i,j] = LED[j,i] = abs(RMST[i] - RMST[j])
      ##    note that here we don't have a baseline group so we define tha abs difference
      ##
      LER <- LED <-  matrix(0, ncol = K, nrow = K)
      for (i in 1:K) {
        for (j in 1:K)
          if (i != j) {
            ## First select data from  the two groups
            idx <- { cluster == ClusterLabel[i] | cluster == ClusterLabel[j]  }
            x   <- data[idx,]
            ##  Create a 0-1 vector, with gr==1 if cluster==ClusterLabel[i]
            gr0  <- cluster[idx]
            gr   <- ifelse(gr0 == ClusterLabel[i], 1, 0)
            u    <- rmst2(time = x$Survival, status = x$Death, arm = gr, tau = TAU)
            
            rmst_i <- u$RMST.arm1$rmst[1]
            rmst_j <- u$RMST.arm0$rmst[1]
            
            LER[i,j]  <- LER[j, i] <- max(rmst_i / rmst_j, rmst_j / rmst_i  )
            LED[i, j] <- LED[j, i] <- abs(rmst_i - rmst_j)
          }
      }
      
      ## index of the upper triangle
      iut <- which(upper.tri(LER, diag = FALSE))
      
      ## metric: min of pairwise LER discrepancy
      ans$LERmin <- min(LER[iut])
      
      ## metric: scaled summed pairwise LER discrepancy
      ans$LERsum <- sum(LER[iut]) / length(iut)
      
      ## metric: min of pairwise LED discrepancy
      ans$LEDmin <- min(LED[iut])
      
      ## metric: scaled summed pairwise LED discrepancy
      ans$LEDsum <- sum(LED[iut]) / length(iut)
      
      
      return(ans)
    }
    
    ###
    
    survPval = function(SurvDat,CLg,nYears=11){
      
      fit <- survfit(Surv(Survival, Death) ~ CLg,data = SurvDat,subset = Survival < (365 * nYears))
      suv <- survminer::ggsurvplot(fit, risk.table = TRUE, risk.table.height = 0.5,
                                   xlim = c(0,4000), break.time.by = 500, pval = TRUE)
      pVal <- survminer:::surv_pvalue(fit,method = "survdiff",data=SurvDat)
      return(list(fit=fit,suv = suv,pVal=pVal))
      
    }
    #########################################################
    ########################################################
    ########################################################
    library(survival)
    survRes = list()
    nYears = 11
    tau = 3724
    lung21<-survivalData1
    dim(survivalData1)
    str(lung21)
    #to remove all noise point
    # lung21<-lung21[order(lung21$C),]
    # l0<-length(which(lung21$C==0))
    # l0
    # dim(lung21)
    # lung212<-lung21[(l0+1):dim(lung21)[1],]
    # dim(lung212)
    #to remove the noise point if there is only one noise point.
    lung21<-lung21[order(lung21$C),]
    l0<-length(which(lung21$C==0))
    l0
    if(l0==1){
      lung21<-lung21[2:dim(lung21)[1],]
    }
    dim(lung21)
    
    l2d1 <- dim(lung21)[1]
    lung2<-lung21[1:l2d1,]
    dim(lung2)
    print('dim lung21 are: ')
    print(dim(lung21))
    
    ClustGroup = lung21[1:l2d1,5:5]
    SurvDat = lung21[1:l2d1,3:4]
    CLg = ClustGroup
    
    LNormDist = c(unlist(survival.curves.separation(SurvDat, CLg,tau)),
                  unlist(rmst.separation(SurvDat, CLg,tau)))
    
    str(LNormDist)
    print('the value is)))))))))))))))))))))))))))))))))))))))))))))))))')
    print(LNormDist[5])


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    
    
    
    # #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    dim(survivalData)
    sd<-as.data.frame(survivalData[,2:3])
    head(sd)
    sdc<-cbind(sd,C)
    dim(sdc)
    sdc_ordered<-sdc[order(sdc$C),]
    head(sdc_ordered)
    zer<-length(which(sdc$C==0))+1
    zer
    sdc_cut<-as.data.frame(sdc_ordered[zer:106,])
    str(sdc_cut)
    head(sdc_cut)
    dim(sdc_cut)
    fitcox <- coxph(Surv(as.numeric(Survival),as.numeric(Death)) ~ as.factor(C), data = (sdc_cut))
    
    summary(fitcox)
    ftest <- cox.zph(fitcox,terms=FALSE)
    ggcoxzph(ftest)
    summary(ftest)
    0.4405+0.2114+0.3869
    summary(fitcox)
    
    ###############################
    summary(ftest)
    summary(fitcox)
    hr_vec<-exp(fitcox$coefficients)
    len_hr<-length(hr_vec)
    for(m in seq(1,len_hr,by=1)){
      if(hr_vec[m]<1){
        hr_vec[m]<-1/hr_vec[m]
      }
    }
    
    
    ratio_vect<-c()
    for(h in seq(1,len_hr,by=1)){
      if(h==len_hr) {
        print("h equals to len_hr")
        next
      }
      
      ratio<-hr_vec[h]/hr_vec[h+1]
      if(ratio<1){
        ratio<-1/ratio
      }
      ratio_vect<-append(ratio_vect,ratio)
    }
    dsd<-c(hr_vec,ratio_vect)
    HR_min<-min(dsd)
    HR_min
    
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    source("clust_sep_2.R")
    xxxx<-clust_sep_2(c(sdc_cut$Survival),c(sdc_cut$C))
    xxxx
    




# 
# 
# 
# 
dim(survivalData)
sd<-as.data.frame(survivalData[,2:3])
head(sd)
sdc<-cbind(sd,C)
head(sdc)
sdc_ordered<-sdc[order(sdc[,3]),]
head(sdc_ordered)
zer<-length(which(sdc[,3]==0))+1
zer
sdc_cut<-as.data.frame(sdc_ordered[zer:587,])
str(sdc_cut)
head(sdc_cut)
fitcox <- coxph(Surv(as.numeric(Survival),as.numeric(Death)) ~ as.factor(C), data = sdc_cut)
#summary(fitcox)
ftest <- cox.zph(fitcox,terms=FALSE)
#SC_BRCA_rank
summary(ftest)
summary(fitcox)
plot(ftest)
ggcoxzph(ftest)#,font.main = 10,font.x = c(8,  "black"),font.y = c(8,  "black"), font.tickslab = c(7, "plain", "black"),caption = "A-Lung")
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#
# ###################################################
#
#pathway Analysis
library(otrimle)
library(survRM2)
library(survminer)
library(survival)
require(org.Hs.eg.db)
library(limma)
require(clusterProfiler)

# source("robust_clustering.R")
# source("robust_scaling.R")
# source("RSC_thresholding.R")
# source("survival_separation_5y.R")
# source("suvival_analysis_on_clusterings.R")
source("gene_mapping.R")
#source("mapping_symbol_to_entrez.R")
source("pathway_analysis.R")
C = c(1, 1, 3, 2, 1, 1, 2, 3, 1, 1, 2, 4, 1, 4, 1, 1, 4, 4, 4, 4, 4, 1, 2, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2,
      1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 2, 2, 2, 2, 3, 1, 2, 3, 1, 1, 2, 1, 1, 1, 1, 1, 1,
      1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1 ,1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1,
      1, 2, 2, 2, 1, 2, 1)
length(C)
lung_mat[1:3,1:3]
red_dat<-lung_mat
dim(red_dat)
pathway_matrix <- cbind(C,red_dat)
dim(pathway_matrix)
pathway_matrix[1:8,1:4]

pathway_df<-as.data.frame(pathway_matrix)
str(pathway_df)
dim(pathway_df)
pathway_df[1:3,1:3]
pathway_df <- pathway_df[order(pathway_df$C),]
pathway_df[1:10,1:3]
length(pathway_df$C)
dim(pathway_df)
#pathway_df<-pathway_df[8:106,]
clusterR<-pathway_df$C
red_data<-as.matrix(pathway_df[,2:12043])
#red_data<-as.matrix(pathway_df)
dim(red_data)
length(clusterR)
dim(red_data)
red_data[1:3,1:3]
#red_data1<-red_data[5:122,2:17900]+0.000000000001
#clusterR1<-clusterR[5:122]
FR = pathway_analysis(cluster=clusterR,red_data)
write.csv(FR$FR, file = "33Lung_Proposed.csv")
#Dotplot of relevant pathways
str(FR$formula_res)
write.table(FR$formula_res,file="rank_brca_res.txt")
write.csv(FR$formula_res,file="rank_brca_res.csv")
formula_res<-as.data.frame(read.csv("rank_brca_res.csv"))
dp2 = dotplot(FR$formula_res, x=~group,showCategory=4,font.size=15) +
ggplot2::facet_grid(~othergroup) +
ggplot2::scale_colour_gradient(low = "black", high = "white") +
ggplot2::theme(axis.text.y=ggplot2::element_text(size=ggplot2::rel(1.8)))+
#ggplot2::ggtitle("sdasdas")+
ggplot2::labs(title='Fig.3A')+
ggplot2::theme(plot.title = element_text(color="black", size=28,hjust =-1.6,vjust=-3))
#
plot(dp2)









#pathway Analysis
library(otrimle) 
library(survRM2)
library(survminer)
library(survival)
require(org.Hs.eg.db)
library(limma)
require(clusterProfiler)


source("mapping_symbol_to_entrez.R")
source("pathway_analysis.R")


red_dat<-lung_mat
dim(red_dat)
pathway_matrix <- cbind(C,red_dat)
dim(pathway_matrix)
pathway_matrix[1:3,1:3]

pathway_df<-as.data.frame(pathway_matrix)
str(pathway_df)
pathway_df[1:3,1:3]
pathway_df <- pathway_df[order(pathway_df$C),]
pathway_df[400:410,1:3]
length(pathway_df$C)
pathway_df[1:3,1:14]
clusterR<-pathway_df$C
dgk<-dim(pathway_df)[2]
red_data<-as.matrix(pathway_df[,2:dgk])
dim(red_data)
red_data[1:3,1:3]
clusterR
FR = pathway_analysis(cluster=clusterR,red_data)
write.csv(FR$FR, file = "Lung_proposed.csv")
#Dotplot of relevant pathways
str(FR$formula_res)
# write.table(FR$formula_res,file="rank_brca_res.txt")
# write.csv(FR$formula_res,file="rank_brca_res.csv")
# formula_res<-as.data.frame(read.csv("rank_brca_res.csv"))
dp2 = dotplot(FR$formula_res, x=~group,showCategory=4,font.size=15) + 
  ggplot2::facet_grid(~othergroup) + 
  ggplot2::scale_colour_gradient(low = "black", high = "white") + 
  ggplot2::theme(axis.text.y=ggplot2::element_text(size=ggplot2::rel(1.8)))+
  #ggplot2::ggtitle("sdasdas")+
  ggplot2::labs(title='Fig.3A')+
  ggplot2::theme(plot.title = element_text(color="black", size=28,hjust =-1.6,vjust=-3))

plot(dp2)
