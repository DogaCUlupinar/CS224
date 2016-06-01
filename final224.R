setwd("/Users/dulupinar/Documents/UCLA/Classes/Spring16/CS229/project")
library(stringr)
library(MASS)
library(ggbiplot)
library(cluster)
library(data.table)
library(matrixStats)
library(nsprcomp)
library(factoextra)
library(clue)
library(scatterplot3d)
### THE FUNCTIONS ###

#E-step function- Performs the expecation step for the EM algorthim for mixture modeling
#it uses current iteration parameters to recompute the missing data which in our case is
#the parameter that indicates which population each individual is from
E.step = function(Xij,Fjk,PIk,Rik,K){
  eps = 1e-6
  num_individuals = dim(Rik)[1]
  num_snps = dim(Fjk)[1]
  for (k in 1:K){
    for ( i in 1:num_individuals){
      if (i %% 100 == 0) print(i)
      sum_indiv = 0
      for (j in 1:num_snps){
        #log probability
        sum_indiv = sum_indiv + (1 - Xij[i,j])*log(1-Fjk[j,k]) + Xij[i,j]*log(Fjk[j,k])
      }
      Rik[i,k] = sum_indiv + log(PIk[1,k])
    }
  }
  
  #using log sum exp trick
  likelihood = sum(apply(Rik, 1, max))
  
  #Change back to likelihood
  maxi = matrix(apply(Rik,1,max))
  Rik = sweep(Rik,1,maxi,'-')
  Rik = exp(Rik)
  Rik = Rik + eps
  Rik = sweep(Rik,1,matrix(apply(Rik,1,sum)),'/')
  
  return(list( Rik = Rik,likelihood = likelihood)) #normailize by sum of each row (k populations)
}

#Mixture function which runs a mixture model on a population given the matrix of individuals by SNPS
# and the expected number of populations
mixture = function(Xij,K){
  #intialize data
  num_snps = ncol(Xij)
  num_individuals = nrow(Xij)
  
  Fjk  = matrix(runif(num_snps*K),num_snps,K)  #MxK
  PIk  = matrix(1/K,1,K)                              #1xK
  
  Rik  = matrix(runif(num_individuals*K),num_individuals,K)
  
  
  
  old_likelihood = -Inf 
  likelihood = 0
  while (abs(abs(likelihood) - abs(old_likelihood)) > 1){
    old_likelihood = likelihood
    
    ## E-step
    My_vars = E.step(Xij,Fjk,PIk,Rik,K)
    likelihood = My_vars$likelihood
    
    ## M-Step
    # PIk
    PIk  = t(matrix(apply(My_vars$Rik,2,mean)))
    
    # Fjk
    Fjk =  sweep(t(Xij) %*% My_vars$Rik,2,t(matrix(apply(My_vars$Rik,2,sum))),'/')
    print(My_vars$likelihood)
    
  }
  return(list(Rik = My_vars$Rik,PIk = PIk,Fjk = Fjk))
}


#adds label data to graphs
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}

#given a k for the kfold it gives that many divisions in data and makes the nth the test data
#indexing starts at 1
kfold = function(data,k = 10,n = 1){
  data_size = length(data)
  test_size = data_size/k
  start = max((n-1)*test_size,1)
  end = min(start + test_size,data_size)
  return(data[start:end])
}

#caculated F1 score
F1score <- function(cluster,populations,verbosity = 0){
  F1_score_total = 0
  pop_table = table(cluster,populations)
  distinct_populations = dimnames(pop_table)[[2]]
  for (population in distinct_populations){
    
    # precision = TP/(fp + tp)
    precision = max(pop_table[,population])/sum(pop_table[which.max(pop_table[,population]),])
    # recall= tp/(tp+fn)
    recall = max(pop_table[,population])/sum(pop_table[,population])
    #f1-score
    F1_score = 2*(precision * recall)/(precision + recall)
    F1_score_total = F1_score_total + F1_score
    
    if (verbosity > 1) print(sprintf("Pop:%s and f1 score: %d", population,F1_score))
  }
  if (verbosity > 0) print(table(cluster,populations))
  return(F1_score_total/length(distinct_populations))
}

### INITALIZE ALL THE DATA ###
CHROMOSOME_NUMBER = "22"
if ("data_init" %in% ls() == FALSE){
  task = "read in data"
  start_readin = Sys.time()
  print("Reintializing Missing Data")
  #read in data
  haploid = fread(sprintf("./data/chr-%s.geno.csv",CHROMOSOME_NUMBER),header=FALSE)
  individuals = read.table(sprintf("./data/chr-%s.ind",CHROMOSOME_NUMBER))
  if("snps" %in% ls() == FALSE) snps = read.table(sprintf("./data/chr-%s.snp",CHROMOSOME_NUMBER))
  
  #process data
  num_individuals = dim(haploid)[2]
  num_snps = dim(snps)[1]
  snp_id = paste("snp",str_split_fixed(snps$V1, ":", 2)[,2],sep = "")
  diploid = haploid[1:num_snps,seq(1,num_individuals,2),with = FALSE] + haploid[1:num_snps,seq(2,num_individuals,2), with = FALSE]
  half_hap = haploid[1:num_snps,seq(1,num_individuals,2),with = FALSE]
  
  colnames(diploid) = as.character(individuals[seq(1,num_individuals,2),1])
  rownames(diploid) = snp_id
  
  colnames(half_hap) = as.character(individuals[seq(1,num_individuals,2),1])
  rownames(half_hap) = snp_id
  
  half_hap = t(half_hap)
  diploid = t(diploid)
  data_init = 1
  total_time = Sys.time() - start_readin
  print(task)
  print(total_time)
  rm(haploid)
}

INTERVAL = 15 #based on the LD that is present according to online research

#split up data
global_populations = str_split_fixed(individuals[seq(1, num_individuals, 2), 3], ":", 2)[,2]
local_populations  = individuals[seq(1, num_individuals, 2), 3]
rows_wanted = which(global_populations != "AFR")
pop_reduce = global_populations[rows_wanted ]

dip = data.matrix(diploid[,seq(1,ncol(diploid),INTERVAL) ])
hip = data.matrix(half_hap[,seq(1,ncol(half_hap),INTERVAL) ])


### COMPUTATION TASKS ###

## pca
task = "compute pca"
start = Sys.time()
pca = prcomp(dip)
total = Sys.time() - start
print(task)
print(total)


## kmeans for whole genome
task = "kmeans for whole genome"
start = Sys.time()
bb = kmeans(diploid,4)
total = Sys.time() - start
print(task)
print(total)

## plain plot
plot((pca$sdev^2/sum(pca$sdev^2))[1:10],type = "o",xlab = "Principal Components",ylab = "Varaince Explained", main = "Scree Plot")
plot(pca$x[,1],pca$x[,2],xlab = "PC1", ylab = "PC2", main = "PC1 vs PC2")
pairs(pca$x[,1:3])
plot(pca$x[,1],pca$x[,2],col= c("black","red","blue","green")[as.factor(global_populations)])

## kmeans for PCA data
task = "kmeans for projected data"
start = Sys.time()
bbb = kmeans(pca$x[,1:4],4,nstart = 5) 
total = Sys.time() - start
print(task)
print(total)

bbb = kmeans(pca$x[,1:5],4,nstart = 5) 
clusplot(pca$x[,2:3],bbb$cluster,color = TRUE, shade = TRUE,lines = FALSE,main = "" )

#write.csv(pca$x[,1:4],sprintf("./data/chr-%s.geno.pca.csv",CHROMOSOME_NUMBER),row.names = FALSE)

## attempt with hierachical clustering did not work
if (FALSE){
  #heirarchical clustering
  dissimilarity = dist(pca$x[,1:5])
  hcluster = hclust(dissimilarity,method = "complete")
  clusMember = cutree(hcluster,4)
  table(clusMember,global_populations)
  hcd = as.dendrogram(hcluster)
  clusDendro = dendrapply(hcd, colLab)
  plot(hcd)
}

## Determining the number of clusters
kmax = 15
wss = sapply(1:kmax,function(k){kmeans(pca$x[,1:5],k,nstart = 10)$tot.withinss})
plot(1:kmax, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
fviz_nbclust(pca$x[,2:3], kmeans, method = "wss")
fviz_nbclust(pca$x[,1:4], hcut, method = "silhouette",hc_method = "complete") #this one gives us 4
fviz_nbclust(pca$x[,2:3], kmeans, method = "silhouette") #this one gives us three


## Sparse pca
pca = nsprcomp(hip,nrestart = 10,em_maxiter = 100,ncomp = 3,k = 300)
reduced_snps = union(union(which(pca$rotation[,1] != 0),which(pca$rotation[,2] != 0)),which(pca$rotation[,3] != 0))

## Mixture modeling
output = mixture(hip[,reduced_snps],K=4)
clustering = apply(output$Rik,1,which.max)
table(clustering,global_populations)
 
#prediction time
population_names = dimnames(table(global_populations))$global_populations 
indicies = sapply(population_names,function(x){which(global_populations == x)})
datum = pca_total
F1score(clustering,global_populations)


## kfold cross validation
average_f1score = 0
num_folds = 10
for (i in c(1:num_folds)){
  test = Reduce(function(x,y){c(x,y)},sapply(indicies,kfold,k=num_folds,n=i))
  train = seq(1,length(global_populations),1)[-test]
  trained_clusters = kmeans(datum[train,],4,nstart = 10)
  predicted = cl_predict(trained_clusters,datum[test,])
  average_f1score = average_f1score + F1score(predicted,global_populations[test])
  print(F1score(predicted,global_populations[test]))
}
print(average_f1score/num_folds)

## all of the pcas into one
pca_22 = read.csv(sprintf("./data/chr-%s.geno.pca.csv",22))
pca_21 = read.csv(sprintf("./data/chr-%s.geno.pca.csv",21))
pca_20 = read.csv(sprintf("./data/chr-%s.geno.pca.csv",20))
pca_total = cbind(pca_22, pca_20)
pca_total = cbind(pca_total,pca_21)
good = c(0.8481792,0.8781718,0.9070475)
bad = c(0.786469,0.786469,0.786469)
dependent = c(20,21,22)
plot(dependent,good,type ="b",ylim = c(.73,.95),col = "blue",ylab = "F1 score",xlab = "Chromosome")
lines(dependent,bad,type = "b", col = "red")
legend(x="topright", legend = c("Using PCA","Whole Chromsome") , col=c("blue","red"), pch=1)

## time benchmark
Pca_seconds = c(281,174,169)
kmeans_seconds = c(714)
plot(dependent,Pca_seconds,type ="b",ylim = c(100,800),col = "blue",ylab = "Run time (seconds)",xlab = "Chromosome")
lines(20,kmeans_seconds,type = "b", col = "red")
legend(x = "topright", legend = c("PCA","Kmeans Whole Chromosome") , col=c("blue","red"), pch=1)
