## Project 4 ancestry mapping

I did project 4, ancestry mapping for the hard difficulty of determining the global ancestry from unknown populations. Here you will find all the code for the project in final224.R and the presentation i gave in discussion Friday May 20th.

## Approach

I used the unsupervised clustering machine learning technique, kmeans to cluster the indivudals into k different populations where our k is 4 {Asian,European,American,African}. In order to do so i needed to reduce the dimensionality of the data. This was accomplished using Principal Component Analysis (PCA) and accounting for LD. I took snp out of every 15 (this was according to data i found online for the size of a typical LD block) and on this reduced dataset i ran a PCA. Looking at the first few PC's i was able to classify the data quite well. This process was repeated for chromosomes 22 ,21 and 20 using the 1000 genomes data set. The PCA's from the different chromosomes can be combine and the kmeans can be run over these new clusters. For more information and graphics please refer to the presentation

## Testing

A k-fold cross validation sheme was used to assess the classification done by my method. I used 10 folds where i used divided the data from each population into to 10 slices and trained on the other 9 and predicted the remaining 1 slice. Accuracy improved with each addition of each new chromosome PCA data set.

## Additional Exploratory work - Sparse PCA

I performed a Sparse PCA to determine the min # of SNPs required (and which snps were required) to get the same classification accuracy as working on the entire Chromosome. Then using these reduced set of SNPS we can try to pinpoint which genes are most expressive in determining the population structure in a group of individuals 

## Additional Exploratory work - Mixture Modeling

In addition i implemented a simple mixture model and ran it on the reduced SNPs data set. It is a simplified version of admixture modeling where I state that the population for an individual is selected from the multinomial distribution governed how the prior P_pop of the percentatge of study individuals from each population.A given SNP for a given individual was drawn using the Bernoiulli distrbution dependent on parameter F, the allele frequency of that gene at that population. I learned the parameters F and P_pop using expectation maximization algorithim. The results were unfortunately pretty bad, with an F1 score of .741 where when kmeans was applied on the same data we get an F1score of  7.45

## The code set up

The code is set up into different regions to perform different tasks. The entire script is not meant to be run all at once. Please refer to the comments in the code to run what you desire.Split into three major sections, function implementations, data initialization, and all machine learning tasks(most of the project here is calling the function i wrote and benchmarking)
