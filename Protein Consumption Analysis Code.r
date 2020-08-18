# The accompanied dataset gives estimates of the average protein consumption (in grams per person per day) 
# from different food sources for the inhabitants of 25 European countries.

#installing all the required packages
install.packages("knitr")
library(knitr)
install.packages("rmarkdown")
library(rmarkdown)
install.packages("ggplot2")
library(ggplot2)
install.packages("factoextra")
library(factoextra)
install.packages("dplyr")
library(dplyr)
install.packages("GGally")
library(GGally)
install.packages("cluster", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(cluster)
library(psych)
install.packages("psych", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")


############# reading the Protein Consumption dataset and attaching it #############
data<-read.csv("Protein_Consumption.csv", fill = TRUE)
#fill = True is added so that if there are rows which have unequal lengths or there is some missing data then it will fill implicitly
attach(data)

##### Check for the dimensions of the dataset attached ###########
dim(data)
#Ans- There are 25 observations and 11 variables
head(data)
tail(data)
#######################################################################

# Q1) Use principal components analysis to investigate the relationships between the countries on the basis of these variables

############## Performing PCA ##############################
View(data)
cor(data[-1])
# Removing the first variable from the dataset as it is a categorical variable
#while the correlation requires quantitative(numerical values)

data_pca<- prcomp(data[,-1],scale=TRUE)
# scale=TRUE:- the variable means are set to 0, and variances are set to 1
data_pca #the components for all the variables are displayed here
summary(data_pca)
# PC1 is able to restore 41% of the total variance, PC2 has 17% of total variance restored, PC3 has 13%, 
#PC4 has 10%, PC5 has almost 7%, PC6 has 4%, PC7 has 3%, PC8 has almost 2% while PC9 has 1% of the total variance restored

#Conclusion: When I add the variances of 7 Principal Components 95% of the total variance has been restored
#i.e. 95% of the estimated protein consumption comes from these 7 components
#To have a clear idea of choosing the principal components lets' draw the scree plot and then take a final call

# sample scores stored in data_pca$x
# singular values (square roots of eigenvalues) stored in data_pca$sdev
# loadings (eigenvectors) are stored in data_pca$rotation
# variable means stored in data_pca$center
# variable standard deviations stored in data_pca$scale
# A table containing eigenvalues and %'s accounted, follows

# Eigenvalues are sdev^2
eigen_data<-data_pca$sdev^2
eigen_data 

#the eigen values have no names to it so we will now assign the names to it
names(eigen_data) <- paste("PC",1:10,sep="")
eigen_data

#Calculating the sum of all the eigen values
sumlambdas<-sum(eigen_data)
sumlambdas

## Proportion of variance explained by each component
propvar <- eigen_data/sumlambdas
propvar

#Calculating the cumulative sum of proportion of the percentage of total variance
cumvar_data <- cumsum(propvar)
cumvar_data

#Putting all these values in a matrix format using the row-wise distribution
matlambdas <- rbind(eigen_data,propvar,cumvar_data)
matlambdas

# Giving apt names to these variables
rownames(matlambdas)<- c("Eigenvalues","Prop. variance","Cum. prop. variance")
matlambdas

#very big values are displayed in each of these components so I am rounding these values till 4 decimal places
round(matlambdas,4)

summary(data_pca)
print(data_pca$rotation)

# Sample scores stored in data_pca$x
data_pca$x

fviz_eig(data_pca)
summary(data_pca)

## Plot a biplot to view components on n-dimensional plane
biplot(data_pca, scale = 0, main = 'Principal componets')

plot(propvar, xlab = 'Prinicipal component',ylab = 'Proporiton of variance explained',type = 'b', main = 'Prop. of Variance')
#The optimum number of components are ~ 6 i.e PC1 : PC6

# cumulative scree plot
plot(cumvar_data,xlab = 'Principal component',ylab = 'Cumulative proportion of variance explained',type = 'b', main = 'Cumulative Prop.of Variance')
#Approx: ~ 92% of the variance is explained by 6 components i.e PC1 to PC6
############## End of PCA ########################
# Conclusion: PCA is a dimension -reduction tool that can be used to reduce a large set
# of variables to a small set that still contains most of the information in the large set.
# From the above scree plot I am concluding that I would want to consider
# the first 6 principal components as they help is restoring 92% of the total variance
# i.e the first 6 components will contribute to maximum average protein consumption (in grams per person per day)
##################################################

#Q2) Carry out cluster analysis to study relation between countries on their diet
########## Clustering Analysis ###################

data1<- read.csv("Protein_Consumption.csv", row.names=1, fill = TRUE)
matstd.prot<-scale(data1)

# Creating a (Euclidean) distance matrix of the standardized data 
dist.prot <- dist(matstd.prot, method="euclidean")

# Invoking hclust command (cluster analysis by single linkage method)      
clusprot.nn <- hclust(dist.prot)

# Plotting vertical dendrogram      
# create extra margin room in the dendrogram, on the bottom (Protein consumption' labels)
par(mar=c(6, 4, 4, 2) + 0.1)
plot(as.dendrogram(clusprot.nn),ylab="Distance between Protein Consumption",ylim=c(0,2.5),main="Dendrogram of protein consumption for inhabitants in Europe")

#####################################################################

# take a random sample of size 15 from a dataset of data1
# sample without replacement
mysample <- data1[sample(1:nrow(data1),15,replace=FALSE),]

# Standardizing the data with scale()
matstd.loan<- scale(mysample)
# Creating a (Euclidean) distance matrix of the standardized data
dist.employ <- dist(matstd.loan, method="euclidean")
# Invoking hclust command (cluster analysis by single linkage method)
clusemploy.nn <- hclust(dist.employ, method = "single")

#Plotting

# Create extra margin room in the dendrogram, on the bottom (loan labels)
par(mar=c(8, 4, 4, 2) + 0.1)
# Object "clusemploy.nn" is converted into a object of class "dendrogram"
# in order to allow better flexibility in the (vertical) dendrogram plotting.
plot(as.dendrogram(clusemploy.nn),ylab="Distance",ylim=c(0,6),
     main="Dendrogram.")

################### Agne Function ############################

# We will use agnes function as it allows us to select option for data standardization, the distance measure and clustering algorithm in one single function

(agn.employ <- agnes(mysample, metric="euclidean", stand=TRUE, method = "single"))

#  Description of cluster merging
agn.employ$merge

#Dendogram
plot(as.dendrogram(agn.employ), xlab= "Distance",xlim=c(8,0),
     horiz = TRUE,main="Dendrogram")

#Interactive Plots
plot(agn.employ,ask=TRUE)

################################################################

######################### K-means #################################
# Standardizing the data with scale()
matstd.employ <- scale(data1)

# K-means, k=2, 3, 4, 5
# Centers (k's) are numbers thus, 10 random sets are chosen

(kmeans2.employ <- kmeans(matstd.employ,2,nstart = 10))

# Computing the percentage of variation accounted for. Two clusters
perc.var.2 <- round(100*(1 - kmeans2.employ$betweenss/kmeans2.employ$totss),1)
names(perc.var.2) <- "Perc. 2 clus"
perc.var.2
fviz_cluster(kmeans2.employ,data=matstd.employ)
# Conclusion: Only 2 clusters are formed but the % of variance is 65%

# Computing the percentage of variation accounted for. Three clusters
(kmeans3.employ <- kmeans(matstd.employ,3,nstart = 10))
perc.var.3 <- round(100*(1 - kmeans3.employ$betweenss/kmeans3.employ$totss),1)
names(perc.var.3) <- "Perc. 3 clus"
perc.var.3
fviz_cluster(kmeans3.employ,data=matstd.employ)
#Conclusion: Three clusters with 52%of the variance is restored and there are three separate groups which are visible

# Computing the percentage of variation accounted for. Four clusters
(kmeans4.employ <- kmeans(matstd.employ,4,nstart = 10))
perc.var.4 <- round(100*(1 - kmeans4.employ$betweenss/kmeans4.employ$totss),1)
names(perc.var.4) <- "Perc. 4 clus"
perc.var.4
fviz_cluster(kmeans4.employ,data=matstd.employ)
# Conclusion: 4 clusters and 44% of the variance is stored in these clusters and there is an overlap

# Computing the percentage of variation accounted for. Five clusters
(kmeans5.employ <- kmeans(matstd.employ,5,nstart = 10))
perc.var.5 <- round(100*(1 - kmeans5.employ$betweenss/kmeans5.employ$totss),1)
names(perc.var.5) <- "Perc. 5 clus"
perc.var.5
fviz_cluster(kmeans5.employ,data=matstd.employ)
#Conclusion: 5 clusters with 37.5% of the total variance is restored and the clusters are overlapping

############# END of Clustering ###########################
# Conclusion: Clustering is an exploratory data analysis technique which helps in identifying
# subgroups within the dataset. When I select 4 clusters there is an overlap between the clusters and
# the percentage of variance restored is also only 44% but when I cluster for 3 there are clearly 
# three separate groups formed and 52% of the total variance is restored.
# Hence, I will chose 3 clusters in our problem statement where those countries in Europe who have 
# similar consumption of protein are placed in the same groups.

################################################################

# Q3) Identify the important factors underlying the observed variables 
# and examine the relationships between the countries with respect to these factors

########### Factor Analysis #######################################
#calculating the correlation matrix for all the numeric data in our dataset

corrm.emp = cor(data1)
corrm.emp
plot(corrm.emp)
#this is the correlation plot

#calculating the PCA and plotting these variances
euroemp_pca <- prcomp(data1, scale=TRUE)
summary(euroemp_pca) 
plot(euroemp_pca)
#looks like Pc1, pc2, pc3,pc4,pc5 restores maximum of variance

# A table containing eigenvalues and %'s accounted, follows. 
# Eigenvalues are the sdev^2
(eigen_euroemp <- round(euroemp_pca$sdev^2,2))
names(eigen_euroemp) <- paste("PC",1:10,sep="")
eigen_euroemp

sumlambdas <- sum(eigen_euroemp)
sumlambdas

propvar <- round(eigen_euroemp/sumlambdas,2)
propvar

cumvar_euroemp <- cumsum(propvar)
cumvar_euroemp

matlambdas <- rbind(eigen_euroemp,propvar,cumvar_euroemp)
matlambdas

rownames(matlambdas) <- c("Eigenvalues","Prop. variance","Cum. prop. variance")
rownames(matlambdas)

eigvec.emp <- euroemp_pca$rotation
print(euroemp_pca)

#Taking the first five PCs to generate linear combinations for all the variables
pcafactors.emp <- eigvec.emp[,1:5]
pcafactors.emp

# Multiplying each column of the eigenvector's matrix by the square-root of the corresponding eigenvalue in order to get the factor loadings
unrot.fact.emp <- sweep(pcafactors.emp,MARGIN=2,euroemp_pca$sdev[1:5],`*`)
unrot.fact.emp

# Computing communalities
communalities.emp <- rowSums(unrot.fact.emp^2)
communalities.emp

# Performing the varimax rotation. The default in the varimax function is norm=TRUE thus, Kaiser normalization is carried out
rot.fact.emp <- varimax(unrot.fact.emp)
View(unrot.fact.emp)
rot.fact.emp

#The print method of varimax omits loadings less than abs(0.1). In order to display all the loadings, it is necessary to ask explicitly the contents of the object $loadings
fact.load.emp <- rot.fact.emp$loadings[1:5,1:5]
fact.load.emp

#Computing the rotated factor scores for the 25 European Countries. 
scale.emp <- scale(data1)
scale.emp
#as.matrix(scale.emp)%*%fact.load.emp%*%solve(t(fact.load.emp)%*%fact.load.emp)

library(psych)
install.packages("psych", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
fit.pc <- principal(data1, nfactors=5, rotate="varimax")
fit.pc
round(fit.pc$values, 3)
fit.pc$loadings
# Loadings with more digits
for (i in c(5,1)) { print(fit.pc$loadings[[1,i]])}
# Communalities
fit.pc$communality 
#Cereals is able to restore 96% of the total variance

# Rotated factor scores
fit.pc$scores
# Play with FA utilities

fa.parallel(data1) # See factor recommendation
fa.plot(fit.pc) # See Correlations within Factors
fa.diagram(fit.pc) # Visualize the relationship
vss(data1) # See Factor recommendations for a simple structure
########### END of FCA ########################################################################
# Conclusion: The goal of FCA is to identify groups items when considered together and explains as much of the observed co-variance as possible. 
# I can see that out of 5 factors taken into consideration we have reduced the factors to two which contains most of the information in the dataset.
# When I consider 5 factors then I will be able to
# restore 87% of the total variance where RC1 contributes for most of the variance followed by RC5 which restores the second highest total variance.
# White meat, pulses, nuts, oilseeds, red meat, eggs, cereals and milk are the most important protein consumed products (in grams per person per day) in Europe
################################################################################################

############ Multinomial Linear Regression ########################
