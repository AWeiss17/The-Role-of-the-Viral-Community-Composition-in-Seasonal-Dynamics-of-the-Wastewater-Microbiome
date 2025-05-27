
########## After Nils approach #################

#####################################
### Load Required Libraries
#####################################
library(lattice)
library(ggplot2)
library(plyr)
library(tidyr)
library(vegan)
library(dplyr)
library(ape)
library(readr)


#####################################
### Step 1: Load tables and set rownames
#####################################
data_pr2_wide_SRR <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/count_data_pr2_wide_SRR.csv")
data_silva_wide_SRR <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/count_data_silva_wide_SRR.csv")

# Load Environmental Data
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")


#####################################
### Step 2: Combine PR2 and SILVA Data (RNA)
#####################################

data_RNA <- merge(
  data_pr2_wide_SRR, data_silva_wide_SRR, 
  by = "...1",  # Merge by first column (sample ID)
  all = TRUE    # Keep all rows, fill missing values with NA
)

#####################################
### Step 3: Set Row Names for Count Tables & Environmental Data
#####################################
# Convert to Data Frame & Set Row Names for RNA Data
data_RNA <- as.data.frame(data_RNA)
rownames(data_RNA) <- data_RNA$...1
data_RNA <- data_RNA[, -1]  # Remove first column after setting row names

# Convert Environmental Data & Set Row Names
Environmental_Data_RNA <- as.data.frame(Environmental_Data_RNA)
rownames(Environmental_Data_RNA) <- Environmental_Data_RNA$SRR_ID
Environmental_Data_RNA <- Environmental_Data_RNA[, -12]


#####################################
### Step 4: Further Environmental Table Preparation
#####################################

# Convert categorical variables to factors
Environmental_Data_RNA$Month <- as.factor(Environmental_Data_RNA$Month)
Environmental_Data_RNA$Season <- as.factor(Environmental_Data_RNA$Season)
Environmental_Data_RNA$Sample_ID <- as.factor(Environmental_Data_RNA$Sample_ID)

envcols_RNA <- colnames(Environmental_Data_RNA)


#####################################
### Step 5: Normalize the Taxa-Data (relative Counts)
#####################################
# Normalize datasets row-wise
data_relative_RNA <- sweep(data_RNA, 1, rowSums(data_RNA), FUN = '/')


#####################################
### Step 6: Combine Count Table with Environmental Data
#####################################
data_relative_RNA_env <- cbind(data_relative_RNA, Environmental_Data_RNA[rownames(data_relative_RNA), ])



#############################################################################################
########################### From here I followed Nils script ################################
#############################################################################################



#####################################
### Step 1: Define some functions
#####################################

`%notin%` <- Negate(`%in%`)

#defining a function which returns significance stars of a given p-value
sign_stars <- function(p){
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  stars=NA
  if(is.numeric(p)){
    if(p>=0 && p<=0.001){
      stars <- "***"
    } else if(p>0.001 && p<=0.01){
      stars <- "**"
    } else if(p>0.01 && p<=0.05){
      stars <- "*"
    } else if(p>0.05 && p<=0.1){
      stars <- "."
    } else if(p>0.01 && p<=1){
      stars <- " "
    } else{
      print(simpleError("p must be 0<=p<=1"))
    }
    if(!is.na(stars)){
      return(stars)
    }
  } else{
    print(simpleError("p must be numeric"))
  }
}

# function to clamp a value between a min and max value, used in color gradient generation
clamp=function(x,minimum=-Inf,maximum=Inf){
  return(min(c(maximum,max(c(x,minimum)))))
}

# Evaluate NMDS goodness
draw_stress_plot=function(){
  par(mfrow = c(1,2))  # Create a two-panel plot
  stressplot(usedNMDS)  # Shepard's stress plot: evaluates NMDS goodness-of-fit
  plot(usedNMDS, display = 'sites', type = 't', main = 'Goodness of fit')  # Ordination diagram
  points(usedNMDS, display = 'sites', cex = goodness(usedNMDS) * 200)  # Add points with size reflecting goodness-of-fit
  par(mfrow = c(1,1))  # Reset to default plotting layout
}

#####################################
### Step 2: Do some Transformations
#####################################

# To handle skewness in taxa abundance data, do transformation.
# Raw relative abundance (from e.g. data_pr2) might have extreme differences between dominant and rare taxa, 
# potentially skewing clustering or ordination
# log Transformation compresses high values more aggressively, while the square-root transformation is less severe
# Both aim to reduce the dominance of high-abundance data.
# Transformations are only applied to taxa data (not Environmental data)
# Transformations are not necessary for NMDS but for other analyses
# NMDS directly uses relative abundances

# log.transform data
logdata_RNA <- data.frame(data_relative_RNA_env[,which(colnames(data_relative_RNA_env) %in% envcols_RNA)],lapply(data_relative_RNA_env[,-which(colnames(data_relative_RNA_env) %in% envcols_RNA)]+1,log))

#squareroot transform data
sqrdata_RNA <- data.frame(data_relative_RNA_env[,which(colnames(data_relative_RNA_env) %in% envcols_RNA)],lapply(data_relative_RNA_env[,-which(colnames(data_relative_RNA_env) %in% envcols_RNA)],sqrt))

#create species data tables
logdata_RNA_spec <- logdata_RNA[,-which(colnames(logdata_RNA) %in% envcols_RNA)]

sqrdata_RNA_spec <- sqrdata_RNA[,-which(colnames(sqrdata_RNA) %in% envcols_RNA)]


#####################################
### Step 3: Dissimilarity Trees (Cluster Dendrograms)
#####################################
# Cluster dendrograms are tree-like visualization that shows how samples cluster together based on their pairwise dissimilarities.
# They allow to investigate how similar or different samples are, and group them into clusters.
# By comparing untransformed, log-transformed, and square-root-transformed data, 
# one can evaluate the effect of tansformations on clustering.
# vegdist() calculates dissimilarities between samples using Bray-Curtis
# hclust() performs hierarchical clustering using "complete linkage" (= defines the distance between two clusters as the max. distance between their individual elements)

# the 3 variants (untransformed, log, sqr) are chosen to evaluate how transformations affect sample clustering
# Dissimilarity tree for untransformed data RNA
plotree <- hclust(vegdist(data_relative_RNA), "complete") 
plot (plotree, cex = .6)
rect.hclust(plotree, 8) # define 8 clusters

# Dissimilarity tree for log data
plotree_log <- hclust(vegdist(logdata_RNA_spec), "complete")
plot (plotree_log, cex = .6)
rect.hclust(plotree_log, 8) # define 8 clusters

# Dissimilarity tree for square-root data
plotree_sqr <- hclust(vegdist(sqrdata_RNA_spec), "complete")
plot (plotree_sqr, cex = .6)
rect.hclust(plotree_sqr, 8) # define 8 clusters



#####################################
### Step 4: Shannon Diversity vs. Temperature of DNA data
#####################################

#plotting a linear regression of diversity over temperature.
plot(diversity(data_relative_RNA)~Environmental_Data_RNA$Temp_manual,xlab="water temperature",ylab="shannon diversity index",main="diversity and temperature in RNA samples")
testmodel <-  lm(diversity(data_relative_RNA)~Environmental_Data_RNA$Temp_manual)
abline(testmodel$coefficients[1],testmodel$coefficients[2])

label=paste0("R_sqrd = ",round(summary(testmodel)[[8]],3),", F(",summary(testmodel)$fstatistic[[2]],", ",summary(testmodel)$fstatistic[[3]],") = ",
             round(summary(testmodel)$fstatistic[[1]],2),", p = ",format(summary(testmodel)[[4]][7],digits=2)," ",sign_stars(summary(testmodel)[[4]][7]),"\n",
             "slope = ",format(summary(testmodel)[[4]][2],digits=3),", p = ",format(summary(testmodel)[[4]][8],digits=2)," ",sign_stars(summary(testmodel)[[4]][8]))

text(mean(Environmental_Data_RNA$Temp_manual),min(diversity(data_relative_RNA))+0.1*(max(diversity(data_relative_RNA))-min(diversity(data_relative_RNA))),label, cex=0.8)
boxplot(diversity(data_relative_RNA))


#summary: diversity not correlated with sampling time
summary(lm(diversity(data_relative_RNA)~as.Date(paste0(Environmental_Data_RNA$Date,"-",Environmental_Data_RNA$Date))))

#calculate NMDS
NMDS.bray_RNA<-metaMDS(data_relative_RNA, dist="bray", k=3, trymax=100,wascores=TRUE, trace=TRUE, zero="add")
NMDS.bray_log_RNA<-metaMDS(logdata_RNA_spec, dist="bray", k=3, trymax=100,wascores=TRUE, trace=TRUE, zero="add")
NMDS.bray_sqr_RNA<-metaMDS(sqrdata_RNA_spec, dist="bray", k=3, trymax=100,wascores=TRUE, trace=TRUE, zero="add")



#with which plot dimensions and for which NMDS should be plotted?
limits <- data.frame(x=c(-0.8,1.2),y=c(-0.6,1))
usedNMDS <- NMDS.bray_RNA

#pdf(paste0("NMDS_stressplot_",datatype,"_01.pdf"),width=14,height=8)
draw_stress_plot()
plot (usedNMDS, display = 'sites', type = 't', main = 'Goodness of fit',cex=0.65) 
# and this adds the points with size reflecting goodness of fit (bigger = worse fit)
abline(v=0)
abline(h=0)
points (usedNMDS, display = 'sites', cex = goodness (usedNMDS)*200) 


#positioning the legend
legendinset=c(0.2,0)
#empty plot background
mds.fig.spe <- ordiplot(usedNMDS, type = "none",xlim=limits$x,ylim=limits$y)
#pdf(paste0("NMDS_spider_",datatype,"_01_temp_season.pdf"),width=10,height=10)
#by temperature and by season
#ordisurf line profile of temperature gradient
ordisurf(NMDS.bray_RNA, Environmental_Data_RNA$Temp_manual,main="",col="forestgreen")



#defining point colors and associated season factor levels
pointcolors <- data.frame(Season=levels(as.factor(Environmental_Data_RNA$Season)),color=c("yellow","green","red","blue")[1:4])
#plots one factor level for each iteration, in the defined colors
for(i in 1:nrow(pointcolors)){
  points(mds.fig.spe, "sites", pch = 19, col = pointcolors$color[i], select = Environmental_Data_RNA$Season == pointcolors$Season [i])
}


#connects points to their group centroids
ordispider(usedNMDS,groups=Environmental_Data_RNA$Season,col=pointcolors$color)
#draws legend of pointcolors
legend("topright", bg="white", inset=legendinset,title="Season", pointcolors$Season,lty=1, lwd=5,col=pointcolors$color,cex=0.8)
#reads stress value from model and plots it
text(labels=paste0("stress value: ",round(usedNMDS$stress,5)),x=0,y=limits$y[1],cex=0.8)
#dev.off()





