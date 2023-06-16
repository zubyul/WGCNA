#___INPUT___#
getwd()
workingDir = "C:/Users/Administrator/Downloads/FemaleLiver-Data" 
#here/downloads; make a string to where the file is

setwd(workingDir)
library(WGCNA) #load WGCNA

options(stringsAsFactors = FALSE)
#factors are invariable, do not use
femData = read.csv("LiverFemale3600.csv") #read CSV
#suggested to examine head, dim, names, etc

datExpr0 = as.data.frame(t(femData[, -c(1:8)]));
#matrix -> data only(transposed/sideways (minus first data columns))
names(datExpr0) = femData$substanceBXH; 
rownames(datExpr0) = names(femData)[-c(1:8)];
#switch names and columns

#__quality__#
gsg = goodSamplesGenes(datExpr0, verbose = 3);
#filters samples with too many missing entries
#verbosity seems to change output column width?
gsg$allOK #are they all all ok?

if(!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
#print the ones being removed & remove

plot(hclust(dist(datExpr0), method ="average"))
sampleTree = hclust(dist(datExpr0), method = "average");
#heirarchical structure (distance size of file, UPGMA method)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
#removed parameters for the plot from the source 
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#title and params

abline(h = 15, col = "red")
#cutting outliers above the line

clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
#cut tree at height 15 to remove branches without min obj # of 10
#just running it gives a table of clusters with branch # below
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#keep


#___CLINICAL TRAIT DATA__#
traitData = read.csv("ClinicalTraits.csv");
dim(traitData)
names(traitData)
#load
allTraits = traitData[, -c(31, 16)];
allTraits = allTraits[, c(2, 11:36) ];
#remove columns that hold information we don't need
dim(allTraits)
names(allTraits)

femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Mice);
#returns a vector of the positions of (first) matches of its first argument in its second.
#matches the mouse names ex. F3_320

datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
#rearrange and name
collectGarbage();
#??


sampleTree2 = hclust(dist(datExpr), method = "average")

traitColors = numbers2colors(datTraits, signed = FALSE);
# Convert traits to a color representation: white means low, red means high, grey means missing entry.
# reflects quantity of biological data about the rats beyond experimental data
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
# Plot the sample dendrogram and the colors underneath.
save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")
