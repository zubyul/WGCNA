#__setup__#
getwd()
workingDir = "C:/Users/Administrator/Downloads/FemaleLiver-Data" 
setwd(workingDir)
library(WGCNA) #load WGCNA
options(stringsAsFactors = FALSE) #factors are invariable, do not use

lnames = load(file = "FemaleLiver-01-dataInput.RData");
#loads datExpr and datTraits from previous run

#__soft thresholding__#
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#choosing a set of soft thresholding powers:create powers
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#genetates a list with powerestimate and fitindices
#"An index of fit is a catch-all term for a variety of methods to tell you how
##well observed data fits a particular probability distribution.
##An index of fit is typically normalized (i.e. units of measurement are
##removed), and the values will usually be between 0 and 1"
#function:using powers candidates analyzes data. verbose to show background activity while processing command
#par(mfrow = c(1,2));
#cex1 = 0.9;
#character expansion factor in order to account for mfrow
#don't know why so I'm removing it
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex = NULL,col="red");
#fit indices as funct of soft-thresholding
abline(h=0.90,col="red")
#establishes r^2 cutoff
#"The left panel shows the scale-free fit index (y-axis) as a function of the soft-thresholding power (x-axis)."
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# "The right panel displays the _mean connectivity(??)_ (degree, y-axis) as a function of the soft-thresholding power (x-axis)".

#you're trying to have least scale effects while still maintaining mean connectivity

#__Network construction and module detection__#
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
# power - 6; min. module size- 30; medium cluster split sensitivity (auto);mergeCutHeight module merging threshold;
#output numeric not color; save topo overlap matrix
##maxBlockSize tells the function how large the largest block can be that the reader's computer can handle
#t if this code were to be used to analyze a data set with more than 5000 probes, the function blockwiseModules will split the data set into several blocks. This will break some of the plotting code below

table(net$colors)
#will produce table of modules in order of descending size (except 0 which is genes with no modules)

#__plotting__#
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
