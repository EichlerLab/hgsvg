library(getopt)
options <- matrix(c("sv", "v", 2, "character",
                    "count", "c", 2, "character",
                    "operation", "o", "op", "character",
                    "sample", "s", 2, "character"), byrow=T, ncol=4)
LABELSIZE=5
args <- getopt(options)
                                        #
#setwd("/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/IlluminaCombined/HG00514")
#args <- data.frame(sv="int_caller_full.INS.bed",count="callers.INS.tab", sample="HG00514", operation="INS")

callTab   <- read.table(as.character(args$sv),header=T,comment.char="")
callNames <- names(callTab)
firstName <- which(callNames == "tEnd.1")[1]+1
lastName  <- which(callNames == "orth_filter.1")[1]-1

countTab  <- read.table(as.character(args$count),header=T,comment.char="")
nSamples  <- lastName - firstName + 1

countNames <- names(countTab)
firstCountName <- which(countNames == "tEnd")[1]+1
lastCountName <- length(countNames)

allCounts <- apply(countTab[,firstCountName: lastCountName],2,sum)


#                                      
# First pca of all values
tpca <- prcomp(t(callTab[,firstName:lastName]), center=T, scale=T)

tabNames <- names(callTab[,firstName:lastName])

# 
# Counts table

tabAllCounts <- sapply(tabNames, function(n) if (is.na(allCounts[n])) { return(0) } else { return(allCounts[n]) } )

library(gdsfmt)
library(ggrepel)

n <- length(colnames(tpca$x))
colnames(tpca$x) <-  paste("c",seq(1,n),sep="")
pcaDF <- as.data.frame(tpca$x)

pcaSummary <- summary(tpca)
require(gridExtra)
pdf(sprintf("MethodPCA.%s.%s.pdf", args$operation, args$sample,sep=""), width=12,height=6)

p1 <- ggplot(pcaDF, aes(x=c1,y=c2)) + geom_point(color = 'black') + geom_text_repel(aes(label = tabNames),size=LABELSIZE)  + xlab(sprintf("PC 1 %2.2f%% variance ",100*pcaSummary$importance[2,1]) ) + ylab(sprintf("PC 2 %2.2f%% variance ", 100*pcaSummary$importance[2,2])) + labs(title=sprintf("%s, Illumina combined %s", args$sample, args$operation)) + theme_bw()

p2 <- ggplot(pcaDF, aes(x=c2,y=c3)) + geom_point(color = 'black') + geom_text_repel(aes(label = tabNames),size=LABELSIZE)  + xlab(sprintf("PC 2 %2.2f%% variance ",100*pcaSummary$importance[2,2]) ) + ylab(sprintf("PC 3 %2.2f%% variance ", 100*pcaSummary$importance[2,3])) + labs(title=sprintf("%s, Illumina combined %s", args$sample, args$operation)) + theme_bw()

grid.arrange(p1,p2,ncol=2)

dev.off()

print("done plotting pca")
library(lsa)
library(RColorBrewer)
library(lattice)
library(proxy)
tabPass <- callTab[which(callTab$orth_filter == "PASS"),]
tabFail <- callTab[which(callTab$orth_filter == "FAIL"),]

cmat <- cosine(as.matrix(tabPass[,firstName:lastName]))
jmat <- dist(t(tabPass[,firstName:lastName]), method="Jaccard", pairwise=T)
reds <- brewer.pal(9, "RdBu")
cr <- colorRamp(reds)
rampCol <- rgb(cr(seq(0,1,by=0.01))/255)
#levelplot(as.matrix(jmat), col.regions=rampCol)
#heatmap(as.matrix(jmat), symm=T, col=rampCol)


jmat <- dist(t(callTab[,firstName:lastName]), method="Jaccard", pairwise=T)
#reds <- brewer.pal(9, "Blues")
cr <- colorRamp(reds)
n <- dim(jmat)[1]
cl <- hclust(dist(as.matrix(jmat)))
pdf(sprintf("Jaccard.%s.%s.pdf",args$operation, args$sample))
levelplot(as.matrix(jmat)[cl$order,cl$order], col.regions=rampCol,scales=list(x=list(rot=90)), xlab="", ylab="", main=sprintf("Jaccard similarity %s %s",args$sample, args$operation) )
dev.off()
            
#hm <- heatmap(as.matrix(jmat), symm=T, col=rampCol,plot=F)


passSum <- apply(tabPass[,firstName:lastName],2,sum)
callsSum <- apply(callTab[,firstName:lastName],2,sum)

pdf(sprintf("MethodsBar.%s.%s.pdf",args$operation, args$sample),width=8,height=4)
bpx <- barplot(rbind(passSum, callsSum-passSum), names.arg=tabNames,col=c("black","red"), main=sprintf("%s %s", args$sample, args$operation), xaxt="n")
mc <- max(callsSum)
text(cex=1, x=bpx-.25, y=-mc*0.15, tabNames, xpd=TRUE, srt=45,pos=1)

legend("topright", legend=c("Confirmed", "Unconfirmed"), pch=22, pt.bg=c("black","red"), pt.cex=2)
dev.off()

passTab <- rbind(passSum, callsSum, 100*passSum/callsSum)
rownames(passTab) <- c("confirmed", "total", "fraction")
methodSummary <- sprintf("MethodSummary.%s.%s.tsv",args$operation,args$sample)

write.table(passTab, methodSummary, sep="\t", quote=F)
#apply(tab[,4:16],2,length)
#
#
#tabins <- read.table("int_margin.INS.bed",header=T,comment.char="")
#
#tpca <- prcomp(t(tab[,4:16]), center=T, scale=T)
#names(tpca)
#library(gdsfmt)
#dim(tpca$x)
#library(ggrepel)
#plot(tpca$x[,1], tpca$x[,2])
#
#names(tab)
#df <- as.data.frame(tpca$x)
#
#
#colnames(tpca$x) <-  paste("c",seq(1,13),sep="")
#
#
#
#names <- colnames(tab)[4:16]
#
#tpca <- prcomp(t(tabPass[,4:16]), center=T, scale=T)
#df <- as.data.frame(tpca$x)
#pdf("MethodPCA.NA1940.pass.pdf")
#
#ggplot(df, aes(x=PC1,y=PC2)) + 
#  geom_point(color = 'black') + 
#  geom_text_repel(aes(label = colnames(tab)[4:16]))  + xlab(sprintf("PC 1 %2.2f%% variance ",100*s$importance[2][1]) ) + ylab(sprintf("PC 2 %2.2f%% variance ", 100*s$importance[2][1])) + labs(title="NA19240, Filtered")
#de.vooff()
#
