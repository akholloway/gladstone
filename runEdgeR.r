#setwd("<PATH>/BMI206_rnaSeq/counts/")
setwd("~/Documents/teaching/BMI206_rnaSeq/BMI206_rnaSeqLab")

library(edgeR)
#edgeRUsersGuide() # this user's guide is SUPER GOOD!

allDat <- read.delim(file="test.counts.txt", header=TRUE, sep="\t", strip.white=TRUE, row.names=1, stringsAsFactors=FALSE)
allDat <- data.matrix(allDat)

# keep "expressed" genes
dat=allDat[(which(apply(allDat,1,sum)>10)),]

# sanity check histogram of read counts
par(mfrow=c(1,2),mar=c(3,3,3,3))
hist(allDat[,1],xlim=c(0,500),nclass=1000)
hist(dat[,1],xlim=c(0,500),ylim=c(0,80),nclass=1000)

# set groups
targets <- data.frame("sample"=colnames(dat),"treat"=c(rep("AD",times=3),rep("CON",times=3)))

# check targets to make sure samples are specified correctly
targets

# edgeR
group=targets$treat
design <- model.matrix(~group)
colnames(design) <- levels(group)

y <- DGEList(counts=dat,group=group)
y <- calcNormFactors(y)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

fit <- edgeR::glmFit(y,design)
lrt <- glmLRT(fit,coef=2)

topTags(lrt)

# put data together
colnames(lrt$fitted.values)=paste(colnames(lrt$fitted.values),".fit",sep="")
fin <- data.frame(dat, lrt$fitted.values, lrt$table, "Adjusted_FDR"=p.adjust(lrt$table$PValue, method="BH"))

mcFin <- apply(fin[,1:12],2,sum)
mcFin=mcFin/100000

# plots to check raw vs. fitted values cpm
#pdf("rawVSfitted_cpm.pdf",height=8,width=12)
par(mfrow=c(2,3),mar=c(5,4,3,1))
for(i in 1:6){
  plot(fin[,i]/mcFin[i],fit$fitted.values[,i]/mcFin[i],main=colnames(dat)[i],xlab="raw",ylab="fitted",pch=16)
}
#dev.off()

# log2 raw vs fitted cpm
#pdf("rawVSfitted_log2_cpm.pdf",height=8,width=12)
par(mfrow=c(2,3),mar=c(5,4,3,1))
for(i in 1:6){
  plot(log(fin[,i]/mcFin[i],2),log(fit$fitted.values[,i]/mcFin[i],2),main=colnames(dat)[i],xlab="log2 raw",ylab="log2 fitted",pch=16)
}
#dev.off()

# log2 FC vs Median WT
wtFitlog2cpm=log(apply(fit$fitted.values[,4:6]/mcFin[4:6],1,median),2)

#pdf("log2FCvsMedianWT.pdf",height=5,width=5)
par(mfrow=c(1,1),mar=c(5,4,3,1))
plot(wtFitlog2cpm,fin$logFC,main="log2FC vs median WT cpm",xlab="log2 median WT cpm",ylab="log2FC",pch=16)
points(wtFitlog2cpm[which(fin$PValue<0.05)],fin$logFC[which(fin$PValue<0.05)],pch=16,col="blue")
points(wtFitlog2cpm[which(fin$Adjusted_FDR<0.05)],fin$logFC[which(fin$Adjusted_FDR<0.5)],pch=16,col="red")
legend("right",pch=16,col=c("red","blue"),legend=c("FDR < 0.05","p-value < 0.05"))
#dev.off()

# write out results
fin=fin[order(fin$PValue),]
write.table(fin, "stats.txt", col.names=NA, row.names=T, quote=F, sep="\t")




################################
#
#  Did we find the right gene?
#
################################

#http://www.news-medical.net/news/20110629/Gene-overexpression-linked-with-Alzheimers-disease-and-Down-syndrome.aspx

