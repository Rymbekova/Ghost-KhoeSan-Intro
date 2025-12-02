## module load R; R --vanilla

tedir<-"/lisc/scratch/admixlab/mktest/"

## preparation
library(biomaRt)
## get all gene IDs
hostname="grch38.ensembl.org";datasetname="hsapiens_gene_ensembl"
enscall = useEnsembl("genes", dataset=datasetname)
allgenes <- biomaRt::getBM(
  attributes = c("entrezgene_id", "external_gene_name","ensembl_gene_id"),
  mart = enscall
)
write.table(allgenes,paste(tedir,"allgeneids.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

# GOs in desired format
gos<-read.table(paste(tedir,"hu_gene2go.gz",sep=""),sep="\t")
gos<-merge(gos,allgenes,by.x=2,by.y=1,all.x=F,all.y=F)
goc<-unique(gos[,3])
allgos<-list()
for (ggo in goc) { 
  gsu<-gos[which(gos[,3]==ggo),,drop=F]
  allgos[[ggo]]<-c(ggo,paste(gsu[1,5],gsu[1,6]),paste(unique(gsu[,10]),collapse=" ")) }
allgos<-do.call(rbind,allgos)
write.table(allgos,paste(tedir,"hu_annotation.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

# gtf in desired format
gtf<-read.table(paste(tedir,"gencode.v49.exons.gtf.gz",sep=""),sep="\t")
gtf[,9]<-gsub("\\.[1-9];",";",gtf[,9])
gs<-do.call(rbind,strsplit(gtf[,9],split=";"))[,1]
gs<-do.call(rbind,strsplit(gs,split=" "))
gs<-paste(gs[,1],' "',gs[,2],'";',sep="")
gtf2<-cbind(gtf[,1:8],gs)
write.table(gtf2,paste(tedir,"gencode.v49.exons.cor.gtf",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

q()
