## module load R; R --vanilla

method=as.character(unlist((commandArgs(TRUE))))
print(method)

if (method=="getdata") {

  dadir<-"/lisc/scratch/admixlab/mk_data/san/"
  library("GenomicRanges")
  options("scipen"=100)
  
  # load fragments & individual list
  alf<-read.table(paste(dadir,"allfrag.bed.gz",sep=""),sep="\t")
  alf<-GRanges( seqnames=alf[,1], ranges=IRanges(start=alf[,2],end=alf[,3]),mcols=alf[,5])
  ksp<-unlist(read.table(paste(dadir,"ksp.txt",sep=""),sep="\t"))[-1]
  
  # load annotated variants (missense and extended set of functional sites)
  anno<-read.table(paste(dadir,"allmissense.txt",sep=""),sep="\t")
  anno<-GRanges( seqnames=anno[,1], ranges=IRanges(start=anno[,3],end=anno[,3]),mcols=anno[,12])
  anno2<-read.table(paste(dadir,"allsimp.txt",sep=""),sep="\t")
  annox<-anno2[which(anno2[,5]%in%c("3_prime_UTR_variant","5_prime_UTR_variant","mature_miRNA_variant","missense_variant","missense_variant,splice_region_variant","protein_altering_variant,incomplete_terminal_codon_variant","splice_acceptor_variant","splice_acceptor_variant,non_coding_transcript_variant","splice_donor_5th_base_variant,intron_variant","splice_donor_5th_base_variant,intron_variant,non_coding_transcript_variant","splice_donor_region_variant,intron_variant","splice_donor_region_variant,intron_variant,non_coding_transcript_variant","splice_donor_variant","splice_donor_variant,non_coding_transcript_variant","splice_polypyrimidine_tract_variant,intron_variant","splice_polypyrimidine_tract_variant,intron_variant,non_coding_transcript_variant","splice_region_variant,3_prime_UTR_variant","splice_region_variant,5_prime_UTR_variant","splice_region_variant,intron_variant","splice_region_variant,intron_variant,non_coding_transcript_variant","splice_region_variant,non_coding_transcript_exon_variant","splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant","splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant,non_coding_transcript_variant","splice_region_variant,synonymous_variant","start_lost","start_lost,splice_region_variant","stop_gained","stop_gained,splice_region_variant","stop_lost","stop_lost,splice_region_variant","stop_retained_variant")),]
  annox<-GRanges( seqnames=annox[,1], ranges=IRanges(start=annox[,3],end=annox[,3]),mcols=annox[,-c(1:3)])
  
  
  # get genes with mutations: either nonsyn (genwithmut) or broader functional (genwithfun)
  ## this loops over the chromosomes and looks up each individual whether it carries the alternative state at each SNV
  ## here I use GRanges to make things fast
  genwithmut<-list()
  genwithfun<-list()
  for (chrom in c(1:22,"X")) { 
    print(chrom)
    ct<-read.table(paste(dadir,"chr",chrom,".gt.gz",sep=""),sep=" ",header=F)[,c(1,3:27)]
    colnames(ct)<-c("pos",ksp)
    ct<-GRanges( seqnames=rep(paste("chr",chrom,sep=""),nrow(ct)), ranges=IRanges(start=ct[,1],end=ct[,1]),mcols=ct[,-1])
    poforch<-list()
    for (ind in ksp) {
      sub<-alf[grep(ind,alf@elementMetadata[,1]),]
      sub<-ct[overlapsAny(ct,sub)]
      sub<-sub[,grep(ind,colnames(sub@elementMetadata))]
      poforch[[ind]]<-sub[grep("0/1|1/1",sub@elementMetadata[,1]),]
    }
    poforch<-reduce(unlist(as(poforch, "GRangesList")))
    panno<-anno[overlapsAny(anno,poforch)]
    genwithmut[[chrom]]<-reduce(panno)
    panno<-annox[overlapsAny(annox,poforch)]
    genwithfun[[chrom]]<-reduce(panno)
  }
  
  # extract the position information for sites and write into table for GOWINDA
  genwithmut2<-as.data.frame(as(genwithmut, "GRangesList"))[,c(3,5)]
  genwithfun2<-as.data.frame(as(genwithfun, "GRangesList"))[,c(3,5)]
  
  write.table(genwithmut2,paste(dadir,"genwithmut.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
  write.table(genwithfun2,paste(dadir,"genwithfun.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

  q()
  }

if (method=="getsig") {
  
  ### second step: after results are there, convert ensg-IDs into Gene Names and get only the significant ones
  options("scipen"=100)
  
  dadir<-"/lisc/scratch/admixlab/mk_data/san/"
  tedir<-"/lisc/scratch/admixlab/mktest/"
  allgenes<-read.table(paste(tedir,"allgeneids.txt",sep=""),sep="\t",header=F,as.is=T)
  getgene<-function(input) { ip=toupper(unlist(strsplit(input,split=",")));op=paste(allgenes[which(allgenes[,3]%in%ip),2],collapse=","); return(op) }
  
  resu<-read.table(paste(dadir,"results_go_gene.txt",sep=""),sep="\t",header=F,as.is=T)
  # get significant categories
  sigres<-resu[which(resu[,5]<0.05),]
  # add a column with gene names
  sigres<-cbind(sigres,unlist(lapply(sigres[,10],getgene)))
  # save the output
  colnames(sigres)<-c("GO_ID","avrg_obs","true_obs","p_val","p_val_cor","uniq_obs","total_obs","total_gene","Category","Ensembl_ID","Gene_name")
  write.table(sigres,file=paste(dadir,"results_go_gene_name.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)

  q()
  }
