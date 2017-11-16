library(rtracklayer)
plot_data<-c()
gtf_data<-c()
file_mat<-matrix(c("genscan_dmoj_clean.txt","genscan_dgri_clean.txt","genscan_dvir_clean.txt","dmoj-all-r1.04.gtf","dgri-all-r1.05.gtf","dvir-all-r1.06.gtf","dmoj","dgri","dvir"),,3)
for( i in 1:3){
  #i<-1
  pure_gene_out<-fread(file_mat[i,1])
  neg_begin<-pure_gene_out[pure_gene_out$S=="-"]$Begin
  neg_end<-pure_gene_out[pure_gene_out$S=="-"]$End
  pure_gene_out[pure_gene_out$S=="-"]$End<-neg_begin
  pure_gene_out[pure_gene_out$S=="-"]$Begin<-neg_end
  promoter_loc<-grep("Prom",pure_gene_out$Type)
  pure_gene_out<-pure_gene_out[-promoter_loc]
  gene_Ex_tab<-read.table(text=sprintf("%.4f", pure_gene_out$Gn.Ex),sep=".")
  pure_gene_out$unique_ID<-paste(pure_gene_out$Scaffold_ID,gene_Ex_tab[,1],sep="_")
  passed_genscan_out<-pure_gene_out[pure_gene_out$P>0.5]
  #passed_genscan_out<-pure_gene_out
  gene_ex_list <- split(passed_genscan_out, passed_genscan_out$unique_ID)
  gene_region_data<-mclapply(gene_ex_list,mc.cores = getOption("mc.cores", 4L),function(x)genscan_gene_region_finder(x))
  gene_region_mat<-do.call(rbind,gene_region_data)
  rownames(gene_region_mat)<-names(gene_region_data)
  gene_region_df<-as.data.frame(as.matrix(gene_region_mat),stringsAsFactors=F)
  colnames(gene_region_df)<-c("scafold_ID", "strand","start","end","unique_ID")
  #gtf_file<-import(file_mat[i,2])
  #gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="gene"]
  #cm_gene_region_df<-gene_region_df[which(gene_region_df$scafold_ID%in%gene_gtf_file@seqnames),]
  #genscan_gene_region_gr<-GRanges(cm_gene_region_df$scafold_ID,strand=cm_gene_region_df$strand,IRanges(as.numeric(cm_gene_region_df$start),as.numeric(cm_gene_region_df$end)),cm_gene_region_df$unique_ID)
  #genscan_gene_overlap<-findOverlaps(genscan_gene_region_gr,gene_gtf_file)
  #write.table(as.data.frame(genscan_gene_region_gr[genscan_gene_overlap@from]),paste(file_mat[i,3],"overlap",".txt",sep="_"),quote=F,col.names = T,row.names=F,sep="\t")
  write.table(gene_region_df,paste(file_mat[i,3],"original",".txt",sep="_"),quote=F,col.names = T,row.names=F,sep="\t")
}
head(plot_data)
colnames(plot_data)<-c("genscan","spec")
colnames(gtf_data)<-c("gtf","spec")
plot(density(log10(as.numeric(gene_region_df$end)-as.numeric(gene_region_df$start))))



dev.list()
pdf("fb_genscan_original_length_distirbution.pdf",height=6,width=8)
dev.list()
gtf_file<-import("dmel-all-r6.18.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="gene"]
plot(density(log10(as.data.frame(gene_gtf_file)$width)),col="blue",ylim=c(0,0.9),xlim=c(1.5,5.5),main="Distirbution of gene length",xlab="log10(gene length)")
gtf_file<-import("dmoj-all-r1.04.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="gene"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="turquoise1")
gtf_file<-import("dgri-all-r1.05.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="gene"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="steelblue1")
gtf_file<-import("dvir-all-r1.06.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="gene"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="navy")

#dmoj<-read.delim("dmoj_original_.txt",stringsAsFactors = F)
#lines(density(log10(dmoj$end-dmoj$start)),col="magenta1")
#dgri<-read.delim("dgri_original_.txt",stringsAsFactors = F)
#lines(density(log10(dgri$end-dgri$start)),col="mediumvioletred")
#dvir<-read.delim("dvir_original_.txt",stringsAsFactors = F)
#lines(density(log10(dvir$end-dvir$start)),col="red")




gtf_file<-import("dmel-all-r6.18.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
plot(density(log10(as.data.frame(gene_gtf_file)$width)),col="blue",ylim=c(0,2),xlim=c(0.5,4.5),main="Distirbution of gene length",xlab="log10(gene length)")
gtf_file<-import("dmoj-all-r1.04.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="turquoise1")
gtf_file<-import("dgri-all-r1.05.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="steelblue1")
gtf_file<-import("dvir-all-r1.06.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="navy")

dmoj<-read.delim("genscan_dmoj_clean.txt",stringsAsFactors = F)
dmoj<-dmoj[dmoj$P>0.5&!is.na(dmoj$P),]
lines(density(log10(abs(dmoj$End-dmoj$Begin))),col="magenta1")

dgri<-read.delim("genscan_dgri_clean.txt",stringsAsFactors = F)
dgri<-dgri[dgri$P>0.5&!is.na(dgri$P),]
lines(density(log10(abs(dgri$End-dgri$Begin))),col="mediumvioletred")
dvir<-read.delim("genscan_dvir_clean.txt",stringsAsFactors = F)
dvir<-dvir[dvir$P>0.5&!is.na(dvir$P),]
lines(density(log10(abs(dvir$End-dvir$Begin))),col="red")

legend("topright", legend=c("fb dmel","fb dmoj","fb dgri","fb dvir", "genscan dmoj","genscan dgri", "genscan dvir"),col=c("blue", "turquoise1","steelblue1","navy","magenta1","mediumvioletred","red"), lty=1, cex=0.8)
dev.off()




#-----------------------exon--------------------
file_mat<-matrix(c("genscan_dmoj_clean.txt","genscan_dgri_clean.txt","genscan_dvir_clean.txt","dmoj-all-r1.04.gtf","dgri-all-r1.05.gtf","dvir-all-r1.06.gtf","dmoj","dgri","dvir"),,3)
i<-1
for(i in 1:3){
  pure_gene_out<-fread(file_mat[i,1])
  neg_begin<-pure_gene_out[pure_gene_out$S=="-"]$Begin
  neg_end<-pure_gene_out[pure_gene_out$S=="-"]$End
  pure_gene_out[pure_gene_out$S=="-"]$End<-neg_begin
  pure_gene_out[pure_gene_out$S=="-"]$Begin<-neg_end
  promoter_loc<-grep("Prom",pure_gene_out$Type)
  pure_gene_out<-pure_gene_out[-promoter_loc]
  gene_Ex_tab<-read.table(text=sprintf("%.4f", pure_gene_out$Gn.Ex),sep=".")
  pure_gene_out$unique_ID<-paste(pure_gene_out$Scaffold_ID,gene_Ex_tab[,1],sep="_")
  passed_genscan_out<-pure_gene_out[pure_gene_out$P>0.5]
  gene_ex_list <- split(passed_genscan_out, passed_genscan_out$unique_ID)
  gene_region_data<-mclapply(gene_ex_list,mc.cores = getOption("mc.cores", 4L),function(x)genscan_gene_region_finder(x))
  gene_region_mat<-do.call(rbind,gene_region_data)
  rownames(gene_region_mat)<-names(gene_region_data)
  gene_region_df<-as.data.frame(as.matrix(gene_region_mat),stringsAsFactors=F)
  colnames(gene_region_df)<-c("scafold_ID", "strand","start","end","unique_ID")
  
  gtf_file<-import(file_mat[i,2])
  gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="gene"]
  cm_gene_region_df<-gene_region_df[which(gene_region_df$scafold_ID%in%gene_gtf_file@seqnames),]
  genscan_gene_region_gr<-GRanges(cm_gene_region_df$scafold_ID,strand=cm_gene_region_df$strand,IRanges(as.numeric(cm_gene_region_df$start),as.numeric(cm_gene_region_df$end)),cm_gene_region_df$unique_ID)
  genscan_gene_overlap<-findOverlaps(genscan_gene_region_gr,gene_gtf_file)
  overlap_gene_ID<-genscan_gene_region_gr[genscan_gene_overlap@from]$`cm_gene_region_df$unique_ID`
  overlap_gene_loc<-which(passed_genscan_out$unique_ID%in%overlap_gene_ID)
  overlap_gene_data<-passed_genscan_out[overlap_gene_loc,]
  overlap_gene_ex_gr<-GRanges(overlap_gene_data$Scaffold_ID,strand=overlap_gene_data$S,IRanges(overlap_gene_data$Begin,overlap_gene_data$End),overlap_gene_data$unique_ID)
  
  exon_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
  exon_gtf_df<-as.data.frame(exon_gtf_file)
  unique_chr_loc<-apply(exon_gtf_df,1,function(x)paste(x[1:3],collapse = "_"))
  dup_loc<-which(duplicated(unique_chr_loc))
  no_dup_exon_gtf_file<-exon_gtf_file[-dup_loc]
  exon_overlap<-findOverlaps(overlap_gene_ex_gr,no_dup_exon_gtf_file)
  write.table(as.data.frame(overlap_gene_ex_gr[exon_overlap@from]),paste(file_mat[i,3],"exon","overlap",".txt",sep="_"),quote=F,col.names = T,row.names=F,sep="\t")
}



dev.list()
pdf("fb_genscan_exon_length_distirbution.pdf",height=6,width=8)
dev.list()
gtf_file<-import("dmel-all-r6.18.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
plot(density(log10(as.data.frame(gene_gtf_file)$width)),col="blue",ylim=c(0,1.3),xlim=c(1,4),main="Distirbution of exon length",xlab="log10(exon length)")
gtf_file<-import("dmoj-all-r1.04.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="turquoise1")
gtf_file<-import("dgri-all-r1.05.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="steelblue1")
gtf_file<-import("dvir-all-r1.06.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="navy")


dmoj<-read.delim("dmoj_exon_overlap_.txt",stringsAsFactors = F)
lines(density(log10(dmoj$width)),col="magenta1")
dgri<-read.delim("dgri_exon_overlap_.txt",stringsAsFactors = F)
lines(density(log10(dgri$width)),col="mediumvioletred")
dvir<-read.delim("dvir_exon_overlap_.txt",stringsAsFactors = F)
lines(density(log10(dvir$width)),col="red")


legend("topright", legend=c("fb dmel","fb dmoj","fb dgri","fb dvir", "genscan dmoj","genscan dgri", "genscan dvir"),
       col=c("blue", "turquoise1","steelblue1","navy","magenta1","mediumvioletred","red"), lty=1, cex=0.8)
dev.off()








#--------------------blastn_length--------------
dev.list()
pdf("fb_blastn_gene_length_distirbution.pdf",height=6,width=8)
dev.list()

gtf_file<-import("dmel-all-r6.18.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="gene"]
plot(density(log10(as.data.frame(gene_gtf_file)$width)),col="blue",ylim=c(0,0.9),xlim=c(1.5,5.5),main="Distirbution of gene length",xlab="log10(gene length)")
gtf_file<-import("dmoj-all-r1.04.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="gene"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="turquoise1")
gtf_file<-import("dgri-all-r1.05.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="gene"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="steelblue1")
gtf_file<-import("dvir-all-r1.06.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="gene"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="navy")

dmoj<-fread("/Users/LawCheukTing/Desktop/gia2_result/blastn_dmoj_gene.gtf")
lines(density(log10(dmoj$V5-dmoj$V4)),col="magenta1")
dgri<-fread("/Users/LawCheukTing/Desktop/gia2_result/blastn_dgri_gene.gtf")
lines(density(log10(dgri$V5-dgri$V4)),col="mediumvioletred")
dvir<-fread("/Users/LawCheukTing/Desktop/gia2_result/blastn_dvir_gene.gtf")
lines(density(log10(dvir$V5-dvir$V4)),col="red")

legend("topright", legend=c("fb dmel","fb dmoj","fb dgri","fb dvir", "blastn dmoj","blastn dgri", "blastn dvir"),
       col=c("blue", "turquoise1","steelblue1","navy","magenta1","mediumvioletred","red"), lty=1, cex=0.8)
dev.off()



#--------------------blastn_exon_length--------------
dev.list()
pdf("fb_blastn_exon_length_distirbution.pdf",height=6,width=8)
dev.list()

gtf_file<-import("dmel-all-r6.18.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
plot(density(log10(as.data.frame(gene_gtf_file)$width)),col="blue",ylim=c(0,1.3),xlim=c(1,4),main="Distirbution of exon length",xlab="log10(exon length)")
gtf_file<-import("dmoj-all-r1.04.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="turquoise1")
gtf_file<-import("dgri-all-r1.05.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="steelblue1")
gtf_file<-import("dvir-all-r1.06.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
lines(density(log10(as.data.frame(gene_gtf_file)$width)),col="navy")


dmoj<-fread("/Users/LawCheukTing/Desktop/gia2_result/blastn_dmoj_exons.gtf")
lines(density(log10(dmoj$V5-dmoj$V4)),col="magenta1")
dgri<-fread("/Users/LawCheukTing/Desktop/gia2_result/blastn_dgri_exons.gtf")
lines(density(log10(dgri$V5-dgri$V4)),col="mediumvioletred")
dvir<-fread("/Users/LawCheukTing/Desktop/gia2_result/blastn_dvir_exons.gtf")
lines(density(log10(dvir$V5-dvir$V4)),col="red")

legend("topright", legend=c("fb dmel","fb dmoj","fb dgri","fb dvir", "blastn dmoj","blastn dgri", "blastn dvir"),
       col=c("blue", "turquoise1","steelblue1","navy","magenta1","mediumvioletred","red"), lty=1, cex=0.8)
dev.off()