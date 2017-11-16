genscan_gene_region_finder<-function(gene_ex_mat){
  #gene_ex_mat<-gene_ex_list[[1]]
  gene_start<-min(gene_ex_mat$Begin)
  gene_end<-max(gene_ex_mat$End)
  c(gene_ex_mat$Scaffold_ID[1],gene_ex_mat$S[1],gene_start,gene_end,gene_ex_mat$unique_ID[1])
}
library(parallel)
pure_gene_out<-fread("/Users/LawCheukTing/Desktop/gia2_result/genscan_dmoj_clean.txt")
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
gtf_file<-import("/Users/LawCheukTing/Desktop/gia2_result/dmoj-all-r1.04.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="gene"]
cm_gene_region_df<-gene_region_df[which(gene_region_df$scafold_ID%in%gene_gtf_file@seqnames),]
genscan_gene_region_gr<-GRanges(cm_gene_region_df$scafold_ID,strand=cm_gene_region_df$strand,IRanges(as.numeric(cm_gene_region_df$start),as.numeric(cm_gene_region_df$end)),cm_gene_region_df$unique_ID)
#tiny_genscan_loc<-which(as.data.frame(genscan_gene_region_gr)$width<100)
#genscan_gene_region_gr<-genscan_gene_region_gr[-tiny_genscan_loc]
genscan_gene_overlap<-findOverlaps(genscan_gene_region_gr,gene_gtf_file)
sum(table(genscan_gene_overlap@to)[which(table(genscan_gene_overlap@to)>1)])
#------------the difference in gene length-------------------------
plot(density(log10(as.data.frame(gene_gtf_file)$width)),col="blue")
lines(density(log10(as.data.frame(genscan_gene_region_gr)$width)),col="red")

plot(log10(as.data.frame(genscan_gene_region_gr[genscan_gene_overlap@from])$width),log10(as.data.frame(gene_gtf_file[genscan_gene_overlap@to])$width),pch=".",xlim=c(2,5.5),ylim=c(2,5.5))


genescan_plot_data<-cbind(log10(as.data.frame(gene_gtf_file[genscan_gene_overlap@to])$width),log10(as.data.frame(genscan_gene_region_gr[genscan_gene_overlap@from])$width))


plot(density(genescan_plot_data[,2]),col="blue")


diamonds

ggplot(diamonds, aes(depth, colour = cut)) +
  geom_density()
