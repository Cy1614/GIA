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
length(unique(exon_overlap@to))
genscan_ex_plot_data<-cbind(log10(as.data.frame(overlap_gene_ex_gr[exon_overlap@from])$width),log10(as.data.frame(no_dup_exon_gtf_file[exon_overlap@to])$width))


gene_ex_list <- split(passed_genscan_out, passed_genscan_out$unique_ID)
gene_region_data<-mclapply(gene_ex_list,mc.cores = getOption("mc.cores", 4L),function(x)genscan_gene_region_finder(x))
gene_region_mat<-do.call(rbind,gene_region_data)
colnames(gene_region_mat)<-c("Scaffold_ID","S","Begin","End","unique_ID")
genscan_gene_gtf_data<-cbind(gene_region_mat[,"Scaffold_ID"],"genscan","gene",gene_region_mat[,"Begin"],gene_region_mat[,"End"],".",gene_region_mat[,"S"],".",gene_region_mat[,"unique_ID"])

write.table(genscan_gene_gtf_data,"/Users/LawCheukTing/Desktop/gia2_result/genscan_dmoj_gene.gtf",sep="\t",quote=F,col.names = F,row.names = F)

genscan_exon_gtf_data<-cbind(passed_genscan_out$Scaffold_ID,"genscan",passed_genscan_out$Type,passed_genscan_out$Begin,passed_genscan_out$End,".",passed_genscan_out$S,".",passed_genscan_out$unique_ID)
write.table(genscan_exon_gtf_data,"/Users/LawCheukTing/Desktop/gia2_result/dmoj_exon.gtf",sep="\t",quote=F,col.names = F,row.names = F)

pure_gene_out[which(pure_gene_out$unique_ID=="scaffold_14906_5410"),]
?countOverlaps
which(overlap_gene_ex_gr@seqnames=="scaffold_6540")
