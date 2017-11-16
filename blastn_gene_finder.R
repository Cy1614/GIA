#----------------------------all function------------------------
transfer2gene <- function(df,conversion_list){
  #df<-blast_data
  #conversion_list<-protein.list
  #head(df)
  #head(df$"Query/Gene_ID")
  df$Gene_ID<-conversion_list[df$"Query/Gene_ID"]
  s <- as.numeric(df$Begin < df$End)
  s[which(s==0)] <- '-'
  s[which(s==1)] <- '+'
  df$S<-s
  ##write.table(df, 'blastn_dmoj_converted.txt',col.names=T,row.names=F,quote=F, sep='\t')
  return(df)
}

blastn_gene_finder<-function(gene_data){
  #gene_data<-blast_gene_ex_list[[2]]
  #gene_data<-blast_gene_ex_list[["FBgn0003256"]]
  gene_data_list<-split(gene_data, gene_data$Scaffold_ID)
  gene_data_list<-lapply(gene_data_list,function(x)x[order(x$Begin)])
  remove_exon_outlier<-lapply(gene_data_list,function(x){
    #x<-gene_data_list[[1]];x
    close_gene_loc<-which(diff(x$Begin)<100000)
    if(length(close_gene_loc)>0){
      #close_gene_real_loc<-sort(unique(c(sapply(close_gene_loc,function(x)c(x,x+1)))))
      Vec <- close_gene_loc
      #Vec<-1:10
      Breaks <- c(0, which(diff(Vec) != 1), length(Vec)) 
      consec_list<-lapply(seq(length(Breaks) - 1), 
                          function(i) Vec[(Breaks[i] + 1):Breaks[i+1]])
      
      #if(is.matrix(consec_list)==T){
      #real_consec_list<-list(c(consec_list,tail(consec_list,1)+1))
      #consec_list<-lapply(consec_list,list)
      # }else{
      real_consec_list<-lapply(consec_list,function(d)c(d,tail(d,1)+1))
      #}
      #as.data.frame(x[consec_list[[2]]])
      consec_score<-sapply(real_consec_list,function(u)sum(x[u]$Bit_score))
      max_consec_score_loc<-which(consec_score==max(consec_score))
      #exon_group<-lapply(consec_list[max_consec_len_loc],function(z)c(unlist(z),unlist(tail(z,1))+1))
      exon_group<-real_consec_list[max_consec_score_loc][1]
      chosen_one<-x[unlist(exon_group)]
      return(chosen_one)
    }
    
    if(length(close_gene_loc)==0){
      temp_data<-x[which(x$Bit_score==max(x$Bit_score))][1]
      return(temp_data)
    }
  }
  )
  non_0_list_loc<-which(sapply(remove_exon_outlier,function(q)length(q))!=0)
  non_0_list_data<-remove_exon_outlier[non_0_list_loc]
  non_0_bit_score<-sapply(non_0_list_data,function(t)sum(t$Bit_score))
  winner_loc<-which(non_0_bit_score==max(non_0_bit_score))
  if(length(winner_loc)>1){
    output<-do.call(rbind,non_0_list_data[winner_loc])
    if(nrow(output)>0){
      output$Gene_ID<-paste(output$Gene_ID,1:length(output$Gene_ID),sep="_")
      output$unique_loc<-paste(output$unique_loc,1:length(output$Gene_ID),sep="_")
    }
  }else{
    output<-as.data.frame(non_0_list_data[winner_loc])
  }
  colnames(output)<-c('unique_ID', 'Scaffold_ID', 'Percent_Identity', 
                      'Alignment_Length', 'Num_Mismathces', 'Num_Gap', 
                      'mRNA_Begin', 'mRNA_End', 'Begin', 
                      'End', 'E_val', 'Bit_score',"Gene_ID","S","unique_loc")
  output
}
integrated_best_gene<-function(blast_data,type){
  #type="protein"
  #blast_data <- fread('tblastn_dmoj.txt')
  names(blast_data) <- c('Query/Gene_ID', 'Scaffold_ID', 'Percent_Identity', 
                         'Alignment_Length', 'Num_Mismathces', 'Num_Gap', 
                         'Start_Query_Pos', 'End_Query_Pos', 'Begin', 
                         'End', 'E_val', 'Bit_score')
  if (type=="mRNA"){
    conversion_list<-gene.list
  }
  if (type=="protein"){
    conversion_list<-protein.list
  }

  blast_data<-transfer2gene(blast_data,conversion_list)
  blastn_neg_begin<-blast_data[which(blast_data$S=="-"),"Begin"]
  blastn_neg_end<-blast_data[which(blast_data$S=="-"),"End"]
  blast_data[which(blast_data$S=="-"),"End"]<-blastn_neg_begin
  blast_data[which(blast_data$S=="-"),"Begin"]<-blastn_neg_end
  eval_cutoff<-calculate_cutoff(blast_data$E_val,drawPlot=F)$absolute
  #eval_cutoff<-0.0001
  blast_data$unique_loc<-paste(blast_data$Scaffold_ID,blast_data$Begin,blast_data$End,blast_data$S,sep="_",blast_data$Gene_ID)
  blast_data<-blast_data[-which(duplicated(blast_data$unique_loc)),]
  passed_blast_df<-blast_data[blast_data$E_val<eval_cutoff,]
  #head(passed_blast_df)
  #head(passed_blast_df)
  blast_gene_ex_list <- split(passed_blast_df, passed_blast_df$Gene_ID)
  #blastn_gene_finder(blast_gene_ex_list[[6]])
  best_genes<-mclapply(blast_gene_ex_list,blastn_gene_finder,mc.cores = getOption("mc.cores", 6L))
  best_genes_mat<-do.call(rbind,best_genes)
  #which.max(as.numeric(best_genes_mat$End)-as.numeric(best_genes_mat$Begin))
  
  return(best_genes_mat)
}

setwd('/Users/LawCheukTing/Desktop/gia2_result/')
library(data.table)
library(GenomicRanges)
library(parallel)
library(rtracklayer)

#----------------------------pre-analysis input---------------------------------------------------
df.dmel.transcript <- readLines('dmel-all-transcript-r6.18.fasta')
df.dmel.translation <- readLines('dmel-all-translation-r6.18.fasta')

mrna.list <- gene.list <- c()
for (i in seq.int(length(df.dmel.transcript))){
  if (substr(df.dmel.transcript[i], 1, 1)=='>'){
    mrna.list <- c(mrna.list, substr(df.dmel.transcript[i], 2, 12))
    ind <- gregexpr('parent', df.dmel.transcript[i])[[1]][1]
    gene.list <- c(gene.list, substr(df.dmel.transcript[i], ind+7, ind+17))
  }
}

names(gene.list) <- mrna.list


protein_fasta_name<-df.dmel.translation[grep(">",df.dmel.translation)]
protein_gene_conversion<-function(fasta_names){
  protein.list <- substr(fasta_names, 2, 12)
  ind <- gregexpr('parent', fasta_names)[[1]][1]
  gene.list <- substr(fasta_names, ind+7, ind+17)
  c(protein.list,gene.list)
}
protein_gene_converted<-do.call(rbind,lapply(protein_fasta_name,protein_gene_conversion))

rownames(protein_gene_converted)<-protein_gene_converted[,2]
protein.list<-protein_gene_converted[,2]
names(protein.list)<-protein_gene_converted[,1]

#---------------------------------------------------analysis-------------------------------
gtf_file<-import("/Users/LawCheukTing/Desktop/gia2_result/dmoj-all-r1.04.gtf")
gene_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="gene"]
exon_gtf_file<-gtf_file[gtf_file@elementMetadata[,"type"]=="exon"]
blast_data <- fread('blastn_dvir.txt')
head(blast_data[-duplicated(blast_data$V9),][order(blast_data[-duplicated(blast_data$V9),"V11"])],100)

best_genes_data<-integrated_best_gene(blast_data,"mRNA")
blast_best_gene_ex_list <- split(best_genes_data, best_genes_data$Gene_ID)
blast_gene_region_data<-mclapply(blast_best_gene_ex_list,mc.cores = getOption("mc.cores", 6L),function(x)genscan_gene_region_finder(x))
blast_gene_region_data<-do.call(rbind,blast_gene_region_data)
blast_gene_region_df<-as.data.frame(blast_gene_region_data)
colnames(blast_gene_region_df)<-c("Scaffold_ID", "S","Begin","End","unique_ID")
#-------------------------gene overlap------------------
blast_gene_gr<-GRanges(blast_gene_region_df$Scaffold_ID,strand=blast_gene_region_df$S,IRanges(as.numeric(as.character(blast_gene_region_df$Begin)),as.numeric(as.character(blast_gene_region_df$End))),rownames(blast_gene_region_df))
blast_overlap<-findOverlaps(blast_gene_gr,gene_gtf_file)
#------------------------exon overlap-------------------
exon_gtf_df<-as.data.frame(exon_gtf_file)
unique_chr_loc<-apply(exon_gtf_df,1,function(x)paste(x[1:3],collapse = "_"))
dup_loc<-which(duplicated(unique_chr_loc))
no_dup_exon_gtf_file<-exon_gtf_file[-dup_loc]
blast_exon_gr<-GRanges(best_genes_data$Scaffold_ID,strand=best_genes_data$S,IRanges(as.numeric(as.character(best_genes_data$Begin)),as.numeric(as.character(best_genes_data$End))),as.character(best_genes_data$unique_ID))
blast_exon_overlap<-findOverlaps(blast_exon_gr,no_dup_exon_gtf_file)
#--------------------------summarise data-----------------
gene_data<-c(length(unique(blast_overlap@from)),length(unique(blast_overlap@to)));gene_data
protein_data<-c(length(unique(blast_exon_overlap@from)),length(unique(blast_exon_overlap@to)));protein_data
#-------------------------gene gtf output------------------
blastn_gene_gtf_data<-cbind(as.character(blast_gene_region_df$Scaffold_ID),"blastn","gene",as.character(blast_gene_region_df$Begin),as.character(blast_gene_region_df$End),".",as.character(blast_gene_region_df$S),".",rownames(blast_gene_region_df))
write.table(blastn_gene_gtf_data,"blastn_dvir_gene.gtf",sep="\t",row.names=F,col.names=F,quote=F)

blastn_exon_gtf_data<-cbind(as.character(best_genes_data$Scaffold_ID),"blastn","exons",as.character(best_genes_data$Begin),as.character(best_genes_data$End),".",as.character(best_genes_data$S),".",best_genes_data$Gene_ID)
write.table(blastn_exon_gtf_data,"blastn_dvir_exons.gtf",sep="\t",row.names=F,col.names=F,quote=F)


#----------------------plot genes graph-------------------
blast_plot_data<-cbind(log10(as.data.frame(blast_gene_gr[blast_overlap@from])$width),log10(as.data.frame(gene_gtf_file[blast_overlap@to])$width))

plot_col_scatter(blast_plot_data)
length(unique(blast_overlap@from))
#----------------------plot exons graph-------------------
blast_exon_plot_data<-cbind(log10(as.data.frame(blast_exon_gr[blast_exon_overlap@from])$width),log10(as.data.frame(no_dup_exon_gtf_file[blast_exon_overlap@to])$width))
plot_col_scatter(blast_exon_plot_data)



