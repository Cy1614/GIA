setwd('/Users/LawCheukTing/Desktop/gia2_result/')
library('data.table')
library("parallel")
df.dmoj.blastn<-fread('blastn_dmoj.txt')
df.dmoj.tblastn <- fread('tblastn_dmoj.txt')
head(df.dmoj.tblastn)
names(df.dmoj.blastn) <- c('Query/Gene_ID', 'Target/Scaffold_ID', 'Percent_Identity', 
                           'Alignment_Length', 'Num_Mismathces', 'Num_Gap', 
                           'Start_Query_Pos', 'End_Query_Pos', 'Start_Target_Pos', 
                           'End_Target_pos', 'E_val', 'Bit_score')

df.dmel.transcript <- readLines('dmel-all-transcript-r6.18.fasta')

mrna.list <- gene.list <- c()
for (i in seq.int(length(df.dmel.transcript))){
  if (substr(df.dmel.transcript[i], 1, 1)=='>'){
    mrna.list <- c(mrna.list, substr(df.dmel.transcript[i], 2, 12))
    ind <- gregexpr('parent', df.dmel.transcript[i])[[1]][1]
    gene.list <- c(gene.list, substr(df.dmel.transcript[i], ind+7, ind+17))
  }
}

names(gene.list) <- mrna.list

transfer.bl.gc <- function(df){
  df$Gene<-gene.list[df$`Query/Gene_ID`]
  s <- as.numeric(df$End_Target_pos > df$Start_Target_Pos)
  s[which(s==0)] <- '-'
  s[which(s==1)] <- '+'
  df$S<-s
  #write.table(df, 'blastn_dmoj_converted.txt',col.names=T,row.names=F,quote=F, sep='\t')
  return(df)
}
blastn_df<-transfer.bl.gc(df.dmoj.blastn)
colnames(blastn_df) <- c('unique_ID', 'Scaffold_ID', 'Percent_Identity', 
                           'Alignment_Length', 'Num_Mismathces', 'Num_Gap', 
                           'mRNA_Begin', 'mRNA_End', 'Begin', 
                           'End', 'E_val', 'Bit_score',"Gene_ID","S")

blastn_neg_begin<-blastn_df[which(blastn_df$S=="-"),"Begin"]
blastn_neg_end<-blastn_df[which(blastn_df$S=="-"),"End"]
blastn_df[which(blastn_df$S=="-"),"End"]<-blastn_neg_begin
blastn_df[which(blastn_df$S=="-"),"Begin"]<-blastn_neg_end
blastn_df$unique_loc<-paste(blastn_df$Scaffold_ID,blastn_df$Begin,blastn_df$End,blastn_df$S,sep="_",blastn_df$Gene_ID)
blastn_df<-blastn_df[-which(duplicated(blastn_df$unique_loc)),]
calculate_cutoff(blastn_df$E_val)
passed_blastn_df<-blastn_df[blastn_df$E_val<0.000133,]

blastn_gene_ex_list <- split(passed_blastn_df, passed_blastn_df$Gene_ID)

blastn_gene_finder<-function(gene_data){
  #gene_data<-blastn_gene_ex_list[[8513]]
  gene_data_list<-split(gene_data, gene_data$Scaffold_ID)
  gene_data_list<-lapply(gene_data_list,function(x)x[order(x$Begin)])
  remove_exon_outlier<-lapply(gene_data_list,function(x){
    #x<-gene_data_list[[1]]
    close_gene_loc<-which(diff(x$Begin)<100000)
    if(length(close_gene_loc)>0){
      Vec <- close_gene_loc
      Breaks <- c(0, which(diff(Vec) != 1), length(Vec)) 
      consec_list<-sapply(seq(length(Breaks) - 1), 
             function(i) Vec[(Breaks[i] + 1):Breaks[i+1]])
      if(is.list(consec_list)==F){
        consec_list<-list(consec_list)
      }
      consec_len<-sapply(consec_list,length)
      max_consec_len_loc<-which(consec_len==max(consec_len))
      exon_group<-lapply(consec_list[max_consec_len_loc],function(z)c(z,tail(z,1)+1))
      chosen_one<-which.max(lapply(exon_group,function(y)sum(x[y]$Bit_score)))
      return(x[unlist(exon_group[chosen_one])])
      }
    if(length(close_gene_loc)==0){
        return(x[which(x$Bit_score==max(x$Bit_score))])
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
best_genes<-mclapply(blastn_gene_ex_list,blastn_gene_finder,mc.cores = getOption("mc.cores", 4L))
best_genes_mat<-do.call(rbind,best_genes)
blastn_gene_ex_list <- split(best_genes_mat, best_genes_mat$Gene_ID)
blastn_gene_region_data<-mclapply(blastn_gene_ex_list,mc.cores = getOption("mc.cores", 4L),function(x)genscan_gene_region_finder(x))
blastn_gene_region_data<-do.call(rbind,blastn_gene_region_data)
blastn_gene_region_df<-as.data.frame(blastn_gene_region_data)
head(sort(as.numeric(blastn_gene_region_data[,4])-as.numeric(blastn_gene_region_data[,3])),1000)
length(unique(rownames(blastn_gene_region_df)))

colnames(blastn_gene_region_df)<-c("scafold_ID", "strand","start","end","unique_ID")
#library(rtracklayer)
gene_gtf_file
blastn_gene_gr<-GRanges(blastn_gene_region_df$scafold_ID,strand=blastn_gene_region_df$strand,IRanges(as.numeric(as.character(blastn_gene_region_df$start)),as.numeric(as.character(blastn_gene_region_df$end))),as.character(blastn_gene_region_df$unique_ID))
head(blastn_gene_region_df)
blastn_gene_gr<-blastn_gene_gr[which(as.data.frame(blastn_gene_gr)$width<500000)]
blastn_overlap<-findOverlaps(blastn_gene_gr,gene_gtf_file)
plot(log10(as.data.frame(blastn_gene_gr[blastn_overlap@from])$width),log10(as.data.frame(gene_gtf_file[blastn_overlap@to])$width),pch=".")





lines(x = c(0,100), y = c(0,100))
hist(log10(as.data.frame(blastn_gene_gr[blastn_overlap@from])$width),breaks=1000)
tail(sort(as.data.frame(blastn_gene_gr[blastn_overlap@from])$width),1000)
blastn_gene_gr[blastn_overlap@from][which(as.data.frame(blastn_gene_gr[blastn_overlap@from])$width==30402611)]


blastn_gene_gr[blastn_overlap@from][which(as.data.frame(blastn_gene_gr[blastn_overlap@from])$width==27882551)]
FBtr0085341
head(passed_blastn_df)
passed_blastn_df[which(passed_blastn_df$unique_ID=="FBtr0085341"),]

which(names(blastn_gene_ex_list)=="FBgn0039621")

which(names(blastn_gene_ex_list)=="FBgn0039621")

head(blastn_gene_ex_list)
good_blastn_overlap<-blastn_overlap[-unwanted_loc]
plot(log10(),log10(as.data.frame(gene_gtf_file[good_blastn_overlap@to])$width),pch=".",xlim=c(2,7.5),ylim=c(2,7.5))
