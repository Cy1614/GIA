library('data.table')
library("parallel")
df.dmoj.blastn<-fread('blastn_dmoj.txt')
df.dmoj.tblastn <- fread('tblastn_dmoj.txt')
df.dvir.blastn <- fread('blastn_dvir.txt')
df.dvir.tblastn <- fread('tblastn_dvir.txt')
df.dgri.blastn <- fread('blastn_dgri.txt')
df.dgri.tblastn <- fread('tblastn_dgri.txt')

names(df.dmoj.blastn) <- c('Query/Gene_ID', 'Target/Scaffold_ID', 'Percent_Identity', 
                           'Alignment_Length', 'Num_Mismathces', 'Num_Gap', 
                           'Start_Query_Pos', 'End_Query_Pos', 'Start_Target_Pos', 
                           'End_Target_pos', 'E_val', 'Bit_score')
names(df.dvir.blastn) <- c('Query/Gene_ID', 'Target/Scaffold_ID', 'Percent_Identity', 
                           'Alignment_Length', 'Num_Mismathces', 'Num_Gap', 
                           'Start_Query_Pos', 'End_Query_Pos', 'Start_Target_Pos', 
                           'End_Target_pos', 'E_val', 'Bit_score')
names(df.dgri.blastn) <- c('Query/Gene_ID', 'Target/Scaffold_ID', 'Percent_Identity', 
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
  #df<-df.dmoj.blastn
  df$Gene<-gene.list[df$`Query/Gene_ID`]
  s <- as.numeric(df$End_Target_pos > df$Start_Target_Pos)
  s[which(s==0)] <- '-'
  s[which(s==1)] <- '+'
  df$S<-s
  df <- data.frame(cbind(df$Gene, df$S, df$`Target/Scaffold_ID`, df$Alignment_Length, 
                         df$Start_Target_Pos, df$End_Target_pos, df$E_val, df$Bit_score))
  return(df)
}

e_val_cutoff.dmoj <- calculate_cutoff(df.dmoj.blastn$E_val)
e_val_cutoff.dvir <- calculate_cutoff(df.dvir.blastn$E_val)
e_val_cutoff.dgri <- calculate_cutoff(df.dgri.blastn$E_val)

df.dmoj.blastn <- df.dmoj.blastn[df.dmoj.blastn$E_val < e_val_cutoff.dmoj]
df.dvir.blastn <- df.dvir.blastn[df.dvir.blastn$E_val < e_val_cutoff.dvir]
df.dgri.blastn <- df.dgri.blastn[df.dgri.blastn$E_val < e_val_cutoff.dgri]

df.dmoj.blastn <- transfer.bl.gc(df.dmoj.blastn)
df.dvir.blastn <- transfer.bl.gc(df.dvir.blastn)
df.dgri.blastn <- transfer.bl.gc(df.dgri.blastn)

colnames(df.dmoj.blastn) <- c("Gene_ID","S", 'Scaffold_ID', 'Alignment_Length', 'Begin', 
                              'End', 'E_val', 'Bit_score')
colnames(df.dvir.blastn) <- c("Gene_ID","S", 'Scaffold_ID', 'Alignment_Length', 'Begin', 
                              'End', 'E_val', 'Bit_score')
colnames(df.dgri.blastn) <- c("Gene_ID","S", 'Scaffold_ID', 'Alignment_Length', 'Begin', 
                              'End', 'E_val', 'Bit_score')

write.table(df.dmoj.blastn, 'blastn_dmoj_converted.txt', col.names=T,row.names=F,quote=F, 
            sep='\t')
write.table(df.dvir.blastn, 'blastn_dvir_converted.txt', col.names=T,row.names=F,quote=F, 
            sep='\t')
write.table(df.dgri.blastn, 'blastn_dgri_converted.txt', col.names=T,row.names=F,quote=F, 
            sep='\t')

#calculate_cutoff(df$E_val)
#dmoj_blastn_df<-df[df$E_val<2/140000000,]

# colnames(df) <- c('unique_ID', 'Scaffold_ID', 'Percent_Identity', 
#                   'Alignment_Length', 'Num_Mismathces', 'Num_Gap', 
#                   'mRNA_Begin', 'mRNA_End', 'Begin', 
#                   'End', 'E_val', 'Bit_score',"Gene_ID","S")
