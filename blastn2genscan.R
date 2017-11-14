setwd('/home/cy302/Genome_Informatics/Assignment_2/')
library('data.table')
df.dmoj.blastn<-fread('blastn_dmoj.txt')
df.dmoj.tblastn <- fread('tblastn_dmoj.txt')
df.dvir.blastn <- fread('blastn_dvir.txt')
df.dvir.tblastn <- fread('tblastn_dvir.txt')
df.dgri.blastn <- fread('blastn_dgri.txt')
df.dgri.blastn <- fread('tblastn_dgri.txt')

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
mat.mrna.gene <- cbind(mrna.list, gene.list)

mrna2gene <- function(l){
  l.alternative <- rep(0, length(l))
  for (i in 1:length(unique(l))){
    s <- which(mrna.list==unique(l)[i])
    # since the mrna list is unqiue
    l.alternative[s] <- gene.list[s]
  }
  
  return(l.alternative)
}

transfer.bl.gc <- function(df){
  query <- mrna2gene(df$`Query/Gene_ID`)
  s <- as.numeric(df$End_Target_pos > df$Start_Target_Pos)
  s[which(s==0)] <- '-'
  s[which(s==1)] <- '+'
  df.2 <- data.frame(cbind(query, df$Start_Target_Pos, 
                           df$End_Target_pos, df$Alignment_Length, df$E_val, df$Bit_score, 
                           df$`Target/Scaffold_ID`))
  names(df.2) <- c('Gene_ID', 'start', 'end', 'length', 'e-val', 'bit_scr', 'scaffold')
  write.table(df.2, 'output.txt',col.names=T,row.names=F,quote=F, sep='\t')
  return(df.2)
}

# i don't actually understand why there are zeros in the output dataframe, i will look into that when i come back tonight, 
# around 9pm, feel free to delete the rows with 0 as their gene IDs and feed the file to your function, you can remove the
# write.table() outside the function to change the name of the output file accordingly, Ta
