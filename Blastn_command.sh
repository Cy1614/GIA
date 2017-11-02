#--------------------------------------------------code for Dicky------------------------------------------------------------
cd /local/data/public/ctl43/gia2/result
export PATH="$PATH:/local/data/public/genome_informatics/programs/ncbi-blast-2.5.0+/bin"
target="/local/data/public/ctl43/gia2/db/dmoj-all-chromosome-r1.04.fasta"

query_tblastn="/local/data/public/ctl43/gia2/query/dmel-all-translation-r6.18.fasta"
query_blastn="/local/data/public/ctl43/gia2/query/dmel-all-transcript-r6.18.fasta"
tblastn -query $query_blastn -subject $target -outfmt "6" -out tblastn_moj.txt -evalue 1e-3
blastn -query $query_blastn -subject $target -task blastn -outfmt "6" -out blastn_moj.txt -evalue 1e-3


#------------------------------ the following clot is for Albert -------------------------------------------------------#
export PATH="$PATH:/local/data/public/genome_informatics_2017/programs/ncbi-blast-2.5.0+/bin"
cd /home/cy302/Genome_Informatics/Assignment_2/alignment
query_c="/home/cy302/Genome_Informatics/Assignment_2/dmel-all-transcript-r6.18.fasta"
query_l="/home/cy302/Genome_Informatics/Assignment_2/dmel-all-translation-r6.18.fasta"
target="/home/cy302/Genome_Informatics/Assignment_2/dvir-all-chromosome-r1.06.fasta"
blastn -query $query_c -subject $target -outfmt "6" -out balstn_dvir.txt -evalue 1e-3
tblastn -query $query_l -subject $target -outfmt "6" -out tbalstn_dvir.txt -evalue 1e-3
#-----------------------------------------------------------------------------------------------------------------------#
