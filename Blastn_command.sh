#--------------------------------------------------code for Dicky------------------------------------------------------------
export PATH="$PATH:/local/data/public/genome_informatics/programs/ncbi-blast-2.5.0+/bin"

target="/local/data/public/ctl43/gia2/db/dmoj-all-chromosome-r1.04.fasta"
#-------------------------------------------------tblastn D.mel against D.moj---------------------------------------------
query="/local/data/public/ctl43/gia2/query/dmel-all-translation-r6.18.fasta"
tblastn -query $query -subject $target -outfmt "6" -out tblastn_moj.txt -evalue 1e-3

#---------------------------------------------------blastn D.mel against D.moj---------------------------------------------
query="/local/data/public/ctl43/gia2/query/dmel-all-transcript-r6.18.fasta"
blastn -query $query -subject $target -task blastn -outfmt "6" -out blastn_moj.txt -evalue 1e-3


#------------------------------ the following clot is for Albert -------------------------------------------------------#
export PATH="$PATH:/local/data/public/genome_informatics_2017/programs/ncbi-blast-2.5.0+/bin"
cd /home/cy302/Genome_Informatics/Assignment_2/alignment
query="/home/cy302/Genome_Informatics/Assignment_2/dmel-all-transcript-r6.18.fasta"
target="/home/cy302/Genome_Informatics/Assignment_2/dvir-all-chromosome-r1.06.fasta"
tblastn -query $query -subject $target -outfmt "6" -out tbalstn_dvir.txt -evalue 1e-3
blastn -query $query -subject $target -outfmt "6" -out balstn_dvir.txt -evalue 1e-3
#-----------------------------------------------------------------------------------------------------------------------#