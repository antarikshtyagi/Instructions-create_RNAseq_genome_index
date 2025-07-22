# Instructions-create_RNAseq_genome_index

module load HISAT2
hisat2-build genome.fa genome.fa
##gff3 to GTF
mlm && conda activate ant_base
gffread my.gff3 -T -o my.gtf (e.g. gffread ../Triasper1_GeneCatalog_20210904.gff3 -T -o genes.gtf)

module load GATK
gatk CreateSequenceDictionary -R genome.fa

grep 'rRNA' genes.gtf  > rRNA.gtf

gtfToGenePred -genePredExt rRNA.gtf rRNA.gtf.genePred \
awk 'BEGIN{OFS="\t"} {print $2, $4, $5, $3, $12}' \
rRNA.gtf.genePred > rRNA.gtf.genePred.shortlist

#combine the files
cat genome.dict rRNA.gtf.genePred.shortlist > rRNA.intervals.txt

##Convert gene annotations from GTF to genePred refFlat format:
##without gene names (geneID only)
Conda activate ant_base
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons genes.gtf /dev/stdout | awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' >genes_refFlat.txt

##with gene namesgt	
~/soft/gtfToGenePred     -genePredExt     -geneNameAsName2     -ignoreGroupsWithoutExons  genes.gtf   /dev/stdout |     awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' >genes_refFlat.txt
module load miniconda
conda activate ant_base

agat_sp_extract_attributes.pl -gff genome.gff \
-att Parent,Dbxref,gene,product -m \
-p CDS -o salmonella_ASM2216v1_gtfinfo.txt

agat_sp_extract_attributes.pl -gff genome.gff -att locus_tag,gene,Dbxref,product -m -p level2 -o salmonella_ASM2216v1_gtfinfo7.txt

#Make gtfinfo

library(dplyr)
library(tidyr)
library(DOSE)
library(GO.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
mm = read.table("projects/DESeq2/gencode.vM37.basic.annotation.gtf.table", header = T)
ids <- bitr(mm$ENSEMBL, fromType="ENSEMBL", toType=c("ENTREZID", "GENENAME"), OrgDb="org.Mm.eg.db")
gene_table <- full_join(mm, ids, by = "ENSEMBL")
gene_table = gene_table[!duplicated(gene_table$GeneVer), ]
write.table(gene_table, file = "projects/DESeq2/GRCm39.vM37.basic.annotation.gtf.info.txt", sep = "\t", row.names = F, quote = F)
Rename the columns to:
GeneVer    GeneID    GeneType    GeneSymbol    ENTREZID    GENENAME

