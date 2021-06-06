#!/usr/bin/env bash

#MAGMA_gene-based association analysis for primary biliary cholangitis (PBC)
#MAGMA Gene set enrichment analysis for PBC
#File name: PBC_GWAS_UKBiobank_summary_final



#DIRECTORY
export MAGMA_DIR=/share/pub/mayl/Sherlock/MAGMA
export DATA=/share/pub/mayl/Sherlock/XiangbingYu_Primary_biliary_cholangitis/MAGMA_test
export OUTPUT=/share/pub/mayl/Sherlock/XiangbingYu_Primary_biliary_cholangitis/MAGMA_test

#Formating
#cd $DATA

#generating a location file including three Columns: SNP, CHR, POS
#  nohup gawk '{print $1, $7, $8}'  ../PBC_GWAS_UKBiobank_summary_final   > PBC_GWAS_UKBiobank_summary_final.hg19.location &

#generating a --pval file including two Columns: SNP, P
#If MAGMA detects a header in the file it will look for SNP IDs and p-values in the SNP and P column respectively. 
#If no header is found it will use the first column for SNP IDs and the second column for p-values.
#   nohup gawk '{print $1, $4}' ../PBC_GWAS_UKBiobank_summary_final > PBC_GWAS_UKBiobank_summary_final.results_Pval &


#MAGMA annotation:

$MAGMA_DIR/magma \
    --snp-loc  $DATA/PBC_GWAS_UKBiobank_summary_final.hg19.location  \
    --annotate window=20,20 --gene-loc $MAGMA_DIR/NCBI37.3.gene.loc \
    --out $OUTPUT/PBC_GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation  


#gene-based association analysi:
$MAGMA_DIR/magma \
    --bfile $MAGMA_DIR/1000G_data/g1000_eur \
    --pval $DATA/PBC_GWAS_UKBiobank_summary_final.results_Pval \
    N=13239 \
    --gene-annot   $OUTPUT/PBC_GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation.genes.annot  \
    --out $OUTPUT/PBC_GWAS_UKBiobank_summary_final.hg19_SNP_Gene_Analysis_P

#Pathway-based analysis
#MAGMA gene set-based association analysis
$MAGMA_DIR/magma \
    --gene-results $OUTPUT/PBC_GWAS_UKBiobank_summary_final.hg19_SNP_Gene_Analysis_P.genes.raw \
    --set-annot $MAGMA_DIR/KEGG_pathway_novel_2021_for_MAGMA.txt \
    --out $OUTPUT/PBC_GWAS_UKBiobank_summary_final.hg19_for_MAGMA_KEGG_2021_set_results  

