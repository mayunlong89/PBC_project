

#1-------------------------------------------------------Monocyte counts
monocyte_GWAS <- read.table("monocyte_ieu-b-31.vcf", header=TRUE)

temp <- monocyte_GWAS[,c(1,2,3,4,5,10)]

library(dplyr)
library(tidyr)

# Split name column into firstname and last name
temp1 <- temp[-1,]
temp1 <- temp1 %>% separate(ieu.b.31, c('ES', 'SE','LP','AF','SS','SNP'),sep =":")
colnames(temp1) <- c('CHR','POS','SNP','REF','ALT','ES', 'SE','LP','AF','SS','ID')
monocyte_gwas_sub <- temp1
temp1 <- c()
temp <- c()

monocyte_gwas_sub1 <- monocyte_gwas_sub[,c(1,2,3,4,5,6,7,8)]

#P value has been de-logged
monocyte_gwas_sub1$LP<-10^(-as.numeric(monocyte_gwas_sub1$LP))
colnames(monocyte_gwas_sub1) <- c('CHR','POS','SNP','REF','ALT','ES', 'SE','Pval')

percent_number <- round(0.25*(length(monocyte_gwas_sub1$SNP)))
ID <- sample(monocyte_gwas_sub1$SNP,percent_number)

monocyte_gwas_sub2 <- monocyte_gwas_sub1[which(monocyte_gwas_sub1$SNP %in% ID),]

write.table(monocyte_gwas_sub2,file="monocyte_gwas_subset1_ieu-b-31.vcf")

length(monocyte_gwas_sub$CHR)

###Gene ID annotate to Gene SYMBOL (Gene name)
ref <- read.table("/Users/mayunlong/Desktop/WMU2024/01-Deep learning_代码/01-组学相关代码/NCBI37.3.gene.loc")
colnames(ref) <- c("ID","CHR","Start","End","Strand","SYMBOL")
setwd("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/GWAS-data/MAGMA_result_fuma/five_autoimmune_disease/")
IBD <- read.table("IBD_magma_genes.out",header = TRUE)
#IBD <- IBD[,-which(colnames(IBD)=="gene_name")]
IBD$gene_name <- ref$SYMBOL[match(IBD$GENE,ref$ID)]
write.table(IBD,file="IBD_magma_genes_new.out",col.names = TRUE, row.names = FALSE,quote = FALSE)


temp2 <- temp1[1:100,]

#This one can be used to generate output
write.table(temp2,file="xxxxxx.vcf",row.names = FALSE,col.names = TRUE,quote = FALSE)

#write.table(temp2,file="xxxxxx.vcf",row.names = TRUE,col.names = TRUE,quote = TRUE)


##read file size:
file_info <- file.info("monocyte_gwas_subset_ieu-b-31.vcf")
file_size <-file_info[["size"]]
file_size/1000000000
 
file_info <- file.info("monocyte_gwas_subset1_ieu-b-31.vcf")
file_size <-file_info[["size"]]
file_size/1000000000

