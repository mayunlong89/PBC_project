#Rscript
#Author: Yunlong Ma
#Date: 2021-03-20
#E-mail: glb-biotech@zju.edu.cn
#100,000 times in silico permutation analysis for risk genes associated with PBC **S-MultiXcan vs. MAGMA vs. S-PrediXcan


#Set the work directory
#setwd("C:\\Users\\Administrator\\Desktop\\06-Simulation_analysis")
setwd("F:\\Desktop\\permutation_analysis")
#setwd("F:\\Desktop\\permutation_analysis\\permut_p_0.05")
set.seed(12345)

#Part I Read data on significant genes and background genes

#Read significant genes of Geneset #1
Sig_1 <- read.table("smultixcan_sig.txt", header=T)
Sig_smultixcan <- Sig_1$Gene_name

#Read background genes of S-MultiXcan-based analysis
Backgroud_1_2<- read.table("smultixcan_background.txt", header=T)
Backgroud_smultixcan_background <- Backgroud_1_2$Gene_name

#Read significant genes from Geneset #2
Sig_2 <- read.table("magma_sig.txt", header=T)
Sig_magma_sig <- Sig_2$Gene_name

#Read background genes of MAGMA-based association analysis 
Backgroud_2<- read.table("magma_background.txt", header=T)
Backgroud_magma_background <- Backgroud_2$Gene_name

#Read significant genes from Geneset #3
Sig_3 <- read.table("liver_sig.txt", header=T)
Sig_liver_sig <- Sig_3$Gene_name

#Read background genes of S-PrediXcan analysis of GWAS summary dataset on PBC with liver eQTL data
Background_3<- read.table("liver_background.txt", header=T)
Background_liver_background <- Background_3$Gene_name

#Read significant genes from Geneset #4
Sig_4 <- read.table("blood_sig.txt", header=T)
Sig_blood_sig <- Sig_4$Gene_name

#Read background genes of S-PrediXcan analysis of GWAS summary dataset on PBC with blood eQTL data
Background_4<- read.table("blood_background.txt", header=T)
Background_blood_background <- Background_4$Gene_name



#Calculate the numebr of genes in each gene set
len_Sig_smultixcan <- length(Sig_smultixcan)
len_Backgroud_smultixcan_background <- length(Backgroud_smultixcan_background)
len_Sig_magma_sig <- length(Sig_magma_sig)
len_Backgroud_magma_background <- length(Backgroud_magma_background)
len_Sig_liver_sig <- length(Sig_liver_sig)
len_Background_liver_background <- length(Background_liver_background)
len_Sig_blood_sig <- length(Sig_blood_sig)
len_Background_blood_background <- length(Background_blood_background)


#Part II establish a function for permutation analysis

#Permutation Function
Permut_analysis <- function(x,y,z){
       
   random_selected_genes <- sample(x,y)
   
   temp <- match(random_selected_genes,z)
   
   random_overlaped_genes <- na.omit(temp)
    
    num<-length(random_overlaped_genes)
   
  return(num)
  
}


#100000 times permutation analysis for S-MultiXcan VS. MAGMA
results_magma_vs_smultixcan <- replicate(100000,Permut_analysis(Backgroud_magma_background,len_Sig_magma_sig,Sig_smultixcan))

#100000 times permutation analysis for S-MultiXcan VS. S-PrediXcan on liver
results_liver_vs_smultixcan <- replicate(100000,Permut_analysis(Background_liver_background,len_Sig_liver_sig,Sig_smultixcan))

#100000 times permutation analysis for S-MultiXcan VS. S-PrediXcan on blood
results_blood_vs_smultixcan <- replicate(100000,Permut_analysis(Background_blood_background,len_Sig_blood_sig,Sig_smultixcan))




#Ploting function
Fig_random <- function(x,y,z){
  
  hist(x, col="red",xlab="Counts of overlapped genes",main=NULL)
  temp1 <- match(y,z)
  Observed_overlaped_genes <- na.omit(temp1)
  Observed_gene_num <- length(Observed_overlaped_genes)
  abline(v=Observed_gene_num,col="darkblue",lty="longdash")
  P_value=length(x[x>Observed_gene_num])/length(x)
  x1= Observed_gene_num
  freq <- table(x)
  y1 = max(freq)
  text(x1,y1,P_value)
  
}


#Visulization for S-MultiXcan VS. MAGMA
Fig_random(results_magma_vs_smultixcan,Sig_magma_sig,Sig_smultixcan)

#Visulization for S-MultiXcan VS. S-PrediXcan on liver
Fig_random(results_liver_vs_smultixcan,Sig_liver_sig,Sig_smultixcan)

#Visulization for S-MultiXcan VS. S-PrediXcan on blood
Fig_random(results_blood_vs_smultixcan,Sig_blood_sig,Sig_smultixcan)
 





#--------------temp---------------------

hist(results_magma_vs_smultixcan, col="#426ab3",xlab="Counts of overlapped genes",main=NULL)
temp1 <- match(Sig_magma_sig,Sig_smultixcan)
Observed_overlaped_genes <- na.omit(temp1)
Observed_gene_num <- length(Observed_overlaped_genes)
abline(v=Observed_gene_num,col="darkblue",lty="longdash")
P_value=length(x[x>Observed_gene_num])/length(x)
x1= Observed_gene_num
freq <- table(x)
y1 = max(freq)
text(x1,y1,P_value)


hist(results_liver_vs_smultixcan, col="#426ab3",xlab="Counts of overlapped genes",main=NULL)
temp1 <- match(Sig_liver_sig,Sig_smultixcan)
Observed_overlaped_genes <- na.omit(temp1)
Observed_gene_num <- length(Observed_overlaped_genes)
abline(v=Observed_gene_num,col="darkblue",lty="longdash")
P_value=length(x[x>Observed_gene_num])/length(x)
x1= Observed_gene_num
freq <- table(x)
y1 = max(freq)
text(x1,y1,P_value)




hist(results_blood_vs_smultixcan, col="#426ab3",xlab="Counts of overlapped genes",main=NULL)
temp1 <- match(Sig_blood_sig,Sig_smultixcan)
Observed_overlaped_genes <- na.omit(temp1)
Observed_gene_num <- length(Observed_overlaped_genes)
abline(v=Observed_gene_num,col="darkblue",lty="longdash")
P_value=length(x[x>Observed_gene_num])/length(x)
x1= Observed_gene_num
freq <- table(x)
y1 = max(freq)
text(x1,y1,P_value)






#End
