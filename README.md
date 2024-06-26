# Title of this project:
Integration of single cell sequencing data and GWAS summary statistics identifies ORMDL3-mediated cholangiocytes associated with primary biliary cholangitis

This project has been done and published in [Xiang et al. Journal of Nanobiotechnology, 2021](https://jnanobiotechnology.biomedcentral.com/articles/10.1186/s12951-021-01154-2).

# Importance: 
Primary biliary cholangitis (PBC) is a classical autoimmune disease, which is highly influenced by genetic determinants. Many genome-wide association studies (GWAS) have reported that numerous genetic loci were significantly associated with PBC susceptibility. However, the effects of genetic determinants on liver cells for PBC remain largely unknown.

 # Objective: 
 To identify genetics-modulated functional liver cell subsets involved in the pathogenesis of PBC. 

 ![Figure 1](https://github.com/mayunlong89/PBC_project/blob/main/figures/Figure%201.jpg)

# Introduction
Primary biliary cholangitis (PBC), which is formally known as primary biliary cirrhosis until 2016, is a rare chronic cholestatic liver disease characterized by progressive autoimmune-mediated destruction of the small intrahepatic biliary epithelial cells. PBC patients suffering from chronic cholestasis can eventually lead to cirrhosis and hepatic failure without effective treatments. Although ursodeoxycholic acid has been used as the first-line therapeutic agent for PBC, there exist 10% to 20% of PBC patients resistant to ursodeoxycholic acid and developing to advanced-stage liver disease. Previous studies have reported that a combination of genetic and environmental risk factors have an important influence on the aetiology of PBC. Hence, understanding the genetic mechanisms of PBC is becoming a great interest, which may promote the development of individualized therapeutic strategy for PBC. 

With the advance of single cell sequencing techniques, researchers have an effective avenue to discover more refined and novel cell populations for complex diseases. An accruing and large number of single cell RNA sequencing (scRNA-seq) studies on autoimmune diseases, including rheumatoid arthritis, inflammatory bowel disease, systemic lupus erythematosus, and Crohn’s disease, have been reported to parse the heterogeneity of cellular subpopulations at unprecedented resolution. In view of no scRNA-seq study was conducted for uncovering human liver cell types implicated in PBC, we constructed a computational framework to identify risk genes whose genetically expressions associated with PBC and pinpoint cell subpopulations implicated in the etiology of PBC.

# Methods
## GWAS summary statistics on primary biliary cholangitis
   In the current investigation, we employed a genome-wide association study (GWAS) summary data on PBC 10 for uncovering novel PBC-associated risk genes and variants. For this GWAS dataset, which is downloaded from the IEU open GWAS project (https://gwas.mrcieu.ac.uk/), there were 2,764 PBC patients and 10,475 healthy controls used for performing meta-analysis of genome-wide association signals based on European ancestry. Human subjects enrolled in this GWAS were approved by the University Health Network Research Ethics Board, the Mayo Clinic Institutional Review Board, Etico Indipendente IRCCS Istituto Clinico Humanitas, UC Davis Institutional Review Board and the Oxford Research Ethics Committee 10. A standard quality control (QC) pipeline was applied to remove low-quality SNPs from further analysis. The software package of MaCH 27 with the reference of HapMap3 CEU + TSI samples was implemented to perform a genome-wide imputation analysis. There were 1,124,241 SNPs with minor allele frequency > 0.005 and imputation quality score R2 > 0.5 in GWAS samples included in follow-up analyses. 

## Bulk-based expression profiles of primary biliary cholangitis-related risk genes
  To validate the functionality of these identified PBC-associated risk genes, we downloaded two bulk-based expression datasets based on liver tissue (accession number: GSE159676) and blood (accession number: GSE119600) from the Gene Expression Omnibus (GEO) database. For the liver dataset of GSE159676, there were six healthy controls and three PBC cases based on fresh frozen tissue obtained from explanted livers or diagnostic liver biopsies. The Affymetrix Human Gene 1.0 st array was leveraged to obtain bulk expression profiles of 17,046 probes. With respect to the blood dataset of GSE119600 44, bulk-based expression profiles were performed using RNA isolated from whole blood samples from 47 healthy controls and 90 PBC patients. The Illumina HumanHT-12 V4.0 expression beadchip was leveraged to obtain blood transcriptome of 47,230 probes. Differential gene expression (DGE) analyses between controls and patients of both datasets were examined by using Student’s T-test, and P value less than or equal to 0.05 was considered to be of significance. 

## Single-cell expression profiles of primary biliary cholangitis-related risk genes

  To explore the functions of PBC-related risk genes in single cells of liver tissue, we downloaded two independent single cell RNA profiles (accession number: GSE93170 and GSE115469) from the GEO database. With regard to the dataset of GSE93170 45, there were clinically and pathologically diagnosed six healthy controls and six PBC patients enrolled with written informed consent. Peripheral CD4+T cells were used to extract total RNA. The Agilent microarray of SurePrint G3 human GE 8×60K microarray kit was leveraged to produce gene expression profiles according to manufacturer’s protocols. The Student’s T-test was used for assessing the difference between PBC and control group. Furthermore, we also carried out a co-expression pattern analysis of these identified risk genes among PBC group and control group to evaluate whether the co-expression link changed due to the disease status. 
  
  As for the GSE115469 dataset 46, five samples from primary liver patients were used for single cell RNA sequencing based on the 10× Genomics Chromium Single Cell Kits. A total of 8,444 parenchymal and non-parenchymal cells have obtained the transcriptional profiles based on the CellRanger analysis pipeline. The raw digital matrix of gene expression (namely UMI counts per gene per cell) was filtered, normalized and clustered. Cell was omitted if it has a very high (>0.5) mitochondrial genome transcript ratio or a very small library size (<1500). Based on these single-cell profiles using the standard Seurat package 47, there were 20 discrete clustered. We used well-known marker genes assigned these clusters into 13 distinct cell subpopulations, including portal endothelial cells, cholangiocytes, non-inflammatory macrophages, T cells, γδT cells, inflammatory monocytes/macrophages, natural killer (NK)-like cells, red blood cells (RBCs), sinusoidal endothelial cells, mature B cells, stellate cells, plasma cells, and hepatocytes. 

# References
1.	Beuers U, Gershwin ME, Gish RG, et al. Changing nomenclature for PBC: From 'cirrhosis' to 'cholangitis'. J Hepatol. 2015;63(5):1285-1287.
2.	Lindor KD, Bowlus CL, Boyer J, Levy C, Mayo M. Primary Biliary Cholangitis: 2018 Practice Guidance from the American Association for the Study of Liver Diseases. Hepatology. 2019;69(1):394-419.
3.	Rong G, Zhong R, Lleo A, et al. Epithelial cell specificity and apotope recognition by serum autoantibodies in primary biliary cirrhosis. Hepatology. 2011;54(1):196-203.
4.	Tanaka A, Leung PSC, Gershwin ME. Evolution of our understanding of PBC. Best Pract Res Clin Gastroenterol. 2018;34-35:3-9.
5.	Lleo A, Marzorati S, Anaya JM, Gershwin ME. Primary biliary cholangitis: a comprehensive overview. Hepatol Int. 2017;11(6):485-499.
6.	Kaplan MM, Gershwin ME. Primary biliary cirrhosis. N Engl J Med. 2005;353(12):1261-1273.
7.	Selmi C, Mayo MJ, Bach N, et al. Primary biliary cirrhosis in monozygotic and dizygotic twins: genetics, epigenetics, and environment. Gastroenterology. 2004;127(2):485-492.
8.	Shin S, Moh IH, Woo YS, et al. Evidence from a familial case suggests maternal inheritance of primary biliary cholangitis. World J Gastroenterol. 2017;23(39):7191-7197.
9.	Tanaka A, Leung PSC, Gershwin ME. The genetics of primary biliary cholangitis. Curr Opin Gastroenterol. 2019;35(2):93-98.
10.	Cordell HJ, Han Y, Mells GF, et al. International genome-wide meta-analysis identifies new primary biliary cirrhosis risk loci and targetable pathogenic pathways. Nat Commun. 2015;6:8019.
11.	Hirschfield GM, Liu X, Xu C, et al. Primary biliary cirrhosis associated with HLA, IL12A, and IL12RB2 variants. N Engl J Med. 2009;360(24):2544-2555.
12.	Kawashima M, Hitomi Y, Aiba Y, et al. Genome-wide association studies identify PRKCB as a novel genetic susceptibility locus for primary biliary cholangitis in the Japanese population. Hum Mol Genet. 2017;26(3):650-659.
13.	Qiu F, Tang R, Zuo X, et al. A genome-wide association study identifies six novel risk loci for primary biliary cholangitis. Nat Commun. 2017;8:14828.
14.	Liu X, Invernizzi P, Lu Y, et al. Genome-wide meta-analyses identify three loci associated with primary biliary cirrhosis. Nat Genet. 2010;42(8):658-660.
15.	Mells GF, Floyd JA, Morley KI, et al. Genome-wide association study identifies 12 new susceptibility loci for primary biliary cirrhosis. Nat Genet. 2011;43(4):329-332.
16.	Juran BD, Hirschfield GM, Invernizzi P, et al. Immunochip analyses identify a novel risk locus for primary biliary cirrhosis at 13q14, multiple independent associations at four established risk loci and epistasis between 1p31 and 7q32 risk variants. Hum Mol Genet. 2012;21(23):5209-5221.
17.	Ma Y, Huang Y, Zhao S, et al. Integrative Genomics Analysis Reveals a 21q22.11 Locus Contributing Risk to COVID-19. Hum Mol Genet. 2021.
18.	Zhu Z, Zhang F, Hu H, et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. Nat Genet. 2016;48(5):481-487.
19.	Hindorff LA, Sethupathy P, Junkins HA, et al. Potential etiologic and functional implications of genome-wide association loci for human diseases and traits. Proc Natl Acad Sci U S A. 2009;106(23):9362-9367.
20.	Calabrese GM, Mesner LD, Stains JP, et al. Integrating GWAS and Co-expression Network Data Identifies Bone Mineral Density Genes SPTBN1 and MARK3 and an Osteoblast Functional Module. Cell Syst. 2017;4(1):46-59.e44.

