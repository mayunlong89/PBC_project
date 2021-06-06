# Title of this project:
Integrated GWAS and scRNA-seq analysis identifies genetics-influenced novel genes and cholangiocytes associated with primary biliary cholangitis susceptibility

# Abstract
Importance: Primary biliary cholangitis (PBC) is a classical autoimmune disease, which is highly influenced by genetic factors. Previous genome-wide association studies (GWAS) have reported a group of genetic loci associated with PBC susceptibility. However, the functional effects of genetic loci on liver cell subpopulations for PBC remain largely unknown.
Objective: To highlight novel risk genes and examine whether these genetic risk genes affect liver cell types for PBC based on integrative genomics analysis.
Design, Setting, and Participants: In this integrative genomics analysis, 13,239 European participants were collected from IEU open GWAS project on PBC. There were 1,124,241 qualified SNPs used for GWAS analysis. Expression quantitative trait loci (eQTL) data across 49 tissues were downloaded from the GTEx database. Two bulk-based and two single cell RNA expression profiles were downloaded from the GEO database. Data collection and analysis were performed from August 2020 to June 2021. 
Main outcomes and measures: Summary statistics from the GWAS study were leveraged to extract genetic association signals. We conducted systematic bioinformatics analyses, such as gene-level genetic association analysis, S-MultiXcan integrative genomics analysis, gene-property analysis, drug-gene interaction, gene-gene interaction analysis, and Rolypoly-based scRNA-seq analysis. 
	Results: Based on comprehensive genomics analysis, we found 29 risk genes were significantly associated with PBC. There were 19 genes, including GSNK2B (P = 7.17×10-19), LY6G5B (P = 2.43×10-16), DDAH2 (P = 1.39×10-15), and C6orf48 (P = 9.57×10-9), identified to be novel risk genes for PBC. Based on disease-based enrichment analysis, these genes were significantly overrepresented in autoimmune disease (P = 8.47×10-8), type I diabetes mellitus (P = 2.04×10-4), and immune system diseases (P = 1.41×10-3). Differential gene expression analysis based on three independent RNA data uncovered 22 of 29 genes (75.86%) were significantly expressed among PBC patients compared with controls. Drug-gene interaction analysis revealed 20 of 29 genes (68.96%) showed enrichments in ten potential druggable gene categories. Based on integrative genomics analysis of GWAS with scRNA-seq data, we found that genetics-influenced changiocytes were significantly associated with PBC, and risk genes of ORMDL3, CTSH, TSFM, and ZCRB1 showed highly specific expression in changiocytes. 
	Conclusions and relevance: In this study, we prioritized 29 genes including 19 novel genes were significantly associated with PBC susceptibility, and found genetics-influenced changiocytes had important roles in the etiology of PBC.

# Introduction
  Primary biliary cholangitis (PBC), which is formally known as primary biliary cirrhosis until 2016 1, is a rare chronic cholestatic liver disease characterized by progressive autoimmune-mediated destruction of the small intrahepatic biliary epithelial cells 2,3. PBC patients suffering from chronic cholestasis can eventually lead to cirrhosis and hepatic failure without effective treatments 2,4,5. Although ursodeoxycholic acid has been used as the first-line therapeutic agent for PBC, there exist 10% to 20% of PBC patients resistant to ursodeoxycholic acid and develop to advanced-stage liver disease 2,6. The aetiology of PBC has been extensively reported to be influenced by a combination of genetic and environmental risk factors 7-9. Therefore, there is a great interest in uncovering the genetic architecture that involved in the pathology of PBC, which may promote the development of effective therapeutics for PBC. 

  With the advance of single cell sequencing techniques including single cell RNA sequencing (scRNA-seq), single cell DNA sequencing (scDNA-seq), and mass cytometry, researchers enable to discover more refined and novel cell types for complex diseases 22. An accruing and large number of studies on autoimmune diseases, including rheumatoid arthritis 23, inflammatory bowel disease 24, systemic lupus erythematosus 25, and Crohn’s disease 26, have been reported to parse heterogeneity of cellular subpopulations at unprecedented resolution. In view of there was no single-cell RNA-seq study performed for uncovering liver cell types associated with PBC, we here performed a comprehensive genomics analysis by integrating multi-omics data, including GWAS summary data, eQTL data, bulk gene expression data, and scRNAS-seq data, to identify novel PBC-associated risk genes and genetically pinpoint liver cell subpopulations implicated in the pathogenesis of PBC.

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
21.	Pavlides JM, Zhu Z, Gratten J, McRae AF, Wray NR, Yang J. Predicting gene targets from integrative analyses of summary data from GWAS and eQTL studies for 28 human complex traits. Genome Med. 2016;8(1):84.
22.	Zhang Y, Ma Y, Huang Y, et al. Benchmarking algorithms for pathway activity transformation of single-cell RNA-seq data. Comput Struct Biotechnol J. 2020;18:2953-2961.
23.	Zhang F, Wei K, Slowikowski K, et al. Defining inflammatory cell states in rheumatoid arthritis joint synovial tissues by integrating single-cell transcriptomics and mass cytometry. Nat Immunol. 2019;20(7):928-942.
24.	Corridoni D, Chapman T, Antanaviciute A, Satsangi J, Simmons A. Inflammatory Bowel Disease Through the Lens of Single-cell RNA-seq Technologies. Inflamm Bowel Dis. 2020;26(11):1658-1668.
25.	Nehar-Belaid D, Hong S, Marches R, et al. Mapping systemic lupus erythematosus heterogeneity at the single-cell level. Nat Immunol. 2020;21(9):1094-1106.
26.	Elmentaite R, Ross ADB, Roberts K, et al. Single-Cell Sequencing of Developing Human Gut Reveals Transcriptional Links to Childhood Crohn's Disease. Dev Cell. 2020;55(6):771-783.e775.
27.	Li Y, Willer CJ, Ding J, Scheet P, Abecasis GR. MaCH: using sequence and genotype data to estimate haplotypes and unobserved genotypes. Genet Epidemiol. 2010;34(8):816-834.
28.	de Leeuw CA, Mooij JM, Heskes T, Posthuma D. MAGMA: generalized gene-set analysis of GWAS data. PLoS Comput Biol. 2015;11(4):e1004219.
29.	Xu M, Li J, Xiao Z, Lou J, Pan X, Ma Y. Integrative genomics analysis identifies promising SNPs and genes implicated in tuberculosis risk based on multiple omics datasets. Aging (Albany NY). 2020;12(19):19173-19220.
30.	Auton A, Brooks LD, Durbin RM, et al. A global reference for human genetic variation. Nature. 2015;526(7571):68-74.
31.	Ma X, Wang P, Xu G, Yu F, Ma Y. Integrative genomics analysis of various omics data and networks identify risk genes and variants vulnerable to childhood-onset asthma. BMC Med Genomics. 2020;13(1):123.
32.	Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research. 2000;28(1):27-30.
33.	Wang J, Duncan D, Shi Z, Zhang B. WEB-based GEne SeT AnaLysis Toolkit (WebGestalt): update 2013. Nucleic Acids Res. 2013;41(Web Server issue):W77-83.
34.	Jourquin J, Duncan D, Shi Z, Zhang B. GLAD4U: deriving and prioritizing gene lists from PubMed literature. BMC Genomics. 2012;13 Suppl 8(Suppl 8):S20.
35.	Gene Ontology Consortium: going forward. Nucleic Acids Res. 2015;43(Database issue):D1049-1056.
36.	Ma Y, Li J, Xu Y, et al. Identification of 34 genes conferring genetic and pharmacological risk for the comorbidity of schizophrenia and smoking behaviors. Aging (Albany NY). 2020;12(3):2169-2225.
37.	Zhong Y, Chen L, Li J, et al. Integration of summary data from GWAS and eQTL studies identified novel risk genes for coronary artery disease. Medicine (Baltimore). 2021;100(11):e24769.
38.	Warde-Farley D, Donaldson SL, Comes O, et al. The GeneMANIA prediction server: biological network integration for gene prioritization and predicting gene function. Nucleic Acids Res. 2010;38(Web Server issue):W214-220.
39.	Shannon P, Markiel A, Ozier O, et al. Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome Res. 2003;13(11):2498-2504.
40.	Ma Y, Li MD. Establishment of a Strong Link Between Smoking and Cancer Pathogenesis through DNA Methylation Analysis. Sci Rep. 2017;7(1):1811.
41.	Barbeira AN, Dickinson SP, Bonazzola R, et al. Exploring the phenotypic consequences of tissue specific gene expression variation inferred from GWAS summary statistics. Nat Commun. 2018;9(1):1825.
42.	The Genotype-Tissue Expression (GTEx) project. Nat Genet. 2013;45(6):580-585.
43.	Barbeira AN, Pividori M, Zheng J, Wheeler HE, Nicolae DL, Im HK. Integrating predicted transcriptome from multiple tissues improves association detection. PLoS Genet. 2019;15(1):e1007889.
44.	Ostrowski J, Goryca K, Lazowska I, et al. Common functional alterations identified in blood transcriptome of autoimmune cholestatic liver and inflammatory bowel diseases. Sci Rep. 2019;9(1):7190.
45.	Nakagawa R, Muroyama R, Saeki C, et al. miR-425 regulates inflammatory cytokine production in CD4(+) T cells via N-Ras upregulation in primary biliary cholangitis. J Hepatol. 2017;66(6):1223-1230.
46.	MacParland SA, Liu JC, Ma XZ, et al. Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations. Nat Commun. 2018;9(1):4383.
47.	Macosko EZ, Basu A, Satija R, et al. Highly Parallel Genome-wide Expression Profiling of Individual Cells Using Nanoliter Droplets. Cell. 2015;161(5):1202-1214.
48.	Watanabe K, Umićević Mirkov M, de Leeuw CA, van den Heuvel MP, Posthuma D. Genetic mapping of cell type specificity for complex traits. Nat Commun. 2019;10(1):3222.
49.	Han X, Wang R, Zhou Y, et al. Mapping the Mouse Cell Atlas by Microwell-Seq. Cell. 2018;172(5):1091-1107.e1017.
50.	Etymologia: Bonferroni correction. Emerg Infect Dis. 2015;21(2):289.
51.	Calderon D, Bhaskar A, Knowles DA, et al. Inferring Relevant Cell Types for Complex Traits by Using Single-Cell Gene Expression. Am J Hum Genet. 2017;101(5):686-699.
52.	Hu HJ, Jin EH, Yim SH, et al. Common variants at the promoter region of the APOM confer a risk of rheumatoid arthritis. Exp Mol Med. 2011;43(11):613-621.
53.	Tomer Y, Dolan LM, Kahaly G, et al. Genome wide identification of new genes and pathways in patients with both autoimmune thyroiditis and type 1 diabetes. J Autoimmun. 2015;60:32-39.
54.	von Spee-Mayer C, Siegert E, Abdirama D, et al. Low-dose interleukin-2 selectively corrects regulatory T cell defects in patients with systemic lupus erythematosus. Ann Rheum Dis. 2016;75(7):1407-1415.
55.	Malik FS, Taplin CE. Insulin therapy in children and adolescents with type 1 diabetes. Paediatr Drugs. 2014;16(2):141-150.
56.	Atkinson MA, Eisenbarth GS, Michels AW. Type 1 diabetes. Lancet. 2014;383(9911):69-82.
57.	Sundin DJ, Wolin MJ. Aldesleukin therapy in HIV-infected patients. Am J Health Syst Pharm. 1998;55(14):1520-1523.
58.	Heine A, Held SA, Daecke SN, et al. The JAK-inhibitor ruxolitinib impairs dendritic cell function in vitro and in vivo. Blood. 2013;122(7):1192-1202.
59.	Damsky W, Peterson D, Ramseier J, et al. The emerging role of Janus kinase inhibitors in the treatment of autoimmune and inflammatory diseases. J Allergy Clin Immunol. 2021;147(3):814-826.
60.	Hiasa Y, Akbar SM, Abe M, Michitaka K, Horiike N, Onji M. Dendritic cell subtypes in autoimmune liver diseases; decreased expression of HLA DR and CD123 on type 2 dendritic cells. Hepatol Res. 2002;22(4):241-249.
61.	Harada K, Shimoda S, Ikeda H, et al. Significance of periductal Langerhans cells and biliary epithelial cell-derived macrophage inflammatory protein-3α in the pathogenesis of primary biliary cirrhosis. Liver Int. 2011;31(2):245-253.
62.	Nakamura M, Nishida N, Kawashima M, et al. Genome-wide association study identifies TNFSF15 and POU2AF1 as susceptibility loci for primary biliary cirrhosis in the Japanese population. Am J Hum Genet. 2012;91(4):721-728.
63.	Hitomi Y, Ueno K, Kawai Y, et al. POGLUT1, the putative effector gene driven by rs2293370 in primary biliary cholangitis susceptibility locus chromosome 3q13.33. Sci Rep. 2019;9(1):102.
64.	Tanaka A, Ohira H, Kikuchi K, et al. Genetic association of Fc receptor-like 3 polymorphisms with susceptibility to primary biliary cirrhosis: ethnic comparative study in Japanese and Italian patients. Tissue Antigens. 2011;77(3):239-243.
65.	Zamanou A, Samiotaki M, Panayotou G, Margaritis L, Lymberi P. Fine specificity and subclasses of IgG anti-actin autoantibodies differ in health and disease. J Autoimmun. 2003;20(4):333-344.
66.	Grayson PC, Kaplan MJ. At the Bench: Neutrophil extracellular traps (NETs) highlight novel aspects of innate immune system involvement in autoimmune diseases. J Leukoc Biol. 2016;99(2):253-264.
67.	Mo X, Guo Y, Qian Q, et al. Mendelian randomization analysis revealed potential causal factors for systemic lupus erythematosus. Immunology. 2020;159(3):279-288.
68.	Chung J, Park ES, Kim D, et al. Thyrotropin modulates interferon-gamma-mediated intercellular adhesion molecule-1 gene expression by inhibiting Janus kinase-1 and signal transducer and activator of transcription-1 activation in thyroid cells. Endocrinology. 2000;141(6):2090-2097.
69.	Albertella MR, Jones H, Thomson W, Olavesen MG, Campbell RD. Localization of eight additional genes in the human major histocompatibility complex, including the gene encoding the casein kinase II beta subunit (CSNK2B). Genomics. 1996;36(2):240-251.
70.	Poirier K, Hubert L, Viot G, et al. CSNK2B splice site mutations in patients cause intellectual disability with or without myoclonic epilepsy. Hum Mutat. 2017;38(8):932-941.
71.	Sakaguchi Y, Uehara T, Suzuki H, Kosaki K, Takenouchi T. Truncating mutation in CSNK2B and myoclonic epilepsy. Hum Mutat. 2017;38(11):1611-1612.
72.	Yang CP, Li X, Wu Y, et al. Comprehensive integrative analyses identify GLT8D1 and CSNK2B as schizophrenia risk genes. Nat Commun. 2018;9(1):838.
73.	Niu HM, Yang P, Chen HH, et al. Comprehensive functional annotation of susceptibility SNPs prioritized 10 genes for schizophrenia. Transl Psychiatry. 2019;9(1):56.
74.	Dimitroulas T, Sandoo A, Hodson J, Smith J, Panoulas VF, Kitas GD. Relationship between dimethylarginine dimethylaminohydrolase gene variants and asymmetric dimethylarginine in patients with rheumatoid arthritis. Atherosclerosis. 2014;237(1):38-44.
75.	Fogarty RD, Abhary S, Javadiyan S, et al. Relationship between DDAH gene variants and serum ADMA level in individuals with type 1 diabetes. J Diabetes Complications. 2012;26(3):195-198.
76.	Siegmund T, Donner H, Braun J, Usadel KH, Badenhoop K. HLA-DMA and HLA-DMB alleles in German patients with type 1 diabetes mellitus. Tissue Antigens. 1999;54(3):291-294.
77.	Morel J, Roch-Bras F, Molinari N, Sany J, Eliaou JF, Combe B. HLA-DMA*0103 and HLA-DMB*0104 alleles as novel prognostic factors in rheumatoid arthritis. Ann Rheum Dis. 2004;63(12):1581-1586.
78.	Yen JH, Chen CJ, Tsai WC, Tsai JJ, Ou TT, Liu HW. HLA-DMA and HLA-DMB genotyping in patients with systemic lupus erythematosus. J Rheumatol. 1999;26(9):1930-1933.
79.	Zwart W, Griekspoor A, Kuijl C, et al. Spatial separation of HLA-DM/HLA-DR interactions within MIIC and phagosome-induced immune escape. Immunity. 2005;22(2):221-233.
80.	Hirschfield GM, Xie G, Lu E, et al. Association of primary biliary cirrhosis with variants in the CLEC16A, SOCS1, SPIB and SIAE immunomodulatory genes. Genes Immun. 2012;13(4):328-335.
81.	Swiecki M, Colonna M. The multifaceted biology of plasmacytoid dendritic cells. Nat Rev Immunol. 2015;15(8):471-485.
82.	Lleo A, Maroni L, Glaser S, Alpini G, Marzioni M. Role of cholangiocytes in primary biliary cirrhosis. Semin Liver Dis. 2014;34(3):273-284.
83.	Banales JM, Huebert RC, Karlsen T, Strazzabosco M, LaRusso NF, Gores GJ. Cholangiocyte pathobiology. Nat Rev Gastroenterol Hepatol. 2019;16(5):269-281.
