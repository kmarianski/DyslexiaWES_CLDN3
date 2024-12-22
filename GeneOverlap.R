library(GeneOverlap)

dyslexia_genes <- read.delim("/512_genes.txt", header = FALSE, stringsAsFactors = FALSE)
andMe23 = c('STAU1', 'CADM2', 'MSL2', 'CSE1L', 'MRM1', 'CRAT', 'PLCL1', 'LSAMP', 'RFTN2', 'NOL4L-DT', 'CDK12', 'PPP2R3A', 
           'ARFGEF2', 'STAG1', 'MED1', 'FBXL20', 'DLAT', 'C1orf87', 'AC011997.1', 'ASXL1', 'PET112', 'CRYAB', 'HBP1', 
           'KIF3B', 'TM9SF4', 'BOLL', 'PGAP3', 'MOB4', 'DIXDC1', 'NRXN1', 'ESRRG', 'ERBB2', 'HSPE1-MOB4', 'COG5', 
           'AMIGO1', 'ZNHIT3', 'COMMD7', 'MITF', 'SYPL2', 'CHST9', 'DHRS11', 'PBXIP1', 'MAPRE1', 'PCCB', 'DOLPP1', 
           'PPP2R4', 'IER5L', 'GGNBP2', 'C20orf112', 'DNMT3B', 'SORT1', 'PSMA5', 'POU6F2', 'PCGF6', 'EFCAB8', 'CCDC171', 
           'USMG5', 'ALG9 (paralogue)', 'ALG9', 'PRKAR2B', 'CYB561D1', 'GPC6', 'PIH1D2', 'ARL14EP', 'HSPD1', 'TAF5', 
           'AUTS2', 'ATXN7L2', 'BRE', 'SF3B1', 'HSPB2', 'MPPED2', 'INA', 'RBFOX1', 'NPM1', 'MYBPHL', 'HSPB2-C11orf52', 
           'GMFG', 'AL136218.1', 'COQ10B', 'DDX27', 'CUX2', 'SORCS3', 'SETDB2', 'R3HCC1L', 'MYO19', 'TANC2', 'SRPK2', 
           'UNC5D', 'CDHR4', 'SGCD', 'PHF11', 'SEMA3F', 'TMEM131', 'SIK2', 'CRTAC1', 'PPP1R1B', 'PANX2', 'LOXL4', 
           'NKAPD1', 'AMT', 'PDE1C', 'MST1', 'FSHB', 'TMEM182', 'FAM160A1', 'MRAS', 'TCAP', 'CDH18', 'CILP2', 'PIGW', 
           'B3GAT1', 'HEATR5B', 'HSPE1', 'PLAGL2', 'TCTA', 'GPR22', 'SKOR2', 'SH2B3', 'PHF2', 'NDUFA13', 'CREBRF', 
           'TBC1D5', 'CCSER1', 'ARPC5L', 'SAMD4B', 'CTC-260F20.3', 'RNF123', 'SDHD', 'YJEFN3', 'TRAIP', 'SCAF1', 
           'GOLGA1', 'SNX29', 'PPP2R1B', 'LRRD1', 'C11orf52', 'FAM120A', 'TUSC3', 'AKAP9', 'TSSK6', 'SELO', 'CALN1', 
           'PTPN21', 'FDXACB1', 'IRF3', 'WDR38', 'ATXN2', 'MTSS1L', 'POF1B', 'PLCH2', 'CYP51A1', 'ACSL6', 'MFSD9', 
           'FNIP1', 'SCN5A', 'ZNF711', 'ZPBP2', 'GLIS3', 'PRMT1', 'MAP1LC3B', 'SPATA7', 'ZNFX1', 'BNIP1', 'CTC-432M15.3', 
           'WASF3', 'MSANTD4', 'PRKCE', 'SLC35G2', 'SNF8', 'BCL2L12', 'DZIP1', 'RPS6KB1')

# If the file is a simple list, the genes will be in the first column
# Convert the data frame to a vector of gene names
dyslexia <- dyslexia_genes[[1]]

# Define the universe of genes (e.g., all genes in the genome)
total_genes <- 20000  # Example: total number of genes in the genome

# Create a GeneOverlap object
go_obj <- newGeneOverlap(andMe23, dyslexia, genome.size = total_genes)

# Test the significance of the overlap
go_obj <- testGeneOverlap(go_obj)

# Print the results
print(go_obj)
