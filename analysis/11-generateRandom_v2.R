
##R version 4.0.4

geneUniverse <- readRDS("outputTabs/geneFiltGenomicRanges_Ensembl_v82_gff3.RDS")
ddd_size <- dim(subset(geneUniverse, ddd==T))[1]
cosmic_size <- dim(subset(geneUniverse, cosmic==T))[1]
dddDominantMOI_size <- dim(subset(geneUniverse, ddd_dominantMOI==T))[1]

#####################
### Null DDD Size ###
#####################

### Sample any gene from the gene universe (also DDD) as gene-set member
set.seed(123)
seeds <- sample(1:1000)
nullDDD <- sapply(seeds, function(x){
  set.seed(x)
  sample(geneUniverse$symbol, ddd_size)
})
colnames(nullDDD) <- paste0("RandomSeed_",as.character(seeds))
nullDDD[1:10,1:5]
# RandomSeed_415 RandomSeed_463 RandomSeed_179 RandomSeed_526 RandomSeed_195
# [1,] "USF1"         "PCDHB6"       "TMEM235"      "LPIN3"        "VIMP"        
# [2,] "LYL1"         "NPM3"         "RPL15"        "CLCN6"        "TMEM240"     
# [3,] "ACAD9"        "PEX14"        "LDOC1L"       "LPAR3"        "DLL4"        
# [4,] "SLC22A18"     "ENDOV"        "ATOH8"        "MCM9"         "SPRN"        
# [5,] "GABRE"        "C19orf47"     "SLIRP"        "SF1"          "PMFBP1"      
# [6,] "KRTDAP"       "TNNC1"        "CXorf56"      "IGFBP2"       "HIST1H3B"    
# [7,] "LRFN5"        "TRIM42"       "KDM6B"        "ECSIT"        "PLSCR2"      
# [8,] "MFSD1"        "ABCB7"        "GCSAML"       "KRBA2"        "FBXL17"      
# [9,] "ZNF470"       "C9orf89"      "BPIFB4"       "EFCAB7"       "PPP1R1C"     
# [10,] "HPS3"         "CLDN16"       "FNDC1"        "LURAP1L"      "COPZ2" 

saveRDS(nullDDD, file="outputTabs/double_downsampling/nullDDD_1000perm.RDS")


##################################
##################################
##################################

###########################
### Null Cosmic-T1 Size ###
###########################

set.seed(124)
seeds <- sample(1:1000)

### Sample any gene from the gene universe (also Cosmic-T1) as gene-set member
nullCosmic <- sapply(seeds, function(x){
  set.seed(x)
  sample(geneUniverse$symbol, cosmic_size)
})
colnames(nullCosmic) <- paste0("RandomSeed_",as.character(seeds))
nullCosmic[1:10,1:5]
# RandomSeed_321 RandomSeed_167 RandomSeed_411 RandomSeed_261 RandomSeed_728
# [1,] "HIST1H2AA"    "C2orf16"      "CAST"         "SLC8A1"       "OR2J2"       
# [2,] "CCDC129"      "TMEM187"      "GMEB1"        "DTX1"         "SLC5A9"      
# [3,] "CYBB"         "SUSD1"        "KIAA0232"     "PDE6A"        "ADAMTS15"    
# [4,] "UBFD1"        "ZMAT5"        "PDE1C"        "TRIM61"       "SUOX"        
# [5,] "DEFB105A"     "C8orf76"      "SS18L1"       "MRPL47"       "ABCF2"       
# [6,] "TUB"          "CORO1A"       "TRIM5"        "LIN52"        "RHD"         
# [7,] "MTA1"         "LIN7C"        "BTG1"         "FAM84B"       "OR1B1"       
# [8,] "GLT8D1"       "DIMT1"        "CPNE9"        "HOXC9"        "DOK5"        
# [9,] "ANKRD30B"     "PLBD2"        "CDK5RAP1"     "MAST4"        "SLC16A1"     
# [10,] "UBA1"         "EFCAB3"       "NUP160"       "FAM134C"      "SLC39A2" 

saveRDS(nullCosmic, file="outputTabs/nullCosmic_1000perm.RDS")




