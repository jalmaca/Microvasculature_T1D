# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 05/23/2022
# R version 4.2.2 (2022-31-10) 'Innocent and Trusting'
install.packages("plyr", dependencies=T)
# Loading packages
suppressWarnings(
  {
    library(ggplot2)
    library(cowplot)
    library(Matrix)
    library(ggridges)
    library(ggrepel)
    library(plyr)
    library(dplyr)
    library(Seurat)
    library(monocle3)
    library(plotly)
    library(clustree)
    library(patchwork)
    library(future)
    library(DoubletFinder)
    library(EnhancedVolcano)
    }
)

# SESSION INFORMATION AND SETUP ####
sessionInfo()
R.Version()

# Set working directory or save as a project in a specific folder in your PC
setwd(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\HPA Refrence)")
set.seed(123)

# Check for Working Directory
getwd()

# Load HPAP data
HPAP <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\HPA Refrence\Data files\T1D_T2D_20220428.rds)")

# Save metadata
#sample_meta <- table(HPAP@meta.data[["sample_id"]])
#write.csv(sample_meta, file=r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\HPA Refrence\Data files\Metadatafiles\meta_data.csv)")

# Remove doublets from data
table(HPAP@meta.data[["scDblFinder.class"]]) # command shows all are singlets and theres no need to remove doublets because there are none
#Idents(HPAP) <- "scDblFinder.class"
#HPAP_single <- subset(HPAP, idents = "singlet")

# Re-classification Jalab
Idents(HPAP) <- "sample_id"
HPAP$cell_typev2.0 <- plyr::mapvalues(
  x= HPAP$sample_id,
  from = c("HPAP021_75162_T1DIslets",
           "HPAP022_75751_Adult-ControlIslets",
           "HPAP023_75798_T1DIslets",
           "HPAP024_75799_Adult-ControlIslets",
           "HPAP026_77407_Adult-ControlIslets",
           "HPAP027_78239_Adult-ControlIslets",
           "HPAP028_78240_T1DIslets",
           "HPAP029_78242_Adult-ControlIslets",
           "HPAP032_79863_T1DIslets",
           "HPAP034_81423_Child-ControlIslets",
           "HPAP035_81832_Adult-ControlIslets",
           "HPAP036_81833_Adult-ControlIslets",
           "HPAP037_82543_Adult-ControlIslets",
           "HPAP038_84764_Child-ControlIslets",
           "HPAP039_86672_Child-ControlIslets",
           "HPAP040_86673_Adult-ControlIslets",
           "HPAP042_86674_Child-ControlIslets",
           "HPAP043_86675_Child-ControlIslets",
           "HPAP044_86676_Child-ControlIslets",
           "HPAP045_86677_Adult-ControlIslets",
           "HPAP047_97169_Child-ControlIslets",
           "HPAP049_97170_Adult-ControlIslets",
           "HPAP050_97171_Adult-ControlIslets",
           "HPAP051_97173_T2DIslets",
           "HPAP052_98611_Adult-ControlIslets",
           "HPAP053_98612_Adult-ControlIslets",
           "HPAP054_98613_Adult-ControlIslets",
           "HPAP055_98614_T1DIslets",
           "HPAP056_98615_Adult-ControlIslets",
           "HPAP057_98616_T2DIslets",
           "HPAP058_98617_T2DIslets",
           "HPAP059_98618_Adult-ControlIslets",
           "HPAP061_99868_T2DIslets",
           "HPAP062_99869_T2DIslets",
           "HPAP063_99870_Adult-ControlIslets",
           "HPAP064_99871_T1DIslets",
           "HPAP065_99872_T2DIslets",
           "HPAP067_99873_Adult-ControlIslets",
           "HPAP069_99874_Adult-ControlIslets",
           "HPAP070_99875_T2DIslets",
           "HPAP071_103407_T1DIslets",
           "HPAP072_103408_pre-t1DIslets",
           "HPAP074_103409_Adult-ControlIslets",
           "HPAP075_103410_Adult-ControlIslets",
           "HPAP077_103411_Adult-ControlIslets",
           "HPAP079_103412_T2DIslets",
           "HPAP080_103413_Adult-ControlIslets",
           "HPAP081_103414_T2DIslets",
           "HPAP082_103415_Adult-ControlIslets",
           "HPAP083_103416_T2DIslets",
           "HPAP084_105149_T1DIslets",
           "HPAP085_103417_T2DIslets",
           "HPAP087_103418_T1DIslets",
           "HPAP088_105150_T2DIslets",
           "HPAP090_105151_T2DIslets",
           "HPAP091_105152_T2DIslets",
           "HPAP092_105153_Adult-ControlIslets",
           "HPAP093_105154_Adult-ControlIslets",
           "HPAP099_106838_Adult-ControlIslets",
           "HPAP100_106839_T2DIslets",
           "HPAP101_106840_Adult-ControlIslets",
           "HPAP103_106841_Adult-ControlIslets",
           "HPAP104_107485_Child-ControlIslets",
           "HPAP105_107486_Adult-ControlIslets",
           "HPAP106_107487_T2DIslets",
           "HPAP107_109359_pre-t1DIslets",
           "HPAP109_109361_T2DIslets"),
  to = c("t1d_long",      #HPAP021_75162_T1DIslets
         "?", #HPAP022_75751_Adult-ControlIslets
         "t1d_long",      #HPAP023_75798_T1DIslets
         "GADA",     #HPAP024_75799_Adult-ControlIslets
         "nd", #HPAP026_77407_Adult-ControlIslets
         "?", #HPAP027_78239_Adult-ControlIslets
         "?", #HPAP028_78240_T1DIslets
         "GADA", #HPAP029_78242_Adult-ControlIslets
         "t1d", #HPAP032_79863_T1DIslets
         "nd", #HPAP034_81423_Child-ControlIslets
         "nd", #HPAP035_81832_Adult-ControlIslets
         "nd", #HPAP036_81833_Adult-ControlIslets
         "nd", #HPAP037_82543_Adult-ControlIslets
         "GADA", #HPAP038_84764_Child-ControlIslets
         "?", #HPAP039_86672_Child-ControlIslets
         "nd", #HPAP040_86673_Adult-ControlIslets
         "?", #HPAP042_86674_Child-ControlIslets
         "?", #HPAP043_86675_Child-ControlIslets
         "?", #HPAP044_86676_Child-ControlIslets
         "GADA", #HPAP045_86677_Adult-ControlIslets
         "nd", #HPAP047_97169_Child-ControlIslets
         "?", #HPAP049_97170_Adult-ControlIslets
         "GADA", #HPAP050_97171_Adult-ControlIslets
         "?", #HPAP051_97173_T2DIslets
         "?", #HPAP052_98611_Adult-ControlIslets
         "?", #HPAP053_98612_Adult-ControlIslets
         "?", #HPAP054_98613_Adult-ControlIslets
         "t1d_long", #HPAP055_98614_T1DIslets
         "?", #HPAP056_98615_Adult-ControlIslets
         "?", #HPAP057_98616_T2DIslets
         "?", #HPAP058_98617_T2DIslets
         "?", #HPAP059_98618_Adult-ControlIslets
         "?", #HPAP061_99868_T2DIslets
         "?", #HPAP062_99869_T2DIslets
         "?", #HPAP063_99870_Adult-ControlIslets
         "t1d", #HPAP064_99871_T1DIslets
         "?", #HPAP065_99872_T2DIslets
         "?", #HPAP067_99873_Adult-ControlIslets
         "?", #HPAP069_99874_Adult-ControlIslets
         "?", #HPAP070_99875_T2DIslets
         "t1d", #HPAP071_103407_T1DIslets
         "GADA", #HPAP072_103408_pre-t1DIslets
         "?", #HPAP074_103409_Adult-ControlIslets
         "?", #HPAP075_103410_Adult-ControlIslets
         "?", #HPAP077_103411_Adult-ControlIslets
         "?", #HPAP079_103412_T2DIslets
         "?", #HPAP080_103413_Adult-ControlIslets
         "?", #HPAP081_103414_T2DIslets
         "nd", #HPAP082_103415_Adult-ControlIslets
         "?", #HPAP083_103416_T2DIslets
         "t1d", #HPAP084_105149_T1DIslets
         "?", #HPAP085_103417_T2DIslets
         "t1d_long", #HPAP087_103418_T1DIslets
         "?", #HPAP088_105150_T2DIslets
         "?", #HPAP090_105151_T2DIslets
         "?", #HPAP091_105152_T2DIslets
         "GADA", #HPAP092_105153_Adult-ControlIslets
         "?", #HPAP093_105154_Adult-ControlIslets
         "nd", #HPAP099_106838_Adult-ControlIslets
         "?", #HPAP100_106839_T2DIslets
         "?", #HPAP101_106840_Adult-ControlIslets
         "?", #HPAP103_106841_Adult-ControlIslets
         "?", #HPAP104_107485_Child-ControlIslets
         "?", #HPAP105_107486_Adult-ControlIslets
         "?", #HPAP106_107487_T2DIslets
         "?", #HPAP107_109359_pre-t1DIslets
         "?" #HPAP109_109361_T2DIslets
         )
  )

# Re-classification Sex
Idents(HPAP) <- "sample_id"
HPAP$sex <- plyr::mapvalues(
  x= HPAP$sample_id,
  from = c("HPAP021_75162_T1DIslets",
           "HPAP022_75751_Adult-ControlIslets",
           "HPAP023_75798_T1DIslets",
           "HPAP024_75799_Adult-ControlIslets",
           "HPAP026_77407_Adult-ControlIslets",
           "HPAP027_78239_Adult-ControlIslets",
           "HPAP028_78240_T1DIslets",
           "HPAP029_78242_Adult-ControlIslets",
           "HPAP032_79863_T1DIslets",
           "HPAP034_81423_Child-ControlIslets",
           "HPAP035_81832_Adult-ControlIslets",
           "HPAP036_81833_Adult-ControlIslets",
           "HPAP037_82543_Adult-ControlIslets",
           "HPAP038_84764_Child-ControlIslets",
           "HPAP039_86672_Child-ControlIslets",
           "HPAP040_86673_Adult-ControlIslets",
           "HPAP042_86674_Child-ControlIslets",
           "HPAP043_86675_Child-ControlIslets",
           "HPAP044_86676_Child-ControlIslets",
           "HPAP045_86677_Adult-ControlIslets",
           "HPAP047_97169_Child-ControlIslets",
           "HPAP049_97170_Adult-ControlIslets",
           "HPAP050_97171_Adult-ControlIslets",
           "HPAP051_97173_T2DIslets",
           "HPAP052_98611_Adult-ControlIslets",
           "HPAP053_98612_Adult-ControlIslets",
           "HPAP054_98613_Adult-ControlIslets",
           "HPAP055_98614_T1DIslets",
           "HPAP056_98615_Adult-ControlIslets",
           "HPAP057_98616_T2DIslets",
           "HPAP058_98617_T2DIslets",
           "HPAP059_98618_Adult-ControlIslets",
           "HPAP061_99868_T2DIslets",
           "HPAP062_99869_T2DIslets",
           "HPAP063_99870_Adult-ControlIslets",
           "HPAP064_99871_T1DIslets",
           "HPAP065_99872_T2DIslets",
           "HPAP067_99873_Adult-ControlIslets",
           "HPAP069_99874_Adult-ControlIslets",
           "HPAP070_99875_T2DIslets",
           "HPAP071_103407_T1DIslets",
           "HPAP072_103408_pre-t1DIslets",
           "HPAP074_103409_Adult-ControlIslets",
           "HPAP075_103410_Adult-ControlIslets",
           "HPAP077_103411_Adult-ControlIslets",
           "HPAP079_103412_T2DIslets",
           "HPAP080_103413_Adult-ControlIslets",
           "HPAP081_103414_T2DIslets",
           "HPAP082_103415_Adult-ControlIslets",
           "HPAP083_103416_T2DIslets",
           "HPAP084_105149_T1DIslets",
           "HPAP085_103417_T2DIslets",
           "HPAP087_103418_T1DIslets",
           "HPAP088_105150_T2DIslets",
           "HPAP090_105151_T2DIslets",
           "HPAP091_105152_T2DIslets",
           "HPAP092_105153_Adult-ControlIslets",
           "HPAP093_105154_Adult-ControlIslets",
           "HPAP099_106838_Adult-ControlIslets",
           "HPAP100_106839_T2DIslets",
           "HPAP101_106840_Adult-ControlIslets",
           "HPAP103_106841_Adult-ControlIslets",
           "HPAP104_107485_Child-ControlIslets",
           "HPAP105_107486_Adult-ControlIslets",
           "HPAP106_107487_T2DIslets",
           "HPAP107_109359_pre-t1DIslets",
           "HPAP109_109361_T2DIslets"),
  to = c("f",      #HPAP021_75162_T1DIslets
         "f", #HPAP022_75751_Adult-ControlIslets
         "f",      #HPAP023_75798_T1DIslets
         "m",     #HPAP024_75799_Adult-ControlIslets
         "m", #HPAP026_77407_Adult-ControlIslets
         "f", #HPAP027_78239_Adult-ControlIslets
         "m", #HPAP028_78240_T1DIslets
         "m", #HPAP029_78242_Adult-ControlIslets
         "f", #HPAP032_79863_T1DIslets
         "m", #HPAP034_81423_Child-ControlIslets
         "m", #HPAP035_81832_Adult-ControlIslets
         "f", #HPAP036_81833_Adult-ControlIslets
         "f", #HPAP037_82543_Adult-ControlIslets
         "m", #HPAP038_84764_Child-ControlIslets
         "f", #HPAP039_86672_Child-ControlIslets
         "m", #HPAP040_86673_Adult-ControlIslets
         "m", #HPAP042_86674_Child-ControlIslets
         "m", #HPAP043_86675_Child-ControlIslets
         "f", #HPAP044_86676_Child-ControlIslets
         "f", #HPAP045_86677_Adult-ControlIslets
         "m", #HPAP047_97169_Child-ControlIslets
         "m", #HPAP049_97170_Adult-ControlIslets
         "f", #HPAP050_97171_Adult-ControlIslets
         "f", #HPAP051_97173_T2DIslets
         "m", #HPAP052_98611_Adult-ControlIslets
         "f", #HPAP053_98612_Adult-ControlIslets
         "f", #HPAP054_98613_Adult-ControlIslets
         "m", #HPAP055_98614_T1DIslets
         "m", #HPAP056_98615_Adult-ControlIslets
         "f", #HPAP057_98616_T2DIslets
         "f", #HPAP058_98617_T2DIslets
         "m", #HPAP059_98618_Adult-ControlIslets
         "f", #HPAP061_99868_T2DIslets
         "m", #HPAP062_99869_T2DIslets
         "f", #HPAP063_99870_Adult-ControlIslets
         "m", #HPAP064_99871_T1DIslets
         "m", #HPAP065_99872_T2DIslets
         "m", #HPAP067_99873_Adult-ControlIslets
         "f", #HPAP069_99874_Adult-ControlIslets
         "m", #HPAP070_99875_T2DIslets
         "f", #HPAP071_103407_T1DIslets
         "m", #HPAP072_103408_pre-t1DIslets
         "f", #HPAP074_103409_Adult-ControlIslets
         "m", #HPAP075_103410_Adult-ControlIslets
         "m", #HPAP077_103411_Adult-ControlIslets
         "f", #HPAP079_103412_T2DIslets
         "m", #HPAP080_103413_Adult-ControlIslets
         "f", #HPAP081_103414_T2DIslets
         "m", #HPAP082_103415_Adult-ControlIslets
         "m", #HPAP083_103416_T2DIslets
         "f", #HPAP084_105149_T1DIslets
         "f", #HPAP085_103417_T2DIslets
         "f", #HPAP087_103418_T1DIslets
         "m", #HPAP088_105150_T2DIslets
         "f", #HPAP090_105151_T2DIslets
         "f", #HPAP091_105152_T2DIslets
         "m", #HPAP092_105153_Adult-ControlIslets
         "m", #HPAP093_105154_Adult-ControlIslets
         "f", #HPAP099_106838_Adult-ControlIslets
         "m", #HPAP100_106839_T2DIslets
         "f", #HPAP101_106840_Adult-ControlIslets
         "f", #HPAP103_106841_Adult-ControlIslets
         "m", #HPAP104_107485_Child-ControlIslets
         "f", #HPAP105_107486_Adult-ControlIslets
         "m", #HPAP106_107487_T2DIslets
         "m", #HPAP107_109359_pre-t1DIslets
         "f" #HPAP109_109361_T2DIslets
  )
)

# Re-classification ancestry
Idents(HPAP) <- "sample_id"
HPAP$ancestry <- plyr::mapvalues(
  x= HPAP$sample_id,
  from = c("HPAP021_75162_T1DIslets",
           "HPAP022_75751_Adult-ControlIslets",
           "HPAP023_75798_T1DIslets",
           "HPAP024_75799_Adult-ControlIslets",
           "HPAP026_77407_Adult-ControlIslets",
           "HPAP027_78239_Adult-ControlIslets",
           "HPAP028_78240_T1DIslets",
           "HPAP029_78242_Adult-ControlIslets",
           "HPAP032_79863_T1DIslets",
           "HPAP034_81423_Child-ControlIslets",
           "HPAP035_81832_Adult-ControlIslets",
           "HPAP036_81833_Adult-ControlIslets",
           "HPAP037_82543_Adult-ControlIslets",
           "HPAP038_84764_Child-ControlIslets",
           "HPAP039_86672_Child-ControlIslets",
           "HPAP040_86673_Adult-ControlIslets",
           "HPAP042_86674_Child-ControlIslets",
           "HPAP043_86675_Child-ControlIslets",
           "HPAP044_86676_Child-ControlIslets",
           "HPAP045_86677_Adult-ControlIslets",
           "HPAP047_97169_Child-ControlIslets",
           "HPAP049_97170_Adult-ControlIslets",
           "HPAP050_97171_Adult-ControlIslets",
           "HPAP051_97173_T2DIslets",
           "HPAP052_98611_Adult-ControlIslets",
           "HPAP053_98612_Adult-ControlIslets",
           "HPAP054_98613_Adult-ControlIslets",
           "HPAP055_98614_T1DIslets",
           "HPAP056_98615_Adult-ControlIslets",
           "HPAP057_98616_T2DIslets",
           "HPAP058_98617_T2DIslets",
           "HPAP059_98618_Adult-ControlIslets",
           "HPAP061_99868_T2DIslets",
           "HPAP062_99869_T2DIslets",
           "HPAP063_99870_Adult-ControlIslets",
           "HPAP064_99871_T1DIslets",
           "HPAP065_99872_T2DIslets",
           "HPAP067_99873_Adult-ControlIslets",
           "HPAP069_99874_Adult-ControlIslets",
           "HPAP070_99875_T2DIslets",
           "HPAP071_103407_T1DIslets",
           "HPAP072_103408_pre-t1DIslets",
           "HPAP074_103409_Adult-ControlIslets",
           "HPAP075_103410_Adult-ControlIslets",
           "HPAP077_103411_Adult-ControlIslets",
           "HPAP079_103412_T2DIslets",
           "HPAP080_103413_Adult-ControlIslets",
           "HPAP081_103414_T2DIslets",
           "HPAP082_103415_Adult-ControlIslets",
           "HPAP083_103416_T2DIslets",
           "HPAP084_105149_T1DIslets",
           "HPAP085_103417_T2DIslets",
           "HPAP087_103418_T1DIslets",
           "HPAP088_105150_T2DIslets",
           "HPAP090_105151_T2DIslets",
           "HPAP091_105152_T2DIslets",
           "HPAP092_105153_Adult-ControlIslets",
           "HPAP093_105154_Adult-ControlIslets",
           "HPAP099_106838_Adult-ControlIslets",
           "HPAP100_106839_T2DIslets",
           "HPAP101_106840_Adult-ControlIslets",
           "HPAP103_106841_Adult-ControlIslets",
           "HPAP104_107485_Child-ControlIslets",
           "HPAP105_107486_Adult-ControlIslets",
           "HPAP106_107487_T2DIslets",
           "HPAP107_109359_pre-t1DIslets",
           "HPAP109_109361_T2DIslets"),
  to = c("w",      #HPAP021_75162_T1DIslets
         "w", #HPAP022_75751_Adult-ControlIslets
         "w",      #HPAP023_75798_T1DIslets
         "w",     #HPAP024_75799_Adult-ControlIslets
         "w", #HPAP026_77407_Adult-ControlIslets
         "w", #HPAP027_78239_Adult-ControlIslets
         "w", #HPAP028_78240_T1DIslets
         "w", #HPAP029_78242_Adult-ControlIslets
         "w", #HPAP032_79863_T1DIslets
         "w", #HPAP034_81423_Child-ControlIslets
         "w", #HPAP035_81832_Adult-ControlIslets
         "w", #HPAP036_81833_Adult-ControlIslets
         "w", #HPAP037_82543_Adult-ControlIslets
         "w", #HPAP038_84764_Child-ControlIslets
         "w", #HPAP039_86672_Child-ControlIslets
         "w", #HPAP040_86673_Adult-ControlIslets
         "w", #HPAP042_86674_Child-ControlIslets
         "h", #HPAP043_86675_Child-ControlIslets
         "w", #HPAP044_86676_Child-ControlIslets
         "w", #HPAP045_86677_Adult-ControlIslets
         "w", #HPAP047_97169_Child-ControlIslets
         "w", #HPAP049_97170_Adult-ControlIslets
         "h", #HPAP050_97171_Adult-ControlIslets
         "aa", #HPAP051_97173_T2DIslets
         "aa", #HPAP052_98611_Adult-ControlIslets
         "w", #HPAP053_98612_Adult-ControlIslets
         "w", #HPAP054_98613_Adult-ControlIslets
         "h", #HPAP055_98614_T1DIslets
         "w", #HPAP056_98615_Adult-ControlIslets
         "w", #HPAP057_98616_T2DIslets
         "aa", #HPAP058_98617_T2DIslets
         "w", #HPAP059_98618_Adult-ControlIslets
         "aa", #HPAP061_99868_T2DIslets
         "w", #HPAP062_99869_T2DIslets
         "w", #HPAP063_99870_Adult-ControlIslets
         "aa", #HPAP064_99871_T1DIslets
         "aa", #HPAP065_99872_T2DIslets
         "h", #HPAP067_99873_Adult-ControlIslets
         "w", #HPAP069_99874_Adult-ControlIslets
         "aa", #HPAP070_99875_T2DIslets
         "w", #HPAP071_103407_T1DIslets
         "h", #HPAP072_103408_pre-t1DIslets
         "w", #HPAP074_103409_Adult-ControlIslets
         "w", #HPAP075_103410_Adult-ControlIslets
         "w", #HPAP077_103411_Adult-ControlIslets
         "h", #HPAP079_103412_T2DIslets
         "aa", #HPAP080_103413_Adult-ControlIslets
         "w", #HPAP081_103414_T2DIslets
         "w", #HPAP082_103415_Adult-ControlIslets
         "aa", #HPAP083_103416_T2DIslets
         "w", #HPAP084_105149_T1DIslets
         "w", #HPAP085_103417_T2DIslets
         "w", #HPAP087_103418_T1DIslets
         "w", #HPAP088_105150_T2DIslets
         "w", #HPAP090_105151_T2DIslets
         "h", #HPAP091_105152_T2DIslets
         "h", #HPAP092_105153_Adult-ControlIslets
         "w", #HPAP093_105154_Adult-ControlIslets
         "h", #HPAP099_106838_Adult-ControlIslets
         "w", #HPAP100_106839_T2DIslets
         "h", #HPAP101_106840_Adult-ControlIslets
         "w", #HPAP103_106841_Adult-ControlIslets
         "h", #HPAP104_107485_Child-ControlIslets
         "h", #HPAP105_107486_Adult-ControlIslets
         "w", #HPAP106_107487_T2DIslets
         "w", #HPAP107_109359_pre-t1DIslets
         "h" #HPAP109_109361_T2DIslets
  )
)

# Re-classification ancestry
Idents(HPAP) <- "sample_id"
HPAP$disease <- plyr::mapvalues(
  x= HPAP$sample_id,
  from = c("HPAP021_75162_T1DIslets",
           "HPAP022_75751_Adult-ControlIslets",
           "HPAP023_75798_T1DIslets",
           "HPAP024_75799_Adult-ControlIslets",
           "HPAP026_77407_Adult-ControlIslets",
           "HPAP027_78239_Adult-ControlIslets",
           "HPAP028_78240_T1DIslets",
           "HPAP029_78242_Adult-ControlIslets",
           "HPAP032_79863_T1DIslets",
           "HPAP034_81423_Child-ControlIslets",
           "HPAP035_81832_Adult-ControlIslets",
           "HPAP036_81833_Adult-ControlIslets",
           "HPAP037_82543_Adult-ControlIslets",
           "HPAP038_84764_Child-ControlIslets",
           "HPAP039_86672_Child-ControlIslets",
           "HPAP040_86673_Adult-ControlIslets",
           "HPAP042_86674_Child-ControlIslets",
           "HPAP043_86675_Child-ControlIslets",
           "HPAP044_86676_Child-ControlIslets",
           "HPAP045_86677_Adult-ControlIslets",
           "HPAP047_97169_Child-ControlIslets",
           "HPAP049_97170_Adult-ControlIslets",
           "HPAP050_97171_Adult-ControlIslets",
           "HPAP051_97173_T2DIslets",
           "HPAP052_98611_Adult-ControlIslets",
           "HPAP053_98612_Adult-ControlIslets",
           "HPAP054_98613_Adult-ControlIslets",
           "HPAP055_98614_T1DIslets",
           "HPAP056_98615_Adult-ControlIslets",
           "HPAP057_98616_T2DIslets",
           "HPAP058_98617_T2DIslets",
           "HPAP059_98618_Adult-ControlIslets",
           "HPAP061_99868_T2DIslets",
           "HPAP062_99869_T2DIslets",
           "HPAP063_99870_Adult-ControlIslets",
           "HPAP064_99871_T1DIslets",
           "HPAP065_99872_T2DIslets",
           "HPAP067_99873_Adult-ControlIslets",
           "HPAP069_99874_Adult-ControlIslets",
           "HPAP070_99875_T2DIslets",
           "HPAP071_103407_T1DIslets",
           "HPAP072_103408_pre-t1DIslets",
           "HPAP074_103409_Adult-ControlIslets",
           "HPAP075_103410_Adult-ControlIslets",
           "HPAP077_103411_Adult-ControlIslets",
           "HPAP079_103412_T2DIslets",
           "HPAP080_103413_Adult-ControlIslets",
           "HPAP081_103414_T2DIslets",
           "HPAP082_103415_Adult-ControlIslets",
           "HPAP083_103416_T2DIslets",
           "HPAP084_105149_T1DIslets",
           "HPAP085_103417_T2DIslets",
           "HPAP087_103418_T1DIslets",
           "HPAP088_105150_T2DIslets",
           "HPAP090_105151_T2DIslets",
           "HPAP091_105152_T2DIslets",
           "HPAP092_105153_Adult-ControlIslets",
           "HPAP093_105154_Adult-ControlIslets",
           "HPAP099_106838_Adult-ControlIslets",
           "HPAP100_106839_T2DIslets",
           "HPAP101_106840_Adult-ControlIslets",
           "HPAP103_106841_Adult-ControlIslets",
           "HPAP104_107485_Child-ControlIslets",
           "HPAP105_107486_Adult-ControlIslets",
           "HPAP106_107487_T2DIslets",
           "HPAP107_109359_pre-t1DIslets",
           "HPAP109_109361_T2DIslets"),
  to = c("t1d",      #HPAP021_75162_T1DIslets
         "nd", #HPAP022_75751_Adult-ControlIslets
         "t1d",      #HPAP023_75798_T1DIslets
         "pre_t1d",     #HPAP024_75799_Adult-ControlIslets
         "nd", #HPAP026_77407_Adult-ControlIslets
         "nd", #HPAP027_78239_Adult-ControlIslets
         "t1d", #HPAP028_78240_T1DIslets
         "pre_t1d", #HPAP029_78242_Adult-ControlIslets
         "t1d", #HPAP032_79863_T1DIslets
         "nd", #HPAP034_81423_Child-ControlIslets
         "nd", #HPAP035_81832_Adult-ControlIslets
         "nd", #HPAP036_81833_Adult-ControlIslets
         "nd", #HPAP037_82543_Adult-ControlIslets
         "pre_t1d", #HPAP038_84764_Child-ControlIslets
         "nd", #HPAP039_86672_Child-ControlIslets
         "nd", #HPAP040_86673_Adult-ControlIslets
         "nd", #HPAP042_86674_Child-ControlIslets
         "nd", #HPAP043_86675_Child-ControlIslets
         "nd", #HPAP044_86676_Child-ControlIslets
         "pre_t1d", #HPAP045_86677_Adult-ControlIslets
         "nd", #HPAP047_97169_Child-ControlIslets
         "pre_t1d", #HPAP049_97170_Adult-ControlIslets
         "pre_t1d", #HPAP050_97171_Adult-ControlIslets
         "t2d", #HPAP051_97173_T2DIslets
         "nd", #HPAP052_98611_Adult-ControlIslets
         "nd", #HPAP053_98612_Adult-ControlIslets
         "nd", #HPAP054_98613_Adult-ControlIslets
         "t1d", #HPAP055_98614_T1DIslets
         "nd", #HPAP056_98615_Adult-ControlIslets
         "t2d", #HPAP057_98616_T2DIslets
         "t2d", #HPAP058_98617_T2DIslets
         "nd", #HPAP059_98618_Adult-ControlIslets
         "t2d", #HPAP061_99868_T2DIslets
         "t2d", #HPAP062_99869_T2DIslets
         "nd", #HPAP063_99870_Adult-ControlIslets
         "t1d", #HPAP064_99871_T1DIslets
         "t2d", #HPAP065_99872_T2DIslets
         "nd", #HPAP067_99873_Adult-ControlIslets
         "nd", #HPAP069_99874_Adult-ControlIslets
         "t2d", #HPAP070_99875_T2DIslets
         "t1d", #HPAP071_103407_T1DIslets
         "pre_t1d", #HPAP072_103408_pre-t1DIslets
         "nd", #HPAP074_103409_Adult-ControlIslets
         "nd", #HPAP075_103410_Adult-ControlIslets
         "nd", #HPAP077_103411_Adult-ControlIslets
         "t2d", #HPAP079_103412_T2DIslets
         "nd", #HPAP080_103413_Adult-ControlIslets
         "t2d", #HPAP081_103414_T2DIslets
         "nd", #HPAP082_103415_Adult-ControlIslets
         "t2d", #HPAP083_103416_T2DIslets
         "t1d", #HPAP084_105149_T1DIslets
         "t2d", #HPAP085_103417_T2DIslets
         "t1d", #HPAP087_103418_T1DIslets
         "t2d", #HPAP088_105150_T2DIslets
         "t2d", #HPAP090_105151_T2DIslets
         "t2d", #HPAP091_105152_T2DIslets
         "pre_t1d", #HPAP092_105153_Adult-ControlIslets
         "nd", #HPAP093_105154_Adult-ControlIslets
         "nd", #HPAP099_106838_Adult-ControlIslets
         "t2d", #HPAP100_106839_T2DIslets
         "nd", #HPAP101_106840_Adult-ControlIslets
         "nd", #HPAP103_106841_Adult-ControlIslets
         "nd", #HPAP104_107485_Child-ControlIslets
         "nd", #HPAP105_107486_Adult-ControlIslets
         "t2d", #HPAP106_107487_T2DIslets
         "t1d", #HPAP107_109359_pre-t1DIslets
         "t2d" #HPAP109_109361_T2DIslets
  )
)

Idents(HPAP) <- "cell_typev2.0"
table(HPAP@meta.data[["cell_typev2.0"]])

# Data set for DO THE DISEASE TYPE SORTING
Idents(HPAP) <- "cell_type"
pancreas.integrated <- subset(HPAP, idents = c("Beta", "Alpha", "Delta", "PP_Gamma", "Epsilon", "Ductal", "Acinar", "Immune", "Stellates_Mesenchymal", "Endothelial"))
pancreas.integrated$cell_type <- pancreas.integrated@active.ident
table(pancreas.integrated@meta.data[["cell_typev2.0"]])
table(pancreas.integrated@meta.data[["cell_type"]])
table(pancreas.integrated@meta.data[["sex"]])
table(pancreas.integrated@meta.data[["ancestry"]])
table(pancreas.integrated@meta.data[["disease"]])

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
my_levels <- c("Beta", "Alpha", "Delta", "PP_Gamma", "Epsilon", "Ductal", "Acinar", "Immune", "Stellates_Mesenchymal", "Endothelial")
pancreas.integrated@meta.data$cell_type <- factor(x = pancreas.integrated@meta.data$cell_type, levels = my_levels)
table(pancreas.integrated@meta.data$cell_type)

# Save file
#saveRDS(pancreas.integrated, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\HPA Refrence\Data files\Seurat Objects\pancreas.integrated.rds)")
#pancreas.integrated <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\HPA Refrence\Data files\Seurat Objects\pancreas.integrated.rds)")

#Create jalab file
Idents(pancreas.integrated) <- "cell_typev2.0"
jalab <- subset(pancreas.integrated, idents = c("nd", "GADA", "t1d"))
table(jalab@meta.data[["cell_typev2.0"]])
Idents(jalab) <- "cell_type"

# Save file
#saveRDS(jalab, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\HPA Refrence\Data files\Seurat Objects\jalab.rds)")
jalab <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\HPA Refrence\Data files\Seurat Objects\jalab.rds)")
pancreas.integrated <- jalab

# Plotting
table(pancreas.integrated@meta.data[["cell_typev2.0"]])
table(pancreas.integrated@meta.data[["cell_type"]])
table(pancreas.integrated@meta.data[["sex"]])
table(pancreas.integrated@meta.data[["ancestry"]])
table(pancreas.integrated@meta.data[["disease"]])
Idents(pancreas.integrated) <- "cell_typev2.0"
DimPlot(pancreas.integrated, reduction = "umap", raster = FALSE)

DefaultAssay(object = pancreas.integrated) <- "RNA"
markers.to.plot <- c("INS", "IAPP", "PDX1", "MAFA", "MAFB", "GCG", "DPP4", "GC", "LEPR", "SST", "GHRL", "FRZB", "PPY", "THSD7A",
                     "CFTR", "MMP7", "CELA2A", "CELA3A", "RGS5", "FABP4", "FMOD", "COL3A1", "ENG", "VWF", "CSF1R", "SDS")

# Re select organized idents for DIMPLOT
#Idents(pancreas.integrated) <- "celltype.sample"
#Idents(pancreas.integrated) <- "cell_type"
#DefaultAssay(object = pancreas.integrated) <- "SCT"
DotPlot(pancreas.integrated,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -2, #minimum level
        col.max = 3,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

# Advanced coding for ggplot2
# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = pancreas.integrated) <- "cell_type"
pancreas.integrated$celltype.disease <- paste(Idents(pancreas.integrated),pancreas.integrated$cell_typev2.0, sep = "_")
table(pancreas.integrated@meta.data$celltype.disease)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("Beta_nd", "Beta_GADA", "Beta_t1d", 
                "Alpha_nd", "Alpha_GADA", "Alpha_t1d",
                "Delta_nd", "Delta_GADA", "Delta_t1d",
                "Epsilon_nd", "Epsilon_GADA", "Epsilon_t1d",
                "PP_Gamma_nd", "PP_Gamma_GADA", "PP_Gamma_t1d",
                "Ductal_nd", "Ductal_GADA", "Ductal_t1d",
                "Acinar_nd", "Acinar_GADA", "Acinar_t1d",
                "Stellates_Mesenchymal_nd", "Stellates_Mesenchymal_GADA", "Stellates_Mesenchymal_t1d",
                "Endothelial_nd", "Endothelial_GADA", "Endothelial_t1d",
                "Immune_nd", "Immune_GADA", "Immune_t1d"
                )
head(pancreas.integrated@meta.data$celltype.disease)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$celltype.disease <- factor(x = pancreas.integrated@meta.data$celltype.disease, levels = my_levels2)
table(pancreas.integrated@meta.data$celltype.disease)

# Re select organized idents
Idents(pancreas.integrated) <- "celltype.disease"
DefaultAssay(object = pancreas.integrated) <- "SCT"

# Visualization
#Dotplot
markers.to.plot <- c("INS", "IAPP", "PDX1", "MAFA", "GCG", "DPP4", "GC", "LEPR", "SST", "GHRL", "FRZB", "PPY", "THSD7A",
                     "CFTR", "MMP7", "CELA2A", "CELA3A", "COL1A1", "PDGFRB", "ENG", "VWF", "CSF1R", "SDS")
DotPlot(pancreas.integrated,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -2, #minimum level
        col.max = 3,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

# Gene List
markers.to.plot <- c("HIF1A", "HIF3A", "HIF1AN", "ARNT", "VHL", 
                     "EGLN2", "EGLN1", "EGLN3", #PHD1, 2 and 3
                     "VEGFA", "VEGFB", "VEGFC", "VEGFD",
                     "NRP1", "NRP2", 
                     "KDR", "FLT1", "FLT4",
                     "PDGFA", "PDGFB", "PDGFC", "PDGFD",
                     "PDGFRA", "PDGFRB", "PDAP1",
                     "ADM", "ANGPTL4", "ANGPT2",
                     "ENPP2", "MIF",
                     "PFKP", "GYS1", "LDHA", "PGK1", "PDK4", "PKM", "PFKFB2", "ACOT7", "ALDOA", "PGAM1", "ENO1", "TPI1", "MRPS17",
                     "SLC2A1", "SLC2A3", 
                     "P4HA1", "P4HA2", "VASN", 
                     "NDRG1", "BNIP3L", "MXI1", "CDKN3",
                     "LOX",
                     "TUBB6",
                     "NOS1", "NOS2", "NOS3", "NOA1",
                     "EDN1", "ECE1", "AGT", "REN", "ACE", "ACE2", "AGTRAP",
                     "SLC18A2", "TPH1",
                     "ADA", "ADA2", "SLC17A9", "SLC35B3", "SLC29A1", "SLC28A2", "ENTPD2", "ENTPD5", "ENTPD3", "ENTPD6", "ENTPD1", "ENTPD4", "NT5E"
                     )
DotPlot(pancreas.integrated,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -2, #minimum level
        col.max = 3,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

# Subset Stellate cells
Idents(pancreas.integrated) <- "cell_type"
pancreas.integrated.stellate <- subset(pancreas.integrated, idents = c("Stellates_Mesenchymal"))
Idents(pancreas.integrated.stellate) <- "cell_typev2.0"
pancreas.integrated.stellate@active.ident <- factor(pancreas.integrated.stellate@active.ident, levels=c("nd", "GADA", "t1d"))

# Gene List
DefaultAssay(pancreas.integrated.stellate) <- "SCT"
markers.to.plot <- c("EDNRA", "EDNRB", "FBXW7", 
                     "AGTR1", "AGTR2", "MAS1",
                      "ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C", "ADRB1", "ADRB2", "ADRB3", 
                     "ADORA2B", "ADORA2A", "ADORA3", "ADORA1", 
                     "P2RY1", "P2RY2", "P2RY4", "P2RY6", "P2RY11", "P2RY12", "P2RY13", "P2RY14", "P2RX1", "P2RX2", "P2RX3", "P2RX4", "P2RX5", "P2RX6", "P2RX7", 
                     "PTGDR", "PTGDR2", "PTGER1", "PTGER2", "PTGER3", "PTGER4", "PTGFR", "PTGIR", "TBXA2R", 
                     "HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F", "HTR2A", "HTR2B", "HTR2C", "HTR4", "HTR5A", "HTR5BP", "HTR6", "HTR7", "HTR3A", "HTR3B", "HTR3C", "HTR3D", "HTR3E" 
                     )
DotPlot(pancreas.integrated.stellate,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -2, #minimum level
        col.max = 3,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

#################
# Re select organized idents
Idents(pancreas.integrated) <- "celltype.disease"
DefaultAssay(object = pancreas.integrated) <- "SCT"

# Subset
pancreas.integrated.plot <- subset(pancreas.integrated, 
                                   idents = c("Beta_nd", "Beta_GADA", "Beta_t1d", 
                                              "Alpha_nd", "Alpha_GADA", "Alpha_t1d",
                                              "Endothelial_nd", "Endothelial_GADA", "Endothelial_t1d",
                                              "Stellates_Mesenchymal_nd", "Stellates_Mesenchymal_GADA", "Stellates_Mesenchymal_t1d"
                                              ))
                                   
table(pancreas.integrated.plot@meta.data$celltype.disease)

# Visualization 1
#Dotplot
markers.to.plot <- c("INS", "IAPP", "PDX1", "MAFA", "GCG", "DPP4", "GC", "COL1A1", "PDGFRB", "ENG", "VWF")
DotPlot(pancreas.integrated.plot,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -1, #minimum level
        col.max = 2,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

# Visualization 2
# Re select organized idents
Idents(pancreas.integrated) <- "celltype.disease"
DefaultAssay(object = pancreas.integrated) <- "SCT"
# Subset
pancreas.integrated.plot <- subset(pancreas.integrated, 
                                   idents = c("Stellates_Mesenchymal_nd", "Stellates_Mesenchymal_GADA", "Stellates_Mesenchymal_t1d"
                                   ))

table(pancreas.integrated.plot@meta.data$celltype.disease)
#Dotplot
markers.to.plot <- c("CSPG4", "PDGFRB", "ACTA2", "COL1A1")
DotPlot(pancreas.integrated.plot,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -1, #minimum level
        col.max = 2,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

#Vlnplot
VlnPlot(
  object = pancreas.integrated.plot,
  features = c("CSPG4", "PDGFRB", "ACTA2", "COL1A1"),
  assay = 'SCT',
  #slot = 'counts',
  ncol = 1,
  cols = c("deeppink3",
                      "coral2",
                      "red4",
                      "orange",
                      "indianred",
                      "darkturquoise",
                      "lightgreen",
                      "violet",
                      "purple4",
                      "red"),
                      #y.max = 3,
  pt.size = 1
)



#UMAP
Idents(pancreas.integrated) <- "cell_typev2.0"
Idents(pancreas.integrated) <- "cell_type"
DimPlot(pancreas.integrated, label = FALSE, ncol = 2,  cols = c("deeppink3",
                                                                "coral2",
                                                                "red4",
                                                                "orange",
                                                                "indianred",
                                                                "darkturquoise",
                                                                "lightgreen",
                                                                "violet",
                                                                "purple4",
                                                                "red"
))

#UMAP of genes in ND only
Idents(pancreas.integrated) <- "cell_typev2.0"
pancreas.integrated.nd <- subset(pancreas.integrated, idents = c("nd"))
FeaturePlot(object = pancreas.integrated.nd,
            features = c("CSPG4", "PDGFRB", "ACTA2", "COL1A1"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE)

Idents(pancreas.integrated) <- "cell_typev2.0"
Idents(pancreas.integrated) <- "cell_type"
Beta <- subset(pancreas.integrated, idents = c("Beta"))
table(Beta@meta.data[["cell_typev2.0"]])
Alpha <- subset(pancreas.integrated, idents = c("Alpha"))
table(Alpha@meta.data[["cell_typev2.0"]])
Stellates_Mesenchymal <- subset(pancreas.integrated, idents = c("Stellates_Mesenchymal"))
table(Stellates_Mesenchymal@meta.data[["cell_typev2.0"]])
Endothelial <- subset(pancreas.integrated, idents = c("Endothelial"))
table(Endothelial@meta.data[["cell_typev2.0"]])

# Betacells
Idents(Beta) <- "cell_typev2.0"
Idents(Beta) <- "cell_type"
VlnPlot(
  object = Beta,
  features = c("FN1"),
  assay = 'SCT',
  #slot = 'counts',
  cols = c("deeppink3",
           "coral2",
           "red4",
           "orange",
           "indianred",
           "darkturquoise",
           "lightgreen",
           "violet",
           "purple4",
           "red"),
  #y.max = 3,
  pt.size = 1
)

# Alpha
Idents(Alpha) <- "cell_typev2.0"
Idents(Alpha) <- "cell_type"
VlnPlot(
  object = Alpha,
  features = c("FN1"),
  assay = 'SCT',
  #slot = 'counts',
  cols = c("deeppink3",
           "coral2",
           "red4",
           "orange",
           "indianred",
           "darkturquoise",
           "lightgreen",
           "violet",
           "purple4",
           "red"),
  #y.max = 3,
  pt.size = 1
)

# Endothelial
Idents(Endothelial) <- "cell_typev2.0"
Idents(Endothelial) <- "cell_type"
VlnPlot(
  object = Endothelial,
  features = c("FN1"),
  assay = 'SCT',
  #slot = 'counts',
  cols = c("deeppink3",
           "coral2",
           "red4",
           "orange",
           "indianred",
           "darkturquoise",
           "lightgreen",
           "violet",
           "purple4",
           "red"),
  #y.max = 3,
  pt.size = 1
)

# Stellates_Mesenchymal
Idents(Stellates_Mesenchymal) <- "cell_typev2.0"
Idents(Stellates_Mesenchymal) <- "cell_type"
VlnPlot(
  object = Stellates_Mesenchymal,
  features = c("FN1"),
  assay = 'SCT',
  #slot = 'counts',
  cols = c("deeppink3",
           "coral2",
           "red4",
           "orange",
           "indianred",
           "darkturquoise",
           "lightgreen",
           "violet",
           "purple4",
           "red"),
  #y.max = 3,
  pt.size = 1
)

# DE testing
Idents(pancreas.integrated) <- "celltype.disease"

# 1.Beta-cells (GADA)
beta.ndvsGADA <- FindMarkers(pancreas.integrated,
                             ident.1 = "Beta_GADA", ident.2 = "Beta_nd", # second Ident should be your control
                             test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                             min.pct = 0,
                             logfc.threshold = 0, 
                             pseudocount.use = 1,
                             assay = 'SCT',
                             verbose = TRUE)
head(beta.ndvsGADA, n = 15)
beta.ndvsGADA %>% top_n(n = 2, wt = avg_log2FC)
write.csv(beta.ndvsGADA, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files\DE analysis\raw R files\beta.ndvsGADA.csv)")


# 2.Alpha-cells (GADA)
alpha.ndvsGADA <- FindMarkers(pancreas.integrated,
                             ident.1 = "Alpha_GADA", ident.2 = "Alpha_nd", # second Ident should be your control 
                             test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                             min.pct = 0,
                             logfc.threshold = 0, 
                             pseudocount.use = 1,
                             assay = 'SCT',
                             verbose = TRUE)
head(alpha.ndvsGADA, n = 15)
alpha.ndvsGADA %>% top_n(n = 2, wt = avg_log2FC)
write.csv(alpha.ndvsGADA, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files\DE analysis\raw R files\alpha.ndvsGADA.csv)")

# 3.Mesenchymal cells (GADA)
Stellates_Mesenchymal.ndvsGADA <- FindMarkers(pancreas.integrated,
                              ident.1 = "Stellates_Mesenchymal_GADA", ident.2 = "Stellates_Mesenchymal_nd", # second Ident should be your control 
                              test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                              min.pct = 0,
                              logfc.threshold = 0, 
                              pseudocount.use = 1,
                              assay = 'SCT',
                              verbose = TRUE)
head(Stellates_Mesenchymal.ndvsGADA, n = 15)
Stellates_Mesenchymal.ndvsGADA %>% top_n(n = 2, wt = avg_log2FC)
write.csv(Stellates_Mesenchymal.ndvsGADA, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files\DE analysis\raw R files\Stellates_Mesenchymal.ndvsGADA.csv)")

# 4.Endothelial cells (GADA)
Endothelial.ndvsGADA <- FindMarkers(pancreas.integrated,
                                              ident.1 = "Endothelial_GADA", ident.2 = "Endothelial_nd", # second Ident should be your control 
                                              test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                              min.pct = 0,
                                              logfc.threshold = 0, 
                                              pseudocount.use = 1,
                                              assay = 'SCT',
                                              verbose = TRUE)
head(Endothelial.ndvsGADA, n = 15)
Endothelial.ndvsGADA %>% top_n(n = 2, wt = avg_log2FC)
write.csv(Endothelial.ndvsGADA, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files\DE analysis\raw R files\Endothelial.ndvsGADA.csv)")

# 5.Beta-cells (T1D)
beta.ndvst1d <- FindMarkers(pancreas.integrated,
                             ident.1 = "Beta_t1d", ident.2 = "Beta_nd", # second Ident should be your control
                             test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                             min.pct = 0,
                             logfc.threshold = 0, 
                             pseudocount.use = 1,
                             assay = 'SCT',
                             verbose = TRUE)
head(beta.ndvst1d, n = 15)
beta.ndvst1d %>% top_n(n = 2, wt = avg_log2FC)
write.csv(beta.ndvst1d, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files\DE analysis\raw R files\beta.ndvst1d.csv)")


# 6.Alpha-cells (T1D)
alpha.ndvst1d <- FindMarkers(pancreas.integrated,
                              ident.1 = "Alpha_t1d", ident.2 = "Alpha_nd", # second Ident should be your control 
                              test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                              min.pct = 0,
                              logfc.threshold = 0, 
                              pseudocount.use = 1,
                              assay = 'SCT',
                              verbose = TRUE)
head(alpha.ndvst1d, n = 15)
alpha.ndvst1d %>% top_n(n = 2, wt = avg_log2FC)
write.csv(alpha.ndvst1d, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files\DE analysis\raw R files\alpha.ndvst1d.csv)")

# 7.Mesenchymal cells (T1D)
Stellates_Mesenchymal.ndvst1d <- FindMarkers(pancreas.integrated,
                                              ident.1 = "Stellates_Mesenchymal_t1d", ident.2 = "Stellates_Mesenchymal_nd", # second Ident should be your control 
                                              test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                              min.pct = 0,
                                              logfc.threshold = 0, 
                                              pseudocount.use = 1,
                                              assay = 'SCT',
                                              verbose = TRUE)
head(Stellates_Mesenchymal.ndvst1d, n = 15)
Stellates_Mesenchymal.ndvst1d %>% top_n(n = 2, wt = avg_log2FC)
write.csv(Stellates_Mesenchymal.ndvst1d, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files\DE analysis\raw R files\Stellates_Mesenchymal.ndvst1d.csv)")

# 8.Endothelial cells (T1D)
Endothelial.ndvst1d <- FindMarkers(pancreas.integrated,
                                    ident.1 = "Endothelial_t1d", ident.2 = "Endothelial_nd", # second Ident should be your control 
                                    test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                    min.pct = 0,
                                    logfc.threshold = 0, 
                                    pseudocount.use = 1,
                                    assay = 'SCT',
                                    verbose = TRUE)
head(Endothelial.ndvst1d, n = 15)
Endothelial.ndvst1d %>% top_n(n = 2, wt = avg_log2FC)
write.csv(Endothelial.ndvst1d, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files\DE analysis\raw R files\Endothelial.ndvst1d.csv)")

# 9.Beta-cells (T1D)
beta.t1dvsGADA <- FindMarkers(pancreas.integrated,
                            ident.1 = "Beta_t1d", ident.2 = "Beta_GADA", # second Ident should be your control
                            test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                            min.pct = 0,
                            logfc.threshold = 0, 
                            pseudocount.use = 1,
                            assay = 'SCT',
                            verbose = TRUE)
head(beta.t1dvsGADA, n = 15)
beta.t1dvsGADA %>% top_n(n = 2, wt = avg_log2FC)
write.csv(beta.t1dvsGADA, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files\DE analysis\raw R files\beta.t1dvsGADA.csv)")


# 10.Alpha-cells (T1D)
alpha.t1dvsGADA <- FindMarkers(pancreas.integrated,
                             ident.1 = "Alpha_t1d", ident.2 = "Alpha_GADA", # second Ident should be your control 
                             test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                             min.pct = 0,
                             logfc.threshold = 0, 
                             pseudocount.use = 1,
                             assay = 'SCT',
                             verbose = TRUE)
head(alpha.t1dvsGADA, n = 15)
alpha.t1dvsGADA %>% top_n(n = 2, wt = avg_log2FC)
write.csv(alpha.t1dvsGADA, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files\DE analysis\raw R files\alpha.t1dvsGADA.csv)")

# 11.Mesenchymal cells (T1D)
Stellates_Mesenchymal.t1dvsGADA <- FindMarkers(pancreas.integrated,
                                             ident.1 = "Stellates_Mesenchymal_t1d", ident.2 = "Stellates_Mesenchymal_GADA", # second Ident should be your control 
                                             test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                             min.pct = 0,
                                             logfc.threshold = 0, 
                                             pseudocount.use = 1,
                                             assay = 'SCT',
                                             verbose = TRUE)
head(Stellates_Mesenchymal.t1dvsGADA, n = 15)
Stellates_Mesenchymal.t1dvsGADA %>% top_n(n = 2, wt = avg_log2FC)
write.csv(Stellates_Mesenchymal.t1dvsGADA, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files\DE analysis\raw R files\Stellates_Mesenchymal.t1dvsGADA.csv)")

# 12.Endothelial cells (T1D)
Endothelial.t1dvsGADA <- FindMarkers(pancreas.integrated,
                                   ident.1 = "Endothelial_t1d", ident.2 = "Endothelial_GADA", # second Ident should be your control 
                                   test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                   min.pct = 0,
                                   logfc.threshold = 0, 
                                   pseudocount.use = 1,
                                   assay = 'SCT',
                                   verbose = TRUE)
head(Endothelial.t1dvsGADA, n = 15)
Endothelial.t1dvsGADA %>% top_n(n = 2, wt = avg_log2FC)
write.csv(Endothelial.t1dvsGADA, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files\DE analysis\raw R files\Endothelial.t1dvsGADA.csv)")

# Whole pancreas
# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = pancreas.integrated) <- "cell_type"
pancreas.integrated$celltype.sample <- paste(Idents(pancreas.integrated),pancreas.integrated$sex, sep = "_")
table(pancreas.integrated@meta.data$celltype.sample)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("Beta_f", "Beta_m", "Alpha_f", "Alpha_m", "Delta_f", "Delta_m", "PP_Gamma_f", "PP_Gamma_m", "Epsilon_f", "Epsilon_m",
                "Ductal_m", "Ductal_f", "Acinar_f", "Acinar_m", "Stellates_Mesenchymal_f", "Stellates_Mesenchymal_m", "Endothelial_f", "Endothelial_m", "Immune_f", "Immune_m"
                )
head(pancreas.integrated@meta.data$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$celltype.sample <- factor(x = pancreas.integrated@meta.data$celltype.sample, levels = my_levels2)
table(pancreas.integrated@meta.data$celltype.sample)

# Re select organized idents
Idents(pancreas.integrated) <- "celltype.sample"
Idents(pancreas.integrated) <- "celltype"
DefaultAssay(object = pancreas.integrated) <- "SCT"
DotPlot(pancreas.integrated,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -2, #minimum level
        col.max = 3,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

# Advanced coding for ggplot2
# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = pancreas.integrated) <- "celltype.sample"
pancreas.integrated$celltype.sample.disease <- paste(Idents(pancreas.integrated),pancreas.integrated$disease, sep = "_")
table(pancreas.integrated@meta.data$celltype.sample.disease)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("Beta_f_nd", "Beta_m_nd", "Beta_f_pre_t1d", "Beta_m_pre_t1d", "Beta_f_t1d", "Beta_m_t1d", "Beta_f_t2d", "Beta_m_t2d",
                "Alpha_f_nd", "Alpha_m_nd", "Alpha_f_pre_t1d", "Alpha_m_pre_t1d", "Alpha_f_t1d", "Alpha_m_t1d", "Alpha_f_t2d", "Alpha_m_t2d",
                "Delta_f_nd", "Delta_m_nd", "Delta_f_pre_t1d", "Delta_m_pre_t1d", "Delta_f_t1d", "Delta_m_t1d",  "Delta_f_t2d", "Delta_m_t2d",
                "PP_Gamma_f_nd", "PP_Gamma_m_nd", "PP_Gamma_f_pre_t1d", "PP_Gamma_m_pre_t1d", "PP_Gamma_f_t1d", "PP_Gamma_m_t1d", "PP_Gamma_f_t2d", "PP_Gamma_m_t2d",
                "Epsilon_f_nd", "Epsilon_m_nd", "Epsilon_f_pre_t1d", "Epsilon_m_pre_t1d", "Epsilon_f_t1d", "Epsilon_m_t1d", "Epsilon_f_t2d", "Epsilon_m_t2d",
                "Ductal_m_nd", "Ductal_f_nd", "Ductal_m_pre_t1d", "Ductal_f_pre_t1d", "Ductal_m_t1d", "Ductal_f_t1d", "Ductal_m_t2d", "Ductal_f_t2d",
                "Acinar_f_nd", "Acinar_m_nd", "Acinar_f_pre_t1d", "Acinar_m_pre_t1d", "Acinar_f_t1d", "Acinar_m_t1d", "Acinar_f_t2d", "Acinar_m_t2d",
                "Stellates_Mesenchymal_f_nd", "Stellates_Mesenchymal_m_nd", "Stellates_Mesenchymal_f_pre_t1d", "Stellates_Mesenchymal_m_pre_t1d", "Stellates_Mesenchymal_f_t1d", "Stellates_Mesenchymal_m_t1d", "Stellates_Mesenchymal_f_t2d", "Stellates_Mesenchymal_m_t2d",
                "Endothelial_f_nd", "Endothelial_m_nd", "Endothelial_f_pre_t1d", "Endothelial_m_pre_t1d", "Endothelial_f_t1d", "Endothelial_m_t1d", "Endothelial_f_t2d", "Endothelial_m_t2d",
                "Immune_f_nd", "Immune_m_nd", "Immune_f_pre_t1d", "Immune_m_pre_t1d", "Immune_f_t1d", "Immune_m_t1d", "Immune_f_t2d", "Immune_m_t2d"
                )
                               
head(pancreas.integrated@meta.data$cell_type)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$celltype.sample.disease <- factor(x = pancreas.integrated@meta.data$celltype.sample.disease, levels = my_levels2)
table(pancreas.integrated@meta.data$celltype.sample.disease)

# Re select organized idents
Idents(pancreas.integrated) <- "celltype.sample.disease"
#Idents(pancreas.integrated) <- "celltype"
DefaultAssay(object = pancreas.integrated) <- "SCT"
DotPlot(pancreas.integrated,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -2, #minimum level
        col.max = 3,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))
FeaturePlot(object = pancreas.integrated,
            features = c("INS", "GCG", "SST", "PPY", "GHRL",
                         "KRT19", "CPA1",
                         "COL1A1", "VWF", "SOX10",
                         "TPSAB1", "SDS", "TRAC"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            #min.cutoff = 0,
            #max.cutoff = 1,
            slot = 'counts',
            order = TRUE)

FeaturePlot(object = pancreas.integrated,
            features = c("PGR", "ESR1", "AR"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE)

# JAlab data analysis
# Observing cells
Idents(jalab) <- "cell_typev2.0"
Idents(jalab) <- "cell_type"
DimPlot(jalab, label = FALSE, ncol = 2,  cols = c("deeppink3",
                                                  "coral2",
                                                  "red4",
                                                  "orange",
                                                  "indianred",
                                                  "darkturquoise",
                                                  "lightgreen",
                                                  "violet",
                                                  "purple4",
                                                  "red"
))


table(jalab@meta.data[["cell_typev2.0"]])
table(jalab@meta.data[["disease"]])
DimPlot(jalab, split.by = "cell_typev2.0", group.by = "cell_type", label = FALSE, ncol = 2,  cols = c("red",
                                                                                                            "red4",
                                                                                                            "orange",
                                                                                                            "lightgoldenrod3",
                                                                                                            "sienna",
                                                                                                            "indianred",
                                                                                                            "orangered1",
                                                                                                            "black",
                                                                                                            "darkturquoise",
                                                                                                            "paleturquoise",
                                                                                                            "lightgreen",
                                                                                                            "springgreen4",
                                                                                                            "darkolivegreen",
                                                                                                            "purple4",
                                                                                                            "purple",
                                                                                                            "deeppink",
                                                                                                            "violetred",
                                                                                                            "violet"
))




DimPlot(HPAP, reduction = "umap")

DimPlot(object, 
        order = c("ident1", "ident2", "ident3", "HPAP"))

################################################
################################################

# Mesenchymal ####
# Look at your default assay
DefaultAssay(object = Stellates_Mesenchymal)

# Change default assay to integrated, save information in the "integrated" assay
DefaultAssay(object = Stellates_Mesenchymal) <- "SCT"

# PCA analysis data will be stored in the "reductions' slot
Stellates_Mesenchymal <- RunPCA(object = Stellates_Mesenchymal, features = VariableFeatures(object = Stellates_Mesenchymal))

# CLUSTER ANALYSIS ####
# Clustering, we run multiple permutations to allow clustree to analyze optimal clustering resolution.
Stellates_Mesenchymal <- FindNeighbors(object = Stellates_Mesenchymal, dims = 1:20)
Stellates_Mesenchymal <- FindClusters(object = Stellates_Mesenchymal, resolution = 0)
Stellates_Mesenchymal <- FindClusters(object = Stellates_Mesenchymal, resolution = 0.1)
Stellates_Mesenchymal <- FindClusters(object = Stellates_Mesenchymal, resolution = 0.2)
Stellates_Mesenchymal <- FindClusters(object = Stellates_Mesenchymal, resolution = 0.3)
Stellates_Mesenchymal <- FindClusters(object = Stellates_Mesenchymal, resolution = 0.4)
Stellates_Mesenchymal <- FindClusters(object = Stellates_Mesenchymal, resolution = 0.5)
Stellates_Mesenchymal <- FindClusters(object = Stellates_Mesenchymal, resolution = 0.6)

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution 
clustree(Stellates_Mesenchymal, prefix = "SCT_snn_res.")

# We run our final, chosen clustering resolution, so that seurat_clusters an be populated for further analysis
Stellates_Mesenchymal <- FindClusters(object = Stellates_Mesenchymal, resolution = 0.1)

# Check data
table(Stellates_Mesenchymal$seurat_clusters)
table(Idents(Stellates_Mesenchymal), Stellates_Mesenchymal$orig.ident)

# NON-LINEAR DIMENSIONALITY REDUCTION ####
# RunUMAP
Stellates_Mesenchymal <- RunUMAP(Stellates_Mesenchymal, dims = 1:30) 

#Plotting
DimPlot(object = Stellates_Mesenchymal, reduction = "umap", pt.size = 2, cols = c("royalblue1",
                                                                                  "red3"
                                                                                  ))

# UMAP
Idents(pancreas.integrated) <- "celltype.disease"
Stellates_Mesenchymal_nd <- subset(pancreas.integrated, idents = c("Stellates_Mesenchymal_nd"))
FeaturePlot(object = Stellates_Mesenchymal, 
            features = c("PDGFRB", "ACTA2", "CSPG4", "COL1A1",
                         "PDGFRA", "RGS5", "MKI67", "ACE2", "EDNRA", "ADRA1A"),
            pt.size = 0.5,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 2,
            order = TRUE)

Idents(Stellates_Mesenchymal) <- "cell_typev2.0"
Idents(Stellates_Mesenchymal) <- "cell_type"
table(Stellates_Mesenchymal@meta.data[["cell_typev2.0"]])

#Rename Idents
Idents(Stellates_Mesenchymal) <- "SCT_snn_res.0.1"
Stellates_Mesenchymal <- RenameIdents(Stellates_Mesenchymal, 
                                      "0" = "Activated-Stellate", 
                                      "1" = "Pericyte",
                                      "2" = "Activated-Stellate", 
                                      "3" = "Activated-Stellate"
)

# Saving this information in the metadata slot
table(Idents(Stellates_Mesenchymal))
Stellates_Mesenchymal$celltype.mes <- Idents(Stellates_Mesenchymal)
table(Stellates_Mesenchymal@meta.data$celltype.mes)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("Activated-Stellate", "Pericyte")
head(Stellates_Mesenchymal@meta.data$celltype.mes)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
Stellates_Mesenchymal@meta.data$celltype.mes <- factor(x = Stellates_Mesenchymal@meta.data$celltype.mes, levels = my_levels)

# Advanced coding for ggplot2
# Stellate cells only
# Create a new metadata slot containing combined info, segregating clusters and samples on basis of SEX
Idents(object = Stellates_Mesenchymal) <- "celltype.mes"
Stellates_Mesenchymal$celltype.mes.sex <- paste(Idents(Stellates_Mesenchymal),Stellates_Mesenchymal$sex, sep = "_")
table(Stellates_Mesenchymal@meta.data$celltype.mes.sex)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("Activated-Stellate_f", "Activated-Stellate_m", "Pericyte_f", "Pericyte_m"
)
table(Stellates_Mesenchymal@meta.data$celltype.mes.sex)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
Stellates_Mesenchymal@meta.data$celltype.mes.sex <- factor(x = Stellates_Mesenchymal@meta.data$celltype.mes.sex, levels = my_levels2)
table(Stellates_Mesenchymal@meta.data$celltype.mes.sex)

# Create a new metadata slot containing combined info, segregating clusters and samples on basis of cell_type2.0
Idents(object = Stellates_Mesenchymal) <- "celltype.mes"
Stellates_Mesenchymal$celltype.mes.celltype <- paste(Idents(Stellates_Mesenchymal),Stellates_Mesenchymal$cell_typev2.0, sep = "_")
table(Stellates_Mesenchymal@meta.data$celltype.mes.celltype)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("Activated-Stellate_nd", "Activated-Stellate_GADA", "Activated-Stellate_t1d", 
                "Pericyte_nd", "Pericyte_GADA", "Pericyte_t1d"
)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
Stellates_Mesenchymal@meta.data$celltype.mes.celltype <- factor(x = Stellates_Mesenchymal@meta.data$celltype.mes.celltype, levels = my_levels2)
table(Stellates_Mesenchymal@meta.data$celltype.mes.celltype)

# Selected genes
markers.to.plot <- c("PDGFRB", "ACTA2", "CSPG4", "COL1A1",
                     "PDGFRA", "RGS5", "MKI67", "ACE2", "EDNRA", "ADRA1A")

# Re select organized idents
Idents(Stellates_Mesenchymal) <- "celltype.mes.celltype"
table(Stellates_Mesenchymal@meta.data$celltype.mes.celltype)
#Idents(Stellates_Mesenchymal) <- "celltype"
DefaultAssay(object = Stellates_Mesenchymal) <- "SCT"
DotPlot(Stellates_Mesenchymal,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -2, #minimum level
        col.max = 3,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))


# New metadata column is not paired, so we need to pair
my_levels3 <- c("nd", "pre_t1d", "t1d"
                )
table(pancreas.integrated.stellate@meta.data$disease)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated.stellate@meta.data$disease <- factor(x = pancreas.integrated.stellate@meta.data$disease, levels = my_levels3)
table(pancreas.integrated.stellate@meta.data$disease)

# Selected genes
markers.to.plot <- c("PDGFRB", "ACTA2", "CSPG4", "MCAM", "RGS5",
                     "PDGFRA", "MYH11", "PXN", "TNS1", "TNS2", "TNS3", "TNS4", "FN1", "POSTN")

markers.to.plot <- c("MAS1", "AGTR2", "AGTR1", "EDNRB", "EDNRA")

# Re select organized idents
Idents(pancreas.integrated.stellate) <- "disease"
table(pancreas.integrated.stellate@meta.data$celltype.mes.celltype)
#Idents(Stellates_Mesenchymal) <- "celltype"
DefaultAssay(object = pancreas.integrated.stellate) <- "SCT"
DotPlot(pancreas.integrated.stellate,  
        dot.scale = 8,
        scale = TRUE,
        col.min = -2, #minimum level
        col.max = 3,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

#ENDOTHELIUM ####
# New metadata column is not paired, so we need to pair
my_levels3 <- c("nd", "pre_t1d", "t1d"
)
table(Endothelial@meta.data$disease)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
Endothelial@meta.data$disease <- factor(x = Endothelial@meta.data$disease, levels = my_levels3)
table(Endothelial@meta.data$disease)
# Selected genes
markers.to.plot <- c("MME", "ECE2", "ECE1", "EDN3", "EDN2", "EDN1")

# Re select organized idents
Idents(Endothelial) <- "disease"
table(Endothelial@meta.data$celltype)
#Idents(Stellates_Mesenchymal) <- "celltype"
DefaultAssay(object = Endothelial) <- "SCT"
DotPlot(Endothelial,  
        dot.scale = 10,
        scale = TRUE,
        col.min = -2, #minimum level
        col.max = 3,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

# DE testing
Idents(Stellates_Mesenchymal) <- "celltype.mes"
table(Idents(Stellates_Mesenchymal))

# 1.Pericytes All cells
PericytevsActivated.Stellate <- FindMarkers(Stellates_Mesenchymal,
                             ident.1 = "Pericyte", ident.2 = "Activated-Stellate", # second Ident should be your control
                             test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                             min.pct = 0,
                             logfc.threshold = 0, 
                             pseudocount.use = 1,
                             assay = 'SCT',
                             verbose = TRUE)
head(PericytevsActivated.Stellate, n = 15)
PericytevsActivated.Stellate %>% top_n(n = 2, wt = avg_log2FC)
write.csv(PericytevsActivated.Stellate, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files_Oct3rd\DE analysis\raw R files\PericytevsActivated.Stellate.csv)")

# 2.Pericytes (GADAvsND)
Idents(Stellates_Mesenchymal) <- "celltype.mes.celltype"
table(Idents(Stellates_Mesenchymal))

Pericyte.GADAvsnd <- FindMarkers(Stellates_Mesenchymal,
                             ident.1 = "Pericyte_GADA", ident.2 = "Pericyte_nd", # second Ident should be your control
                             test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                             min.pct = 0,
                             logfc.threshold = 0, 
                             pseudocount.use = 1,
                             assay = 'SCT',
                             verbose = TRUE)
head(Pericyte.GADAvsnd, n = 15)
Pericyte.GADAvsnd %>% top_n(n = 2, wt = avg_log2FC)
write.csv(Pericyte.GADAvsnd, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files_Oct3rd\DE analysis\raw R files\Pericyte.GADAvsnd.csv)")


# 3.Pericytes (T1DvsGADA)
Pericyte.T1DvsGADA <- FindMarkers(Stellates_Mesenchymal,
                                 ident.1 = "Pericyte_t1d", ident.2 = "Pericyte_GADA", # second Ident should be your control
                                 test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                 min.pct = 0,
                                 logfc.threshold = 0, 
                                 pseudocount.use = 1,
                                 assay = 'SCT',
                                 verbose = TRUE)
head(Pericyte.T1DvsGADA, n = 15)
Pericyte.T1DvsGADA %>% top_n(n = 2, wt = avg_log2FC)
write.csv(Pericyte.T1DvsGADA, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files_Oct3rd\DE analysis\raw R files\Pericyte.T1DvsGADA.csv)")

# 3.Pericytes (T1DvsND)
Pericyte.T1Dvsnd <- FindMarkers(Stellates_Mesenchymal,
                                 ident.1 = "Pericyte_t1d", ident.2 = "Pericyte_nd", # second Ident should be your control
                                 test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                 min.pct = 0,
                                 logfc.threshold = 0, 
                                 pseudocount.use = 1,
                                 assay = 'SCT',
                                 verbose = TRUE)
head(Pericyte.T1Dvsnd, n = 15)
Pericyte.T1Dvsnd %>% top_n(n = 2, wt = avg_log2FC)
write.csv(Pericyte.T1Dvsnd, r"(C:\Users\mqadir\Box\Fahd shared to JA+LG\Data Files_Oct3rd\DE analysis\raw R files\Pericyte.T1Dvsnd.csv)")

#########################################################
########################################################

VlnPlot(
  object = Stellates_Mesenchymal,
  features = c("ACTA2"),
  assay = 'SCT',
  #slot = 'counts',
  cols = c("deeppink3",
           "coral2",
           "red4",
           "orange",
           "indianred",
           "darkturquoise",
           "lightgreen",
           "violet",
           "purple4",
           "red"),
  #y.max = 3,
  pt.size = 1
)

FeaturePlot(object = Stellates_Mesenchymal,
            features = c("MYH11"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE)


# Validating Data congruence
pancreas.integrated$checking <- paste(pancreas.integrated$sample_id, pancreas.integrated$cell_typev2.0, sep = "_")
table(pancreas.integrated$checking)

# Output
 HPAP024_75799_Adult-ControlIslets_GADA    HPAP026_77407_Adult-ControlIslets_nd  HPAP029_78242_Adult-ControlIslets_GADA             HPAP032_79863_T1DIslets_t1d 
                                   1860                                     545                                    2466                                    2214 
   HPAP034_81423_Child-ControlIslets_nd    HPAP035_81832_Adult-ControlIslets_nd    HPAP036_81833_Adult-ControlIslets_nd    HPAP037_82543_Adult-ControlIslets_nd 
                                    342                                    1975                                    1880                                    3559 
 HPAP038_84764_Child-ControlIslets_GADA    HPAP040_86673_Adult-ControlIslets_nd  HPAP045_86677_Adult-ControlIslets_GADA    HPAP047_97169_Child-ControlIslets_nd 
                                   2187                                    2080                                    3150                                    2360 
 HPAP050_97171_Adult-ControlIslets_GADA             HPAP064_99871_T1DIslets_t1d            HPAP071_103407_T1DIslets_t1d       HPAP072_103408_pre-t1DIslets_GADA 
                                   2200                                    2016                                    2337                                    1132 
  HPAP082_103415_Adult-ControlIslets_nd            HPAP084_105149_T1DIslets_t1d HPAP092_105153_Adult-ControlIslets_GADA   HPAP099_106838_Adult-ControlIslets_nd 
                                   1690                                    2246                                    3789    

# Output is correct based on manual checking with HPAP on 07/05/2023

# Session info
> sessionInfo()
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggseqlogo_0.1                     BiocParallel_1.32.5               cicero_1.3.9                      Gviz_1.40.1                      
 [5] SeuratWrappers_0.3.1              chromVAR_1.5.0                    motifmatchr_1.18.0                BSgenome.Hsapiens.UCSC.hg38_1.4.4
 [9] BSgenome_1.64.0                   rtracklayer_1.56.1                Biostrings_2.64.1                 XVector_0.38.0                   
[13] TFBSTools_1.34.0                  JASPAR2020_0.99.10                qs_0.25.5                         R.utils_2.12.2                   
[17] R.oo_1.25.0                       R.methodsS3_1.8.2                 devtools_2.4.5                    usethis_2.1.6                    
[21] ggVennDiagram_1.2.2               ggvenn_0.1.9                      DropletUtils_1.16.0               Nebulosa_1.6.0                   
[25] scCustomize_1.1.1                 circlize_0.4.15                   ComplexHeatmap_2.12.1             viridis_0.6.2                    
[29] viridisLite_0.4.2                 EnrichmentBrowser_2.26.0          graph_1.74.0                      escape_1.6.0                     
[33] dittoSeq_1.8.1                    DOSE_3.22.1                       clusterProfiler_4.4.4             MeSHDbi_1.32.0                   
[37] AnnotationHub_3.4.0               BiocFileCache_2.4.0               dbplyr_2.3.0                      org.Hs.eg.db_3.15.0              
[41] GOSemSim_2.22.0                   glmGamPoi_1.8.0                   EnhancedVolcano_1.14.0            DoubletFinder_2.0.3              
[45] future_1.32.0                     patchwork_1.1.2                   clustree_0.5.0                    ggraph_2.1.0                     
[49] plotly_4.10.1                     EnsDb.Hsapiens.v86_2.99.0         ensembldb_2.20.2                  AnnotationFilter_1.20.0          
[53] GenomicFeatures_1.48.4            AnnotationDbi_1.58.0              scDblFinder_1.10.0                Signac_1.10.0                    
[57] harmony_0.1.1                     monocle3_1.3.1                    SingleCellExperiment_1.20.0       SeuratObject_4.1.3               
[61] Seurat_4.3.0                      reticulate_1.28                   data.table_1.14.6                 forcats_1.0.0                    
[65] purrr_1.0.1                       readr_2.1.4                       tidyr_1.3.0                       tibble_3.1.8                     
[69] tidyverse_1.3.2                   dplyr_1.1.0                       ggridges_0.5.4                    Matrix_1.5-3                     
[73] cowplot_1.1.1                     Rcpp_1.0.10                       SoupX_1.6.2                       hdf5r_1.3.8                      
[77] stringr_1.5.0                     leiden_0.4.3                      ggrepel_0.9.3                     ggplot2_3.4.2                    
[81] DESeq2_1.36.0                     SummarizedExperiment_1.28.0       Biobase_2.58.0                    MatrixGenerics_1.10.0            
[85] matrixStats_0.63.0                GenomicRanges_1.50.2              GenomeInfoDb_1.37.1               IRanges_2.32.0                   
[89] S4Vectors_0.36.1                  BiocGenerics_0.44.0              

loaded via a namespace (and not attached):
  [1] KEGGREST_1.36.3               genefilter_1.78.0             locfit_1.5-9.7                remotes_2.4.2                 BiocVersion_3.15.2           
  [6] lattice_0.20-45               paletteer_1.5.0               spatstat.utils_3.0-1          vctrs_0.5.2                   utf8_1.2.3                   
 [11] blob_1.2.4                    withr_2.5.0                   foreign_0.8-84                ggalt_0.4.0                   readxl_1.4.2                 
 [16] lifecycle_1.0.3               cellranger_1.1.0              munsell_0.5.0                 ScaledMatrix_1.6.0            codetools_0.2-19             
 [21] lmtest_0.9-40                 msigdbr_7.5.1                 limma_3.54.1                  DO.db_2.9                     magick_2.7.4                 
 [26] annotate_1.74.0               parallelly_1.36.0             fs_1.6.1                      fastmatch_1.1-3               metapod_1.4.0                
 [31] Rtsne_0.16                    biovizBase_1.44.0             stringi_1.7.12                RcppRoll_0.3.0                sctransform_0.3.5.9002       
 [36] polyclip_1.10-4               rhdf5filters_1.10.0           yulab.utils_0.0.6             goftest_1.2-3                 cluster_2.1.4                
 [41] ape_5.6-2                     TFMPvalue_0.0.9               pkgconfig_2.0.3               pheatmap_1.0.12               prettyunits_1.1.1            
 [46] sparseMatrixStats_1.10.0      googledrive_2.1.1             lubridate_1.9.2               timechange_0.2.0              httr_1.4.6                   
 [51] igraph_1.3.5                  treeio_1.20.2                 progress_1.2.2                GetoptLong_1.0.5              terra_1.7-3                  
 [56] beachmat_2.14.0               graphlayouts_0.8.4            haven_2.5.1                   snakecase_0.11.0              ggfun_0.1.1                  
 [61] htmltools_0.5.4               miniUI_0.1.1.1                yaml_2.3.7                    pillar_1.9.0                  later_1.3.0                  
 [66] fitdistrplus_1.1-11           glue_1.6.2                    DBI_1.1.3                     plyr_1.8.8                    foreach_1.5.2                
 [71] ProtGenerics_1.28.0           gtable_0.3.3                  rsvd_1.0.5                    caTools_1.18.2                GlobalOptions_0.1.2          
 [76] latticeExtra_0.6-30           extrafont_0.19                fastmap_1.1.0                 broom_1.0.5                   checkmate_2.1.0              
 [81] promises_1.2.0.1              ggforce_0.4.1                 hms_1.1.3                     png_0.1-8                     ash_1.0-15                   
 [86] clue_0.3-64                   ggtree_3.4.4                  spatstat.explore_3.0-6        lazyeval_0.2.2                Formula_1.2-5                
 [91] profvis_0.3.7                 extrafontdb_1.0               crayon_1.5.2                  reprex_2.0.2                  boot_1.3-28.1                
 [96] tidyselect_1.2.0              xfun_0.39                     ks_1.14.0                     BiocSingular_1.14.0           VariantAnnotation_1.42.1     
[101] splines_4.2.2                 survival_3.5-0                RVenn_1.1.0                   rappdirs_0.3.3                xgboost_1.7.3.1              
[106] bit64_4.0.5                   modelr_0.1.11                 CNEr_1.32.0                   jpeg_0.1-10                   stringfish_0.15.7            
[111] VGAM_1.1-8                    htmlTable_2.4.1               xtable_1.8-4                  googlesheets4_1.1.1           DT_0.28                      
[116] cachem_1.0.6                  DelayedArray_0.24.0           abind_1.4-5                   mime_0.12                     nabor_0.5.0                  
[121] rjson_0.2.21                  aplot_0.1.10                  processx_3.8.0                spatstat.sparse_3.0-0         tools_4.2.2                  
[126] cli_3.6.0                     magrittr_2.0.3                proxy_0.4-27                  dichromat_2.0-0.1             future.apply_1.11.0          
[131] ggplotify_0.1.1               DelayedMatrixStats_1.20.0     ggbeeswarm_0.7.2              assertthat_0.2.1              qvalue_2.28.0                
[136] fgsea_1.22.0                  ggprism_1.0.4                 janitor_2.2.0                 HDF5Array_1.26.0              ica_1.0-3                    
[141] pbapply_1.7-2                 ggrastr_1.0.2                 scuttle_1.8.4                 tweenr_2.0.2                  Rgraphviz_2.40.0             
[146] zlibbioc_1.44.0               restfulr_0.0.15               RApiSerialize_0.1.2           biomaRt_2.52.0                shadowtext_0.1.2             
[151] tzdb_0.3.0                    geneplotter_1.74.0            ps_1.7.2                      fansi_1.0.4                   tidygraph_1.2.3              
[156] GSEABase_1.58.0               UCell_2.0.1                   tensor_1.5                    ROCR_1.0-11                   KernSmooth_2.23-20           
[161] backports_1.4.1               interp_1.1-3                  farver_2.1.1                  bit_4.0.5                     Rsamtools_2.12.0             
[166] proj4_1.0-12                  RANN_2.6.1                    shiny_1.7.4                   BiocIO_1.6.0                  scattermore_0.8              
[171] scatterpie_0.2.1              RcppAnnoy_0.0.20              maps_3.4.1                    downloader_0.4                KEGGgraph_1.56.0             
[176] rstudioapi_0.14               minqa_1.2.5                   iterators_1.0.14              Rhdf5lib_1.20.0               spatstat.geom_3.0-6          
[181] nlme_3.1-162                  DirichletMultinomial_1.38.0   shape_1.4.6                   gtools_3.9.4                  beeswarm_0.4.0               
[186] rematch2_2.1.2                sf_1.0-9                      listenv_0.9.0                 reshape2_1.4.4                rhdf5_2.42.0                 
[191] gargle_1.5.1                  GSVA_1.44.5                   generics_0.1.3                colorspace_2.1-0              base64enc_0.1-3              
[196] XML_3.99-0.13                 pkgbuild_1.4.2                e1071_1.7-13                  spatstat.data_3.0-1           sp_1.6-0                     
[201] RColorBrewer_1.1-3            dqrng_0.3.0                   GenomeInfoDbData_1.2.9        progressr_0.13.0              memoise_2.0.1                
[206] knitr_1.43                    doParallel_1.0.17             vipor_0.4.5                   httpuv_1.6.8                  class_7.3-21                 
[211] irlba_2.3.5.1                 Rttf2pt1_1.3.12               BiocManager_1.30.21           classInt_0.4-8                seqLogo_1.62.0               
[216] pkgload_1.3.2                 jsonlite_1.8.4                Hmisc_4.7-2                   babelgene_22.9                digest_0.6.31                
[221] poweRlaw_0.70.6               rprojroot_2.0.3               bitops_1.0-7                  here_1.0.1                    RSQLite_2.2.20               
[226] globals_0.16.2                compiler_4.2.2                nnet_7.3-18                   statmod_1.5.0                 scran_1.24.1                 
[231] zoo_1.8-11                    pracma_2.4.2                  interactiveDisplayBase_1.34.0 gridGraphics_0.5-1            rlang_1.1.1                  
[236] urlchecker_1.0.1              nloptr_2.0.3                  uwot_0.1.14                   sessioninfo_1.2.2             rvest_1.0.3                  
[241] htmlwidgets_1.6.2             mvtnorm_1.1-3                 labeling_0.4.2                callr_3.7.3                   Cairo_1.6-0                  
[246] curl_5.0.0                    scater_1.24.0                 parallel_4.2.2                BiocNeighbors_1.16.0          edgeR_3.38.4                 
[251] filelock_1.0.2                scales_1.2.1                  RcppParallel_5.1.7            enrichplot_1.19.2             lme4_1.1-31                  
[256] deldir_1.0-6                  gridExtra_2.3                 bluster_1.6.0                 RCurl_1.98-1.10               GO.db_3.15.0                 
[261] MASS_7.3-58.2                 ellipsis_0.3.2                tidytree_0.4.2                spatstat.random_3.1-3         xml2_1.3.3                   
[266] rpart_4.1.19                  R6_2.5.1                      mclust_6.0.0                  GenomicAlignments_1.32.1      units_0.8-1  















