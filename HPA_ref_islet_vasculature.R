# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 05/23/2022
# R version 3.6.2 (2019-12-12) 'Dark and Stormy Night'
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

















