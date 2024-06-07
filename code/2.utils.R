#HGSC somatic landscape utilities file
#Author: Katherine Lawson-Michod

#Change directory
directory <- "~/Somatic Landscape Paper/hgsc_somatic_landscape_paper" #set working directory
library(tidyverse)
library(maftools)
library(dplyr)
library(ggplot2)
library(stringr)
library(caroline)
library(readr)

#########################################################DEFINITIONS#########################################################

#Define filters based on TCGA criteria for impact or consequence, these are types of variants that should be excluded in TCGA comparisons because they were exclude from TCGA analyses
variant_classification_filter <- c("Intron", "3'Flank", "3'UTR", "5'Flank", "5'UTR", "Splice_Region")

#Define filtered for false positives, these are genes that are notorious for false positives including olfactory receptors, list defined based on Hugo Symbol and RS number
#These orfs and rs will be filtered out when looking at the commonly mutated genes and the significantly mutated genes defined by mutsig 
bad_orfs_genes_1 <- c("AC008686.1", "AC118281.1", "ADGRL1", "AHNAK2", "ANKRD34A", "ANKRD36", "AP3D1", "BAGE2", "BIRC6", "CCDC187", "COBL", "COL3A1", 
                     "DNAH8", "EVL", "FAM90A22P", "FLG", "FRG2C", "FSIP2", "GJC2", "GOLGA6L10", "GOLGA6L2", "GOLGA6L20", "GOLGA6L9", "GSTT4", "GXYLT1",
                     "HLA-A", "HLA-B", "HLA-C", "HLA-DQA2", "HLA-DQB2", "HLA-DRB5", "HRNR", "IGFN1", "IGHV2-5", "IGHV4-31", "IGHV4-59", "IGHV4-61", 
                     "IGLV5-45", "KIAA0355", "KIAA1109", "KIF17", "KRBA1", "LRP1B", "LRRC37A3", "MADCAM1", "MAGEC1", "MRTFA", "MUC12", "MUC16", "MUC17", 
                     "MUC19", "MUC21", "MUC3A", "MUC4", "MUC5B", "MUC6", "NBPF9", "NPIPB15", "NPIPB6", "NXNL1", "OGFR", "OLFML2A", "OR10A2", "OR10A3", 
                     "OR10A4", "OR10A5", "OR10A6", "OR10A7", "OR10AA1P", "OR10AB1P", "OR10AC1", "OR10AD1", "OR10AE1P", "OR10AE3P", "OR10AF1P", "OR10AG1", 
                     "OR10AH1P", "OR10AK1P", "OR10B1P", "OR10C1", "OR10D1P", "OR10D3", "OR10D4P", "OR10D5P", "OR10G1P", "OR10G2", "OR10G3", "OR10G4", "OR10G5P", 
                     "OR10G6", "OR10G7", "OR10G8", "OR10G9", "OR10H1", "OR10H2", "OR10H3", "OR10H4", "OR10H5", "OR10J1", "OR10J2P", "OR10J3", "OR10J4", "OR10J5", 
                     "OR10J6P", "OR10J7P", "OR10J8P", "OR10J9P", "OR10K1", "OR10K2", "OR10N1P", "OR10P1", "OR10Q1", "OR10Q2P", "OR10R1P", "OR10R2", "OR10R3P", "OR10S1", 
                     "OR10T1P", "OR10T2", "OR10U1P", "OR10V1", "OR10V2P", "OR10V3P", "OR10V7P", "OR10W1", "OR10X1", "OR10Y1P", "OR10Z1", "OR11A1", "OR11G1P", "OR11G2", 
                     "OR11H1", "OR11H12", "OR11H13P", "OR11H2", "OR11H3P", "OR11H4", "OR11H5P", "OR11H6", "OR11H7", "OR11I1P", "OR11J1P", "OR11J2P", "OR11J5P", "OR11K1P")
bad_orfs_genes_2 <- c("OR11K2P", "OR11L1", "OR11M1P", "OR11N1P", "OR11P1P", "OR11Q1P", "OR12D1", "OR12D2", "OR12D3", "OR13A1", "OR13C1P", "OR13C2", "OR13C3",
                     "OR13C4", "OR13C5", "OR13C6P", "OR13C7", "OR13C8", "OR13C9", "OR13D1", "OR13D2P", "OR13D3P", "OR13E1P", "OR13F1", "OR13G1", "OR13H1", 
                     "OR13I1P", "OR13J1", "OR13K1P", "OR13Z1P", "OR13Z2P", "OR13Z3P", "OR14A16", "OR14A2", "OR14C36", "OR14I1", "OR14J1", "OR14K1", 
                     "OR14L1P", "OR1A1", "OR1A2", "OR1AA1P", "OR1AB1P", "OR1AC1P", "OR1B1", "OR1C1", "OR1D2", "OR1D3P", "OR1D4", "OR1D5", "OR1E1", "OR1E2", 
                     "OR1E3", "OR1F1", "OR1F12", "OR1F2P", "OR1G1", "OR1H1P", "OR1I1", "OR1J1", "OR1J2", "OR1J4", "OR1K1", "OR1L1", "OR1L3", "OR1L4", 
                     "OR1L6", "OR1L8", "OR1M1", "OR1M4P", "OR1N1", "OR1N2", "OR1P1", "OR1Q1", "OR1R1P", "OR1S1", "OR1S2", "OR1X1P", "OR1X5P", "OR2A1", 
                     "OR2A12", "OR2A13P", "OR2A14", "OR2A15P", "OR2A2", "OR2A20P", "OR2A25", "OR2A3P", "OR2A4", "OR2A41P", "OR2A42", "OR2A5", "OR2A7", 
                     "OR2A9P", "OR2AD1P", "OR2AE1", "OR2AF1P", "OR2AG1", "OR2AG2", "OR2AH1P", "OR2AI1P", "OR2AJ1", "OR2AK2", "OR2AL1P", "OR2AM1P", 
                     "OR2AO1P", "OR2AP1", "OR2AQ1P", "OR2AS1P", "OR2AS2P", "OR2AT1P", "OR2AT2P", "OR2AT4", "OR2B11", "OR2B2", "OR2B3", "OR2B4P", "OR2B6", 
                     "OR2B7P", "OR2B8P", "OR2BH1P", "OR2C1", "OR2C3", "OR2D2", "OR2D3", "OR2E1P", "OR2F1", "OR2F2", "OR2G1P", "OR2G2", "OR2G3", "OR2G6", 
                     "OR2H1", "OR2H2", "OR2H4P", "OR2H5P", "OR2I1P", "OR2J1", "OR2J2", "OR2J3", "OR2J4P", "OR2K2", "OR2L13", "OR2L1P", "OR2L2", "OR2L3", 
                     "OR2L5", "OR2L6P", "OR2L8", "OR2L9P", "OR2M1P", "OR2M2", "OR2M3", "OR2M4", "OR2M5", "OR2M7", "OR2N1P", "OR2P1P", "OR2Q1P", "OR2R1P")
bad_orfs_genes_3 <- c("OR2S1P", "OR2S2", "OR2T1", "OR2T10", "OR2T11", "OR2T12", "OR2T2", "OR2T27", "OR2T29", "OR2T3", "OR2T32P", "OR2T33", "OR2T34", 
                     "OR2T35", "OR2T4", "OR2T5", "OR2T6", "OR2T7", "OR2T8", "OR2U1P", "OR2U2P", "OR2V1", "OR2V2", "OR2W1", "OR2W2P", "OR2W3", "OR2W4P",
                     "OR2W5", "OR2W6P", "OR2X1P", "OR2Y1", "OR2Z1", "OR3A1", "OR3A2", "OR3A3", "OR3A4P", "OR3B1P", "OR3D1P", "OR4A10P", "OR4A11P", 
                     "OR4A12P", "OR4A13P", "OR4A14P", "OR4A15", "OR4A16", "OR4A17P", "OR4A18P", "OR4A19P", "OR4A1P", "OR4A21P", "OR4A2P", "OR4A3P", 
                     "OR4A40P", "OR4A41P", "OR4A42P", "OR4A43P", "OR4A44P", "OR4A45P", "OR4A46P", "OR4A47", "OR4A48P", "OR4A49P", "OR4A4P", "OR4A5", 
                     "OR4A50P", "OR4A6P", "OR4A7P", "OR4A8", "OR4A9P", "OR4B1", "OR4B2P", "OR4C10P", "OR4C11", "OR4C12", "OR4C13", "OR4C14P", "OR4C15",
                     "OR4C16", "OR4C1P", "OR4C2P", "OR4C3", "OR4C45", "OR4C46", "OR4C48P", "OR4C49P", "OR4C4P", "OR4C5", "OR4C50P", "OR4C6", "OR4C7P", 
                     "OR4C9P", "OR4D1", "OR4D10", "OR4D11", "OR4D12P", "OR4D2", "OR4D5", "OR4D6", "OR4D7P", "OR4D8P", "OR4D9", "OR4E1", "OR4E2", "OR4F13P", 
                     "OR4F14P", "OR4F15", "OR4F16", "OR4F17", "OR4F1P", "OR4F21", "OR4F28P", "OR4F29", "OR4F2P", "OR4F3", "OR4F4", "OR4F5", "OR4F6", 
                     "OR4F7P", "OR4F8P", "OR4G11P", "OR4G1P", "OR4G2P", "OR4G3P", "OR4G4P", "OR4G6P", "OR4H12P", "OR4H6P", "OR4K1", "OR4K11P", "OR4K12P", 
                     "OR4K13", "OR4K14", "OR4K15", "OR4K16P", "OR4K17", "OR4K2", "OR4K3", "OR4K4P", "OR4K5", "OR4K6P", "OR4K7P", "OR4K8P", "OR4L1", 
                     "OR4M1", "OR4M2", "OR4N1P", "OR4N2", "OR4N3P", "OR4N4", "OR4N5", "OR4P1P", "OR4P4", "OR4Q1P", "OR4Q2", "OR4Q3", "OR4R1P", "OR4R2P", 
                     "OR4R3P", "OR4S1", "OR4S2", "OR4T1P", "OR4U1P", "OR4V1P", "OR4W1P", "OR4X1", "OR4X2", "OR4X7P", "OR51A10P", "OR51A1P", "OR51A2", 
                     "OR51A3P", "OR51A4", "OR51A5P", "OR51A6P", "OR51A7", "OR51A8P", "OR51A9P", "OR51AB1P", "OR51B2", "OR51B3P", "OR51B4", "OR51B5",
                     "OR51B6", "OR51B8P", "OR51C1P", "OR51C4P", "OR51D1", "OR51E1", "OR51E2", "OR51F1", "OR51F2", "OR51F3P", "OR51F4P", "OR51F5P", 
                     "OR51G1", "OR51G2", "OR51H1", "OR51H2P", "OR51I1", "OR51I2", "OR51J1", "OR51K1P", "OR51L1", "OR51M1", "OR51N1P", "OR51P1P", "OR51Q1",
                     "OR51R1P", "OR51S1", "OR51T1", "OR51V1", "OR52A1", "OR52A4P", "OR52A5", "OR52B1P", "OR52B2", "OR52B3P", "OR52B4", "OR52B5P", "OR52B6", 
                     "OR52D1", "OR52E1", "OR52E2", "OR52E3P", "OR52E4", "OR52E5", "OR52E6", "OR52E7P", "OR52E8", "OR52H1", "OR52H2P", "OR52I1", "OR52I2", 
                     "OR52J1P", "OR52J2P", "OR52J3", "OR52K1", "OR52K2", "OR52K3P", "OR52L1", "OR52L2P", "OR52M1", "OR52M2P", "OR52N1", "OR52N2", "OR52N3P")
bad_orfs_genes_4 <- c("OR52N4", "OR52N5", "OR52P1P", "OR52P2P", "OR52Q1P", "OR52R1", "OR52S1P", "OR52T1P", "OR52U1P", "OR52V1P", "OR52W1", "OR52X1P", "OR52Y1P", "OR52Z1", "OR55B1P", "OR56A1", "OR56A3", "OR56A4", "OR56A5", "OR56A7P", "OR56B1", "OR56B2P", "OR56B3P", "OR56B4", "OR5A1", "OR5A2", "OR5AC1", "OR5AC2", "OR5AC4P", "OR5AH1P", "OR5AK1P", "OR5AK2", "OR5AK3P", "OR5AK4P", "OR5AL1", "OR5AL2P", "OR5AM1P", "OR5AN1", "OR5AN2P", "OR5AO1P", "OR5AP1P", "OR5AP2", "OR5AQ1P", "OR5AR1", "OR5AS1", "OR5AU1", "OR5AW1P", "OR5AZ1P", "OR5B10P", "OR5B12", "OR5B15P", "OR5B17", "OR5B19P", "OR5B1P", "OR5B2", "OR5B21", "OR5B3", "OR5BA1P", "OR5BB1P", "OR5BC1P", "OR5BD1P", "OR5BE1P", "OR5BH1P", "OR5BJ1P", "OR5BK1P", "OR5BL1P", "OR5BM1P", "OR5BN1P", "OR5BN2P", "OR5BP1P", "OR5BQ1P", "OR5BR1P", "OR5BS1P", "OR5BT1P", "OR5C1", "OR5D13", "OR5D14", "OR5D15P", "OR5D16", "OR5D17P", "OR5D18", "OR5D2P", "OR5D3P", "OR5E1P", "OR5F1", "OR5F2P", "OR5G1P", "OR5G3", "OR5G4P", "OR5G5P", "OR5H1", "OR5H14", "OR5H15", "OR5H2", "OR5H3P", "OR5H4P", "OR5H5P", "OR5H6", "OR5H7P", "OR5H8", "OR5I1", "OR5J1P", "OR5J2", "OR5J7P", "OR5K1", "OR5K2", "OR5K3", "OR5K4", "OR5L1", "OR5L2", "OR5M1", "OR5M10", "OR5M11", "OR5M12P","OR5M13P")
bad_orfs_genes_5 <- c("CGREF1", "TTN", "OR5M14P", "OR5M2P", "OR5M3", "OR5M4P", "OR5M5P", "OR5M6P", "OR5M7P", "OR5M8", "OR5M9", "OR5P1P", "OR5P2", "OR5P3", "OR5P4P", "OR5R1", "OR5S1P", "OR5T1", "OR5T2", "OR5T3", "OR5V1", "OR5W1P", "OR5W2", "OR6A2", "OR6B1", "OR6B2", "OR6B3", "OR6C1", "OR6C2", "OR6C3", "OR6C4", "OR6C5P", "OR6C6", "OR6C64P", "OR6C65", "OR6C66P", "OR6C68", "OR6C69P", "OR6C70", "OR6C71P", "OR6C72P", "OR6C73P", "OR6C74", "OR6C75", "OR6C76", "OR6C7P", "OR6D1P", "OR6E1P", "OR6F1", "OR6J1", "OR6K1P", "OR6K2", "OR6K3", "OR6K4P", "OR6K5P", "OR6K6", "OR6L1P", "OR6L2P", "OR6M1", "OR6M2P", "OR6M3P", "OR6N1", "OR6N2", "OR6P1", "OR6Q1", "OR6R1P", "OR6R2P", "OR6S1", "OR6T1", "OR6U2P", "OR6V1", "OR6W1P", "OR6X1", "OR6Y1", "OR7A10", "OR7A11P", "OR7A15P", "OR7A17", "OR7A18P", "OR7A19P", "OR7A1P", "OR7A2P", "OR7A3P", "OR7A5", "OR7A8P", "OR7C1", "OR7C2", "OR7D11P", "OR7D1P", "OR7D2", "OR7D4", "OR7E100P", "OR7E101P", "OR7E102P", "OR7E104P", "OR7E105P", "OR7E106P", "OR7E108P", "OR7E109P", "OR7E10P", "OR7E110P", "OR7E111P", "OR7E115P", "OR7E116P", "OR7E117P", "OR7E11P", "OR7E121P", "OR7E122P", "OR7E125P", "OR7E126P", "OR7E128P", "OR7E129P", "OR7E12P", "OR7E130P", "OR7E136P", "OR7E13P", "OR7E140P", "OR7E145P", "OR7E148P", "OR7E149P", "OR7E14P", "OR7E154P", "OR7E155P", "OR7E156P", "OR7E157P", "OR7E158P", "OR7E159P", "OR7E15P", "OR7E160P", "OR7E161P", "OR7E162P", "OR7E163P", "OR7E16P", "OR7E18P", "OR7E19P", "OR7E1P", "OR7E21P", "OR7E22P", "OR7E23P", "OR7E24", "OR7E25P", "OR7E26P", "OR7E28P", "OR7E29P", "OR7E2P", "OR7E31P", "OR7E33P", "OR7E35P", "OR7E36P", "OR7E37P", "OR7E38P", "OR7E39P", "OR7E41P", "OR7E43P", "OR7E46P", "OR7E47P", "OR7E4P", "OR7E53P", "OR7E55P", "OR7E59P", "OR7E5P", "OR7E62P", "OR7E66P", "OR7E7P", "OR7E83P", "OR7E84P", "OR7E85P", "OR7E86P", "OR7E87P", "OR7E89P", "OR7E8P", "OR7E90P", "OR7E91P", "OR7E93P", "OR7E94P", "OR7E96P", "OR7E97P", "OR7E99P", "OR7G1", "OR7G15P", "OR7G2", "OR7G3", "OR7H1P", "OR7H2P", "OR7K1P", "OR7L1P", "OR7M1P", "OR8A1", "OR8A2P", "OR8A3P", "OR8B10P", "OR8B12", "OR8B1P", "OR8B2", "OR8B3", "OR8B4", "OR8B5P", "OR8B6P", "OR8B7P", "OR8B8", "OR8B9P", "OR8C1P", "OR8D1", "OR8D2", "OR8D4", "OR8F1P", "OR8G1", "OR8G2P", "OR8G3P", "OR8G5", "OR8G7P", "OR8H1", "OR8H2", "OR8H3", "OR8I1P", "OR8I2", "OR8I4P", "OR8J1", "OR8J2", "OR8J3", "OR8K1", "OR8K2P", "OR8K3", "OR8K4P", "OR8K5", "OR8L1P", "OR8Q1P", "OR8R1P", "OR8S1", "OR8S21P", "OR8T1P", "OR8U1", "OR8U8", "OR8U9", "OR8V1P", "OR8X1P", "OR9A1P", "OR9A2", "OR9A3P", "OR9A4", "OR9G1", "OR9G2P", "OR9G3P", "OR9G4", "OR9G9", "OR9H1P", "OR9I1", "OR9I2P", "OR9I3P", "OR9K1P", "OR9K2", "OR9L1P", "OR9M1P", "OR9N1P", "OR9P1P", "OR9Q1", "OR9Q2", "OR9R1P", "OR9S24P", "PABPC3", "PHLDB3", "PKHD1L1", "PLIN4", "PRAMEF33", "PRSS1", "RAPGEF2", "RBM33", "RNF44", "RP1L1", "RTL1", "SCN7A", "SLC25A5", "SLMAP", "TAL1", "TAS2R19", "TAS2R30", "TAS2R31", "TAS2R43", "TFAP4", "TOPAZ1", "TRBV3-1", "TRBV4-2", "VPS13A", "VPS13D", "ZMYM3", "ZNF575", "ZXDA")
bad_orfs_genes <- c(bad_orfs_genes_1, bad_orfs_genes_2, bad_orfs_genes_3, bad_orfs_genes_4, bad_orfs_genes_5)
rm(bad_orfs_genes_1, bad_orfs_genes_2, bad_orfs_genes_3, bad_orfs_genes_4, bad_orfs_genes_5)
bad_orfs_rs <- c("rs1010423", "rs1010424", "rs1031321426", "rs1048801", "rs10536197", "rs10907643", "rs11075798", "rs111329392", "rs112207161", "rs112372129", "rs11261835", "rs1126537", "rs1127430", "rs11276076",             "rs113022012", "rs113027958", "rs1130466", "rs113202486", "rs11322783", "rs113264216", "rs113322110", "rs1135840", "rs11537002", "rs11548735", "rs1158446987", "rs117387753", "rs11754361", "rs11756038",
                 "rs11756039", "rs1178200495", "rs11782383", "rs1178430419", "rs1182056382", "rs1187987064", "rs1188045977", "rs1189430066", "rs1196807095", "rs1198056423", "rs1212108435", "rs12199346", "rs1240508430",
                 "rs12466789", "rs1253314241", "rs12627379", "rs12628603", "rs1263810", "rs1281607798", "rs1283383837", "rs12895357", "rs1305272775", "rs13054858", "rs1311546360", "rs1323591162", "rs1366091859",
                 "rs1366527604", "rs1375638824", "rs1384446739", "rs139106189", "rs1411096827", "rs142024473", "rs1422380205", "rs1424766495", "rs144160315", "rs144403657", "rs1451772", "rs145554461", "rs1469302458",
                 "rs1469676804", "rs1472236985", "rs1479173405", "rs148013251", "rs1491233894", "rs149795204", "rs149815285", "rs150690007", "rs1553961403", "rs1554844974", "rs1555889552", "rs1557586047", "rs1628172",
                 "rs1650017", "rs1735169", "rs1736772", "rs17667531", "rs1800734", "rs1834697", "rs1993829", "rs199522572", "rs199641655", "rs199726057", "rs199736618", "rs199938577", "rs199941915", "rs199993063",
                 "rs199995136", "rs200027230", "rs200032924", "rs200044038", "rs200075071", "rs200110372", "rs200207579", "rs200253954", "rs2003432", "rs200604108", "rs200684350", "rs200862596", "rs200981173",
                 "rs2012761", "rs201411101", "rs201518955", "rs201754859", "rs201899861", "rs201988838", "rs202003805", "rs202149101", "rs2036553", "rs2072032", "rs2073525", "rs2075310", "rs2225470", "rs2236293",
                 "rs2240906", "rs224331", "rs2285888", "rs2338626", "rs2338627", "rs2340550", "rs2550232", "rs2669761", "rs2743967", "rs2770151", "rs2798901", "rs28376898", "rs28559695", "rs2943510", "rs3006286",
                 "rs3009004", "rs3060668", "rs3077646", "rs3085220", "rs3106194", "rs3124747", "rs3218384", "rs33917318", "rs34200684", "rs34628045", "rs35295676", "rs35896628", "rs35993975", "rs36219868", "rs36657",
                 "rs368768685", "rs369823368", "rs370384298", "rs373032", "rs3753270", "rs3762739", "rs377599213", "rs3803738", "rs3810366", "rs3814960", "rs3825942", "rs3828410", "rs3830472", "rs3830809", "rs383369",
                 "rs3837575", "rs4575819", "rs4782272", "rs4879782", "rs5027409", "rs512770", "rs528635259", "rs529014403", "rs532845932", "rs533172496", "rs539345546", "rs546394862", "rs549446", "rs555930275",
                 "rs557008123", "rs55745992", "rs561234742", "rs56147919", "rs565068798", "rs575085", "rs5845253", "rs58500707", "rs60308484", "rs61430934", "rs62089225", "rs62399429", "rs62536540", "rs642223",
                 "rs6599528", "rs66628201", "rs66721522", "rs67495004", "rs6782766", "rs688976", "rs71082910", "rs71119069", "rs71265055", "rs7148408", "rs72344675", "rs72486387", "rs72487164", "rs7255721", "rs72847018",
                 "rs73011014", "rs73069586", "rs73069587", "rs73112044", "rs73112046", "rs731170", "rs73126218", "rs73212641", "rs73269903", "rs73465610", "rs73920826", "rs745614033", "rs745817378", "rs745873123",
                 "rs746307676", "rs746804941", "rs747141768", "rs747391", "rs748260316", "rs749969491", "rs750338758", "rs751112872", "rs751186374", "rs7513079", "rs75275897", "rs753316486", "rs753994746", "rs75444444",
                 "rs755706336", "rs756349688", "rs756366019", "rs75659869", "rs757300194", "rs757723965", "rs758599717", "rs76014288", "rs76065392", "rs76100089", "rs76291283", "rs763551272", "rs763997493", "rs7640178",
                 "rs765451626", "rs765819981", "rs76666113", "rs767064362", "rs767462612", "rs768409846", "rs77324563", "rs773963042", "rs77453597", "rs775584235", "rs776107371", "rs776456215", "rs77797695",
                 "rs778823914", "rs778867790", "rs779610049", "rs780538719", "rs78071782", "rs781577243", "rs781699454", "rs78214734", "rs78553419", "rs78720771", "rs7900838", "rs79202331", "rs79489020", "rs8017603",
                 "rs80207539", "rs80323191", "rs8074547", "rs865863555", "rs866674389", "rs868981506", "rs878856926", "rs878862357", "rs878943546", "rs879053212", "rs879095499", "rs879235545", "rs879250200", "rs9328936",
                 "rs9422022", "rs9603837", "rs9605146", "rs9610841", "rs979369776", "rs144848")

#Define variable classifications for MAF file set up for tissue grant somatic variants
maf_factor_vars <- c('Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'Chromosome', 'Strand', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'Exon_Number', 'all_effects', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'STRAND_VEP', 'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE', 'CANONICAL', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'RefSeq', 'SIFT', 'PolyPhen', 'INTRON', 'DOMAINS', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'IMPACT', 'VARIANT_CLASS', 'HGVS_OFFSET', 'PHENO', 'GENE_PHENO', 'Variant_Classification', 'FILTER', 'flanking_bps', 'vcf_id', 'vcf_pos', 'NCBI_Build', 'Variant_Type', 'Tumor_Sample_Barcode', 'Sequence_Source', 'Sequencer', 'Start_Position', 'End_Position', 'TSL')
maf_numeric_vars <- c('vcf_qual', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF', 'T_DP', 'T_AF', 'N_DP', 'N_AF', 'BKZ', 'DISTANCE', 'AF', 'AFR_AF', 'AMR_AF', 'ASN_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF')
remove_columns <- c('BAM_File', 'dbSNP_Val_Status', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Matched_Norm_Sample_UUID', 'Mutation_Status', 'Score', 'Sequencing_Phase', 'Tumor_Sample_UUID', 'Tumor_Validation_Allele1', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'Tumor_Validation_Allele2', 'Verification_Status', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_ref_count', 'n_alt_count', 'ALLELE_NUM', 'EXON', 'MINIMISED', 'PICK')

#Define tcga major gene list and tissue grant top 10 significantly mutated genes list
tcga_genes <- c("TP53", "CSMD3", "FAT3", "NF1", "BRCA2", "CDK12", "BRCA1", "GABRA6", "RB1") #Major genes identified by TCGA
tissue_top10_filtered <- c("TP53", "CGREF1", "CLEC18B", "RRP7A", "NPC1L1", "KRAS", "ALG1", "PABPC1", "BRCA1", "PDZK1", "RB1")

#Define genes to exclude (either false positive in manual review, or not covered by TCGA)
genes_exclude <- c('ALG1', 'LILRA5', 'HLA-DRB1', 'LILRA1', 'KRT18', 'PRSS3', 'ZNF714', 'CNTNAP5', 'HUWE1', 'NEB', 'PAK3', 'SPDYE5', 'ZNF729')

#Read in dictionary for annotating common snps with ALFA frequency - this dictionary contains dbSNP_RS numbers and AF reported in African and overall populations by ALFA
dbSNP_alfa_dict <- readxl::read_xlsx(file.path(directory, "reference_data/dbSNP_alfa_dict.xlsx"))

#Read in list of SUIDs that failed sample qc 
tissue_fail_suids <- read.delim(file.path(directory, "reference_data/fail_qc_suids_somatic_landscape.txt")) %>% pull(tissue_fail_suids) %>% as.character()
tissue_pass_suids <- read.delim(file.path(directory, "reference_data/pass_qc_suids_somatic_landscape.txt")) %>% pull(tissue_pass_suids) %>% as.character()
length(tissue_pass_suids)

all_suids <- append(tissue_fail_suids, tissue_pass_suids)
df <- as.data.frame(all_suids)

#########################################################SCHILDKRAUT-B#########################################################

#Read in somatic mutations and clinical data 
hgsc_38_high_confidence <- read.delim(file.path(directory, "data/somatic/schildkraut/hgsc_38_high_confidence_processed.maf"), header=TRUE, comment.char="#") %>% filter(SUID %in% tissue_pass_suids)
tissue_epi_cosmic <- readxl::read_xlsx(file.path(directory, "data/epi/schildkraut/cosmic_epidata_10112023.xlsx")) %>% filter(., suid %in% tissue_pass_suids)
stage_1_suids <- tissue_epi_cosmic %>% filter(stage == 1) %>% pull(suid)
length(unique(hgsc_38_high_confidence$SUID)) #20 Stage I suids
hgsc_38_high_confidence_exclude_stage_1 <- hgsc_38_high_confidence %>% filter(!SUID %in% stage_1_suids) #exclude stage I SUIDs
length(unique(hgsc_38_high_confidence_exclude_stage_1$SUID))

tissue_epi_cosmic <- readxl::read_xlsx(file.path(directory, "data/epi/schildkraut/cosmic_epidata_10112023.xlsx"))
summary(as.factor(tissue_epi_cosmic$site))
#Define SUIDs for cases that have high-grade endometriod histotype and are missing TP53 mutation and exclude stage I cases
endo_case_ids <- tissue_epi_cosmic %>% filter(., histology == '2') %>% pull(suid) %>% unique()
tissue_black_endo_somatic_mutations <- hgsc_38_high_confidence %>% filter(SUID %in% endo_case_ids)
endo_with_tp53_ids <- tissue_black_endo_somatic_mutations %>% filter(Hugo_Symbol == "TP53") %>% pull(SUID) %>% unique()
endo_without_tp53_ids <- setdiff(tissue_black_endo_somatic_mutations$SUID, endo_with_tp53_ids) #3 cases has high-grade endometriod histology and TP53 mutations
hgsc_38_high_confidence_exclude_stage_1_and_endo <- hgsc_38_high_confidence_exclude_stage_1 %>% filter(!SUID %in% endo_without_tp53_ids) #Remove endo cases without TP53 mutations from epi data 
tissue_black_suids_final <- unique(hgsc_38_high_confidence_exclude_stage_1_and_endo$SUID) 
length(unique(hgsc_38_high_confidence_exclude_stage_1_and_endo$SUID)) #final analytic cases 191
tissue_pass_suids_df <- as.data.frame(tissue_pass_suids)
tissue_pass_suids_exclude_endo <- tissue_pass_suids_df %>% filter(!tissue_pass_suids %in% endo_without_tp53_ids)
writexl::write_xlsx(tissue_pass_suids_exclude_endo, file.path(directory, "reference_data/pass_qc_suids_not_endo.xlsx"))
tissue_pass_suids_exclude_endo_stage1 <- tissue_pass_suids_exclude_endo %>% filter(tissue_pass_suids %in% tissue_black_suids_final)
writexl::write_xlsx(tissue_pass_suids_exclude_endo_stage1, file.path(directory, "reference_data/pass_qc_suids_not_endo_stage_1.xlsx"))

# Filter somatic mutation data to include only the final 191 analytic samples and format data frame
tissue_black_somatic_mutations <- hgsc_38_high_confidence_exclude_stage_1_and_endo %>% 
  filter(SUID %in% tissue_black_suids_final & !Hugo_Symbol %in% genes_exclude) %>%
  mutate_at(maf_factor_vars, factor) %>%
  mutate_at(maf_numeric_vars, as.numeric) %>%
  mutate(dbSNP_RS = str_replace(dbSNP_RS, "rs", ""))
tissue_black_epi_final <- tissue_epi_cosmic %>% filter(suid %in% tissue_black_suids_final)

# Filter somatic mutation data to exclude common SNPs defined by ALFA and gnomAD in addition to bad orfs list from Aaron
tissue_black_somatic_mutations <- merge(tissue_black_somatic_mutations, dbSNP_alfa_dict, all.x = T)
tissue_AFR_AF_Fail <- tissue_black_somatic_mutations %>% filter(gnomAD_AFR_AF >= 0.01 | AFR_AF >= 0.01 | AF_alt_ALFA_african >= 0.01)
tissue_black_somatic_mutations_filter1 <- anti_join(tissue_black_somatic_mutations, tissue_AFR_AF_Fail)
tissue_black_snp_filtered_somatic_mutations <- tissue_black_somatic_mutations_filter1 %>% filter(!Hugo_Symbol %in% bad_orfs_genes)
summary(tissue_black_snp_filtered_somatic_mutations$Variant_Classification)
tissue_black_tcga_filtered_somatic_mutations <- tissue_black_snp_filtered_somatic_mutations %>% filter(!Variant_Classification %in% variant_classification_filter) #Apply TCGA criteria for variant classification (ie consequence) - basically filtering all noncoding 
tissue_black_somatic_mutations <- tissue_black_tcga_filtered_somatic_mutations
length(unique(tissue_black_somatic_mutations$SUID))
write.delim(tissue_black_somatic_mutations, file.path(directory, "data/somatic/schildkraut/hgsc_38_high_confidence_processed_filtered.maf"))
tissue_black_somatic_tcga_mutations_maf <- maftools::read.maf(file.path(directory, "data/somatic/schildkraut/hgsc_38_high_confidence_processed_filtered.maf"))

#########################################################TCGA#########################################################

#read in TCGA somatic variant file and recode sample IDs to match clinical data file
tcga_somatic_mutations <- read.delim(file.path(directory, "data/somatic/tcga/TCGA_somatic_mutations.maf"))
tcga_somatic_mutations <- tcga_somatic_mutations %>% mutate(case_submitter_id = str_remove(Tumor_Sample_Barcode, "-01.-.*"))
tcga_wes_cases <- c(unique(tcga_somatic_mutations$case_submitter_id)) #define list of tcga WES ids that will be used to subset clinical data to cases that have tumor WES

#read in clinical TCGA data, remove duplicates, and filter to cases with tumor WES 
tcga_clinical <- read_tsv(file.path(directory, "data/epi/tcga/clinical.tsv"))
tcga_clinical$duplicate <- as.factor(duplicated(tcga_clinical$case_submitter_id)) #define indicator for dropping duplicates
tcga_clinical_clean <- tcga_clinical %>% filter(duplicate == TRUE) #drop duplicates
tcga_clinical_wes <- tcga_clinical_clean %>% filter(., case_submitter_id %in% tcga_wes_cases) #316 observations

#subset to white cases and cases younger than 79 to match inclusion for Schildkraut-B
tcga_white_cases <- tcga_clinical_wes %>% filter(., race == "white")
tcga_white_cases_ids <- c(unique(tcga_white_cases$case_submitter_id))
tcga_white_cases <- tcga_white_cases %>% 
  mutate(age_at_diagnosis = recode(age_at_diagnosis, "'--" = NA_character_)) %>% 
  mutate_at(vars('age_at_diagnosis'), as.numeric) %>% 
  mutate('Age at diagnosis' = round(age_at_diagnosis/365, 0)) 
tcga_white_cases_under_79 <- tcga_white_cases %>% filter(., `Age at diagnosis` < 79 | is.na(`Age at diagnosis`))

tcga_white_cases_under_79 %>% group_by(ethnicity) %>% tally()
tcga_white_cases_under_79_ids <- c(unique(tcga_white_cases_under_79$case_submitter_id))
tcga_white_somatic_mutations <- tcga_somatic_mutations %>% filter(., case_submitter_id %in% tcga_white_cases_under_79_ids) #filter somatic data to only include individuals under 79
tcga_white_somatic_mutations <- tcga_white_somatic_mutations %>% mutate_all(., as.factor)
tcga_white_somatic_mutations <- tcga_white_somatic_mutations %>% mutate(dbSNP_RS = str_replace(dbSNP_RS, "rs", ""))
tcga_white_somatic_mutations <- merge(tcga_white_somatic_mutations, dbSNP_alfa_dict, all.x = T)
tcga_white_somatic_mutations_filtered <- tcga_white_somatic_mutations %>% filter(!AF_alt_ALFA_total >= 0.01 | is.na(AF_alt_ALFA_total)) #filter common SNPs
tcga_white_somatic_mutations_filtered <- tcga_white_somatic_mutations %>% filter(., !Hugo_Symbol %in% bad_orfs_genes)
tcga_white_somatic_mutations_filtered <- tcga_white_somatic_mutations_filtered %>% filter(., !Hugo_Symbol %in% genes_exclude)
tcga_white_somatic_mutations <- tcga_white_somatic_mutations_filtered %>% filter(!Variant_Classification %in% variant_classification_filter) #filter variant classifications that were not concordant across TCGA and Schildkraut
write.delim(tcga_white_somatic_mutations, file.path(directory, "data/somatic/tcga/tcga_white_somatic_mutations_filtered.maf"))
tcga_white_somatic_mutations_maf <- maftools::read.maf(file.path(directory, "data/somatic/tcga/tcga_white_somatic_mutations_filtered.maf"))

rm(tcga_clinical, tcga_clinical_clean, tcga_clinical_wes, tcga_somatic_mutations, tcga_wes_cases)
##########################################################MUTSIG2CV#########################################################

cosmic_known <- c("Age_B", "MMRD_B", "POLE_B", "HRD_B", "BER_B", "Chemotherapy_B", "Immunosuppressants_B", "APOBEC_B", "Tobacco_B", "UV_B", "Artifact_B", "Lymphoid_B")

#Read in final MutSig2CV results and counts file
mutsig2cv_counts <- readxl::read_xlsx(file.path(directory, "mutsig2cv/results/mutsig2cv_results_snp_filtered_FINAL.xlsx"))
mutsig2cv_top20 <- head(mutsig2cv_counts, 10)
mutsig2cv_top <- mutsig2cv_top20 %>% filter(!gene %in% genes_exclude)
tissue_genes <- unique(mutsig2cv_top$gene)

MMR <- hgsc_38_high_confidence %>% filter(Hugo_Symbol %in% c("MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2")) %>% group_by(CLIN_SIG) %>% tally()
hgsc_38_high_confidence_exclude_stage_1_and_endo %>% filter(Hugo_Symbol %in% c("MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2") & CLIN_SIG == "uncertain_significance")

tissue_epi_cosmic_unknown <- tissue_epi_cosmic %>% filter(is.na(neoadj_treat)) %>% select(suid, ocwaaid, stage, neoadj_treat)
writexl::write_xlsx(tissue_epi_cosmic_unknown, "/Users/kayleighlawson-michod/Library/CloudStorage/OneDrive-UniversityofUtah/Somatic Landscape Paper/hgsc_somatic_landscape_paper/neoadj_unknown_suids.xlsx")


### Look at TMB by the MMRD signature yes/no
all_results <- readxl::read_xlsx(file.path(directory, "results/all_results.xlsx"))
all_results %>% filter(Study == "Schildkraut-B") %>% pull(MMRD) %>% count()
mmrd_suid_schildkrautB <- all_results %>% filter(MMRD == "1" & Study == "Schildkraut-B") %>% pull(sample_id)
mmrd_suid_tcgaW <- all_results %>% filter(MMRD == "1" & Study == "TCGA-W") %>% pull(sample_id)

schildkrautB_somatic_mmrd <- tissue_black_somatic_mutations %>% filter(SUID %in% mmrd_suid_schildkrautB)
schildkrautB_somatic_NOmmrd <- tissue_black_somatic_mutations %>% filter(!SUID %in% mmrd_suid_schildkrautB)

write.delim(schildkrautB_somatic_mmrd, file.path(directory, "data/somatic/schildkraut/tissue_black_mmrd.maf"))
write.delim(schildkrautB_somatic_NOmmrd, file.path(directory, "data/somatic/schildkraut/tissue_black_NOmmrd.maf"))


