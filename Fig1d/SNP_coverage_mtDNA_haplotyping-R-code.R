setwd("C://Bioinformatics/Heck_USDA_Pathology/Dcitri_mtDNA_haplotyping/Examining_the_effect_of_kmer_size/ska2/ska_mapping/Falk_dataset")
library(BioCircos)
library(RColorBrewer)
snippy <- read.table("Falk-k15-m90-mapped-OM181945.1-snps-tab-for-R-input-snp-proportions.txt", sep = "\t", stringsAsFactors = F, na.strings = "NA")

colnames(snippy) <- c("Strain_ID", "Genus", "Chrom", "Pos", "Type", "Ref", "Alt", "Evidence", "FType", "Strand", "NT_pos", "AA_pos", "Effect", "Locus_tag", "Gene", "Product")

snippy[is.na(snippy)] <- ''

sniptest <- snippy
newprotsamps <- c("KY426014","KY426015","MF181946","MF181947","metasample","MF614804","MF614805","MF614806","MF614807","MF614808","MF614809","MF614810","MF614811","MF614812","MF614813","MF614814","MF614815","MF614816","MF614817","MF614818","MF614819","MF614820","MF614821","MF614822","MF614823","MF614824","MF614825","MF614826","MF614827","MF614828", "OM181946", "MF181947")
sniptest$protocol <- ifelse(sniptest$Strain_ID %in% newprotsamps, 'new_protocol', 'existing_protocol')

# use a high coverage CLas strain KY426014 as the example to plot

sniplen <- length(sniptest[which(sniptest$Genus == "Dcitri" & sniptest$Strain_ID == "metasample"),]$Pos)

# create a new column in sniptest to hold the read depth reported in the "Evidence"
# column by snippy
Evidence_depth <- c()
for(i in sniptest$Evidence){
  Evidence_depth <- c(Evidence_depth, as.numeric(strsplit(strsplit(strsplit(strsplit(i, " ")[[1]], " ")[[1]],":")[[1]], " ")[[2]]))
  print(as.numeric(strsplit(strsplit(strsplit(strsplit(i, " ")[[1]], " ")[[1]],":")[[1]], " ")[[2]]))
}

sniptest$SNP_depth <- Evidence_depth
# create a variable to change shape of SNP value if SNP inside a coding sequence
# or intergenic region
#sniptest$shape <- ifelse(sniptest$FType == "CDS", 'circle', 'rect')

# BioCircos plots

# genome track
Dcit <- list('OM181945.1' = 15039)

Dcit_mtDNA_CDS_arc_labels_weitz = c('ND3', 'COX3', 'ATP6', 'ATP8', 'COX2', 'COX1', 'ND2', 'ND1', 'CYTB', 'ND6', 'ND4L', 'ND4', 'ND5')
Dcit_mtDNA_CDS_arc_labels = c('ND3', 'COX3', 'ATP6', 'ATP8', 'COX2', 'COX1', 'ND2', 'ND1', 'CYTB', 'ND6', 'ND4L', 'ND4', 'ND5')
Dcit_mtDNA_CDS_arc_chroms = c('OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1')
Dcit_mtDNA_CDS_arc_start <- c(352, 776, 1558, 2226, 2508, 3236, 4951, 9091, 10085, 11221, 11847, 12116, 13420)
Dcit_mtDNA_CDS_arc_stop <- c(702, 1558, 2232, 2378, 3171, 4765, 5922, 10002, 11227, 11703, 12122, 13360, 15037)

Dcit_mtDNA_tRNA_arc_labels_weitz = c('tRNA-Phe', 'tRNA-Glu', 'tRNA-Ser', 'tRNA-Asn', 'tRNA-Arg', 'tRNA-Ala', 'tRNA-Gly', 'tRNA-Asp', 'tRNA-Lys', 'tRNA-Leu', 'tRNA-Tyr', 'tRNA-Cys', 'tRNA-Trp', 'tRNA-Met', 'tRNA-Gln', 'tRNA-Ile', 'tRNA-Val', 'tRNA-Leu', 'tRNA-Ser', 'tRNA-Pro', 'tRNA-Thr', 'tRNA-His')
Dcit_mtDNA_tRNA_arc_labels = c('tRNA-Phe', 'tRNA-Glu', 'tRNA-Ser', 'tRNA-Asn', 'tRNA-Arg', 'tRNA-Ala', 'tRNA-Gly', 'tRNA-Asp', 'tRNA-Lys', 'tRNA-Leu', 'tRNA-Tyr', 'tRNA-Cys', 'tRNA-Trp', 'tRNA-Met', 'tRNA-Gln', 'tRNA-Ile', 'tRNA-Val', 'tRNA-Leu', 'tRNA-Ser', 'tRNA-Pro', 'tRNA-Thr', 'tRNA-His')
Dcit_mtDNA_tRNA_arc_chroms = c('OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1')
Dcit_mtDNA_tRNA_arc_start <- c(3, 50, 114, 168, 233, 294, 703, 2379, 2438, 3172, 4769, 4831, 4892, 5923, 5988, 6061, 7827, 9027, 10024, 11706, 11772, 13360)
Dcit_mtDNA_tRNA_arc_stop <- c(58, 112, 167, 231, 293, 353, 761, 2439, 2507, 3232, 4830, 4890, 4952, 5987, 6054, 6123, 7890, 9090, 10086, 11771, 11832, 13419)

Dcit_mtDNA_rRNA_arc_labels_weitz = c('12S_rRNA', '16S_rRNA')
Dcit_mtDNA_rRNA_arc_labels = c('12S_rRNA', '16S_rRNA')
Dcit_mtDNA_rRNA_arc_chroms = c('OM181945.1', 'OM181945.1')
Dcit_mtDNA_rRNA_arc_start <- c(7069, 7891)
Dcit_mtDNA_rRNA_arc_stop <- c(7826, 9026)

Dcit_mtDNA_intergenic_arc_labels_weitz = c('int1', 'int2', 'int3', 'int4', 'int5', 'int6', 'int7', 'int8', 'int9', 'int10', 'int11', 'int12')
Dcit_mtDNA_intergenic_arc_labels = c('int1', 'int2', 'int3', 'int4', 'int5', 'int6', 'int7', 'int8', 'int9', 'int10', 'int11', 'int12')
Dcit_mtDNA_intergenic_arc_chroms = c('OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1', 'OM181945.1')
Dcit_mtDNA_intergenic_arc_start <- c(113, 232, 762, 3233, 4766, 4891, 6055, 6124, 10003, 11704, 11833, 15038)
Dcit_mtDNA_intergenic_arc_stop <- c(113, 232, 775, 3235, 4768, 4891, 6060, 7068, 10023, 11705, 11846, 15039)

CDStrack = BioCircosArcTrack('Dcit_CDS', Dcit_mtDNA_CDS_arc_chroms, Dcit_mtDNA_CDS_arc_start, Dcit_mtDNA_CDS_arc_stop,
                               minRadius = 0.45, maxRadius = 1, opacities = 0.25, labels = Dcit_mtDNA_CDS_arc_labels, colors = c("#6DCBBA", "#6DCBBA", "#6DCBBA", "#6DCBBA", "#6DCBBA", "#6DCBBA", "#6DCBBA", "#6DCBBA", "#6DCBBA", "#6DCBBA", "#6DCBBA", "#35B5A5", "#35B5A5"))
tRNAtrack = BioCircosArcTrack('Dcit_tRNA', Dcit_mtDNA_tRNA_arc_chroms, Dcit_mtDNA_tRNA_arc_start, Dcit_mtDNA_tRNA_arc_stop,
                             minRadius = 0.45, maxRadius = 1, opacities = 0.8, labels = Dcit_mtDNA_tRNA_arc_labels, colors = "#FFFF92")
rRNAtrack = BioCircosArcTrack('Dcit_rRNA', Dcit_mtDNA_rRNA_arc_chroms, Dcit_mtDNA_rRNA_arc_start, Dcit_mtDNA_rRNA_arc_stop,
                             minRadius = 0.45, maxRadius = 1, opacities = 0.8, labels = Dcit_mtDNA_rRNA_arc_labels, colors = "#FD6C5C")
Intergenic_track = BioCircosArcTrack('Dcit_int', Dcit_mtDNA_intergenic_arc_chroms, Dcit_mtDNA_intergenic_arc_start, Dcit_mtDNA_intergenic_arc_stop,
                              minRadius = 0.45, maxRadius = 1, opacities = 0.8, labels = Dcit_mtDNA_intergenic_arc_labels, colors = "#B3ACDC")

#tracklist = CDStrack

# create a link track to indicate gene names of SNPs
link1 <- rep('OM181945.1', sniplen)
link2 <- rep('OM181945.1', sniplen)
linkpos1 <- sniptest[which(sniptest$Genus == "Dcitri" & sniptest$Strain_ID == "metasample"),]$Pos
linkpos2 <- sniptest[which(sniptest$Genus == "Dcitri" & sniptest$Strain_ID == "metasample"),]$Pos + 1
linklabs <- sniptest[which(sniptest$Genus == "Dcitri" & sniptest$Strain_ID == "metasample"),]$Gene

trackback = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 0.45,
                                     borderSize = 0, fillColors = "#EEFFEE") 


trackback2 = BioCircosBackgroundTrack("mySNPBackgroundTrack", minRadius = 0.45, maxRadius = 1,
                                      borderSize = 0.1, borderColors = "black", fillColors = "white") 

snpColors <- brewer.pal(length(levels(as.factor(sniptest[which(sniptest$Genus == "Dcitri" & sniptest$Strain_ID == "metasample"),]$Type))),"Dark2")
names(snpColors) <- levels(as.factor(sniptest[which(sniptest$Genus == "Dcitri" & sniptest$Strain_ID == "metasample"),]$Type))
snpColors2 <- snpColors[sniptest[which(sniptest$Genus == "Dcitri" & sniptest$Strain_ID == "metasample"),]$Type]

tracklist = BioCircosLinkTrack('snps', link1, link2,
                                           linkpos1, link2, linkpos1, linkpos2,
                                           maxRadius = 0.3, #labels = linklabs,
                                           axisColor = 'white', labelColor = as.vector(snpColors2))

snplist = trackback2 + BioCircosSNPTrack('snptrack', link1, linkpos1, 
                                         sniptest[which(sniptest$Genus == "Dcitri" & sniptest$Strain_ID == "metasample"),]$SNP_depth, 
                                         colors = "gray20", 
                                         #labels = sniptest[which(sniptest$Genus == "Dcitri" & sniptest$Strain_ID == "metasample"),]$Product, 
                                         minradius = 0.45, 
                                         maxRadius = 0.99, 
                                         shape = 'circle',
                                         size = 2.6)

tracklist <- tracklist + snplist
tracklist <- tracklist + BioCircosTextTrack('Genus', '\'Diaphorina citri\'',
                                            x = -0.45, y = -1.32, color = 'black', weight = 'bold')

tracklist <- tracklist + CDStrack + tRNAtrack + rRNAtrack + Intergenic_track
#tracklist <- tracklist + BioCircosLineTrack('gcskew', gcskew_cov$Contig, gcskew_cov$Position,
#             gcskew_cov$GC_skew, color = '#40B9D4', width = 1, maxRadius = 1, minRadius = 0.9)

BioCircos(genome = Dcit, genomeFillColor = c("grey95", "grey95"), genomeTicksScale = 1e3, displayGenomeBorder = T,
          genomeBorderColor = 'black', zoom = T,  genomeBorderSize = 0.2,
          genomeLabelDx = 3.1, genomeLabelDy = 25, genomeLabelOrientation = -1,
          genomeTicksLen = 3.5, genomeTicksTextSize = "1.5em", tracklist)

#################################################################################################
## Start Thrips

setwd("C://Bioinformatics/Heck_USDA_Pathology/Dcitri_mtDNA_haplotyping/Examining_the_effect_of_kmer_size/ska2/Rplots")

snippy <- read.table("Thrips-k13-m50-mapped-OP546494.1-snps-tab-for-R-input-snp-proportions.txt", sep = "\t", stringsAsFactors = F, na.strings = "NA")

colnames(snippy) <- c("Strain_ID", "Genus", "Chrom", "Pos", "Type", "Ref", "Alt", "Evidence", "FType", "Strand", "NT_pos", "AA_pos", "Effect", "Locus_tag", "Gene", "Product")

snippy[is.na(snippy)] <- ''

sniptest <- snippy
newprotsamps <- c("KY426014","KY426015","MF181946","MF181947","metasample","MF614804","MF614805","MF614806","MF614807","MF614808","MF614809","MF614810","MF614811","MF614812","MF614813","MF614814","MF614815","MF614816","MF614817","MF614818","MF614819","MF614820","MF614821","MF614822","MF614823","MF614824","MF614825","MF614826","MF614827","MF614828", "OM181946", "MF181947")
sniptest$protocol <- ifelse(sniptest$Strain_ID %in% newprotsamps, 'new_protocol', 'existing_protocol')

# use a high coverage CLas strain KY426014 as the example to plot

sniplen <- length(sniptest[which(sniptest$Genus == "Fintonsa" & sniptest$Strain_ID == "metasample"),]$Pos)

# create a new column in sniptest to hold the read depth reported in the "Evidence"
# column by snippy
Evidence_depth <- c()
for(i in sniptest$Evidence){
  Evidence_depth <- c(Evidence_depth, as.numeric(strsplit(strsplit(strsplit(strsplit(i, " ")[[1]], " ")[[1]],":")[[1]], " ")[[2]]))
  print(as.numeric(strsplit(strsplit(strsplit(strsplit(i, " ")[[1]], " ")[[1]],":")[[1]], " ")[[2]]))
}

sniptest$SNP_depth <- Evidence_depth
# create a variable to change shape of SNP value if SNP inside a coding sequence
# or intergenic region
#sniptest$shape <- ifelse(sniptest$FType == "CDS", 'circle', 'rect')

# BioCircos plots

# genome track
Fint <- list('OP546494.1' = 15220)

Fint_mtDNA_CDS_arc_labels_weitz = c('nad2','nad1','atp8','atp6','nad5','nad4','nad4l','nad6','cox1','nad3','cox2','cox3','cob')
Fint_mtDNA_CDS_arc_labels = c('nad2','nad1','atp8','atp6','nad5','nad4','nad4l','nad6','cox1','nad3','cox2','cox3','cob')
Fint_mtDNA_CDS_arc_chroms = c('OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1')
Fint_mtDNA_CDS_arc_start <- c(159,1197,2933,3147,4349,6086,7385,7761,9471,11027,11445,12385,13792)
Fint_mtDNA_CDS_arc_stop <- c(1131,2121,3130,3842,6022,7433,7663,8237,11024,11377,12104,13173,14904)

Fint_mtDNA_tRNA_arc_labels_weitz = c('tRNA-Pro','tRNA-Tyr','tRNA-Trp','tRNA-Met','tRNA-Ala','tRNA-Phe','tRNA-Asn','tRNA-Glu','tRNA-Ser','tRNA-Leu','tRNA-His','tRNA-Cys','tRNA-Val','tRNA-Ser','tRNA-Leu','tRNA-Asp','tRNA-Arg','tRNA-Gly','tRNA-Lys','tRNA-Ile','tRNA-Thr','tRNA-Gln')
Fint_mtDNA_tRNA_arc_labels = c('tRNA-Pro','tRNA-Tyr','tRNA-Trp','tRNA-Met','tRNA-Ala','tRNA-Phe','tRNA-Asn','tRNA-Glu','tRNA-Ser','tRNA-Leu','tRNA-His','tRNA-Cys','tRNA-Val','tRNA-Ser','tRNA-Leu','tRNA-Asp','tRNA-Arg','tRNA-Gly','tRNA-Lys','tRNA-Ile','tRNA-Thr','tRNA-Gln')
Fint_mtDNA_tRNA_arc_chroms = c('OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1')
Fint_mtDNA_tRNA_arc_start <- c(1,64,1132,2122,2185,2247,3851,3912,3974,4029,6023,7683,8270,9403,11379,12108,12175,12244,12306,13191,13275,14916)
Fint_mtDNA_tRNA_arc_stop <- c(64,125,1196,2183,2246,2310,3913,3974,4029,4092,6085,7743,8324,9469,11444,12174,12240,12306,12367,13256,13336,14983)

Fint_mtDNA_rRNA_arc_labels_weitz = c('12S_rRNA', '16S_rRNA')
Fint_mtDNA_rRNA_arc_labels = c('12S_rRNA', '16S_rRNA')
Fint_mtDNA_rRNA_arc_chroms = c('OP546494.1', 'OP546494.1')
Fint_mtDNA_rRNA_arc_start <- c(2311, 8325)
Fint_mtDNA_rRNA_arc_stop <- c(2932, 9402)

Fint_mtDNA_intergenic_arc_labels_weitz = c('int1', 'int2', 'int3', 'int4', 'int5', 'int6', 'int7', 'int8', 'int9', 'int10', 'int11', 'int12', 'int13', 'int14', 'int15', 'int16', 'int17', 'int18', 'int19')
Fint_mtDNA_intergenic_arc_labels = c('int1', 'int2', 'int3', 'int4', 'int5', 'int6', 'int7', 'int8', 'int9', 'int10', 'int11', 'int12', 'int13', 'int14', 'int15', 'int16', 'int17', 'int18', 'int19')
Fint_mtDNA_intergenic_arc_chroms = c('OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1', 'OP546494.1')
Fint_mtDNA_intergenic_arc_start <- c(126, 2184, 3131, 3843, 4093, 7664, 7744, 8238, 9470, 11025, 11378, 12105, 12241, 12368, 13174, 13257, 13337, 14905, 14984)
Fint_mtDNA_intergenic_arc_stop <- c(158, 2184, 3146, 3850, 4348, 7682, 7760, 8269, 9470, 11026, 11378, 12107, 12243, 12384, 13190, 13274, 13791, 14915, 15220)


CDStrack = BioCircosArcTrack('Fint_CDS', Fint_mtDNA_CDS_arc_chroms, Fint_mtDNA_CDS_arc_start, Fint_mtDNA_CDS_arc_stop,
                             minRadius = 0.45, maxRadius = 1, opacities = 0.25, labels = Fint_mtDNA_CDS_arc_labels, colors = c("#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#35B5A5","#35B5A5"))
tRNAtrack = BioCircosArcTrack('Fint_tRNA', Fint_mtDNA_tRNA_arc_chroms, Fint_mtDNA_tRNA_arc_start, Fint_mtDNA_tRNA_arc_stop,
                              minRadius = 0.45, maxRadius = 1, opacities = 0.8, labels = Fint_mtDNA_tRNA_arc_labels, colors = "#FFFF92")
rRNAtrack = BioCircosArcTrack('Fint_rRNA', Fint_mtDNA_rRNA_arc_chroms, Fint_mtDNA_rRNA_arc_start, Fint_mtDNA_rRNA_arc_stop,
                              minRadius = 0.45, maxRadius = 1, opacities = 0.8, labels = Fint_mtDNA_rRNA_arc_labels, colors = "#FD6C5C")
Intergenic_track = BioCircosArcTrack('Fint_int', Fint_mtDNA_intergenic_arc_chroms, Fint_mtDNA_intergenic_arc_start, Fint_mtDNA_intergenic_arc_stop,
                                     minRadius = 0.45, maxRadius = 1, opacities = 0.8, labels = Fint_mtDNA_intergenic_arc_labels, colors = "#B3ACDC")
#tracklist = CDStrack

# create a link track to indicate gene names of SNPs
link1 <- rep('OP546494.1', sniplen)
link2 <- rep('OP546494.1', sniplen)
linkpos1 <- sniptest[which(sniptest$Genus == "Fintonsa" & sniptest$Strain_ID == "metasample"),]$Pos
linkpos2 <- sniptest[which(sniptest$Genus == "Fintonsa" & sniptest$Strain_ID == "metasample"),]$Pos + 1
linklabs <- sniptest[which(sniptest$Genus == "Fintonsa" & sniptest$Strain_ID == "metasample"),]$Gene

trackback = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 0.46,
                                     borderSize = 0, fillColors = "#EEFFEE") 


trackback2 = BioCircosBackgroundTrack("mySNPBackgroundTrack", minRadius = 0.45, maxRadius = 1,
                                      borderSize = 0.1, borderColors = "black", fillColors = "snow") 

snpColors <- brewer.pal(length(levels(as.factor(sniptest[which(sniptest$Genus == "Fintonsa" & sniptest$Strain_ID == "metasample"),]$Type))),"Dark2")
names(snpColors) <- levels(as.factor(sniptest[which(sniptest$Genus == "Fintonsa" & sniptest$Strain_ID == "metasample"),]$Type))
snpColors2 <- snpColors[sniptest[which(sniptest$Genus == "Fintonsa" & sniptest$Strain_ID == "metasample"),]$Type]

tracklist = BioCircosLinkTrack('snps', link1, link2,
                               linkpos1, link2, linkpos1, linkpos2,
                               maxRadius = 0.3, #labels = linklabs,
                               axisColor = 'white', labelColor = as.vector(snpColors2))

snplist = trackback2 + BioCircosSNPTrack('snptrack', link1, linkpos1, 
                                         sniptest[which(sniptest$Genus == "Fintonsa" & sniptest$Strain_ID == "metasample"),]$SNP_depth, 
                                         colors = "gray20", 
                                         labels = sniptest[which(sniptest$Genus == "Fintonsa" & sniptest$Strain_ID == "metasample"),]$Product, 
                                         minradius = 0.46, 
                                         maxRadius = 0.99, 
                                         shape = 'circle',
                                         size = 2)

tracklist <- tracklist + snplist
tracklist <- tracklist + BioCircosTextTrack('Genus', '\'Frankliniella intonsa\'',
                                            x = -0.45, y = -1.32, color = 'black', weight = 'bold')

tracklist <- tracklist + CDStrack + tRNAtrack + rRNAtrack + Intergenic_track
#tracklist <- tracklist + BioCircosLineTrack('gcskew', gcskew_cov$Contig, gcskew_cov$Position,
#             gcskew_cov$GC_skew, color = '#40B9D4', width = 1, maxRadius = 1, minRadius = 0.9)

BioCircos(genome = Fint, genomeFillColor = c("grey95", "grey95"), genomeTicksScale = 1e3, displayGenomeBorder = T,
          genomeBorderColor = 'black', zoom = T,  genomeBorderSize = 0.2,
          genomeLabelDx = 3.1, genomeLabelDy = 25, genomeLabelOrientation = -1,
          genomeTicksLen = 3.5, genomeTicksTextSize = "1.5em", tracklist)

#################################################################################################
## Start Ovis

snippy <- read.table("Ovis-k15-m90-mapped-JN181255.1-snps-tab-for-R-input-snp-proportions.txt", sep = "\t", stringsAsFactors = F, na.strings = "NA")

colnames(snippy) <- c("Strain_ID", "Genus", "Chrom", "Pos", "Type", "Ref", "Alt", "Evidence", "FType", "Strand", "NT_pos", "AA_pos", "Effect", "Locus_tag", "Gene", "Product")

snippy[is.na(snippy)] <- ''

sniptest <- snippy
newprotsamps <- c("KY426014","KY426015","MF181946","MF181947","metasample","MF614804","MF614805","MF614806","MF614807","MF614808","MF614809","MF614810","MF614811","MF614812","MF614813","MF614814","MF614815","MF614816","MF614817","MF614818","MF614819","MF614820","MF614821","MF614822","MF614823","MF614824","MF614825","MF614826","MF614827","MF614828", "OM181946", "MF181947")
sniptest$protocol <- ifelse(sniptest$Strain_ID %in% newprotsamps, 'new_protocol', 'existing_protocol')

# use a high coverage CLas strain KY426014 as the example to plot

sniplen <- length(sniptest[which(sniptest$Genus == "Ovis" & sniptest$Strain_ID == "metasample"),]$Pos)

# create a new column in sniptest to hold the read depth reported in the "Evidence"
# column by snippy
Evidence_depth <- c()
for(i in sniptest$Evidence){
  Evidence_depth <- c(Evidence_depth, as.numeric(strsplit(strsplit(strsplit(strsplit(i, " ")[[1]], " ")[[1]],":")[[1]], " ")[[2]]))
  print(as.numeric(strsplit(strsplit(strsplit(strsplit(i, " ")[[1]], " ")[[1]],":")[[1]], " ")[[2]]))
}

sniptest$SNP_depth <- Evidence_depth
# create a variable to change shape of SNP value if SNP inside a coding sequence
# or intergenic region
#sniptest$shape <- ifelse(sniptest$FType == "CDS", 'circle', 'rect')

# BioCircos plots

# genome track
Ovis <- list('JN181255.1' = 16463)

Ovis_mtDNA_CDS_arc_labels_weitz = c('ND1','ND2','COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CYTB')
Ovis_mtDNA_CDS_arc_labels = c('ND1','ND2','COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CYTB')
Ovis_mtDNA_CDS_arc_chroms = c('JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1')
Ovis_mtDNA_CDS_arc_start <- c(2743,3908,5330,7017,7773,7934,8614,9467,9883,10173,11751,13555,14156)
Ovis_mtDNA_CDS_arc_stop <- c(3698,4949,6874,7700,7973,8614,9397,9812,10179,11550,13571,14082,15295)

Ovis_mtDNA_tRNA_arc_labels_weitz = c('tRNA-Phe','tRNA-Val','tRNA-Leu','tRNA-Ile','tRNA-Gln','tRNA-Met','tRNA-Trp','tRNA-Ala','tRNA-Asn','tRNA-Cys','tRNA-Tyr','tRNA-Ser','tRNA-Asp','tRNA-Lys','tRNA-Gly','tRNA-Arg','tRNA-His','tRNA-Ser','tRNA-Leu','tRNA-Glu','tRNA-Thr','tRNA-Pro')
Ovis_mtDNA_tRNA_arc_labels = c('tRNA-Phe','tRNA-Val','tRNA-Leu','tRNA-Ile','tRNA-Gln','tRNA-Met','tRNA-Trp','tRNA-Ala','tRNA-Asn','tRNA-Cys','tRNA-Tyr','tRNA-Ser','tRNA-Asp','tRNA-Lys','tRNA-Gly','tRNA-Arg','tRNA-His','tRNA-Ser','tRNA-Leu','tRNA-Glu','tRNA-Thr','tRNA-Pro')
Ovis_mtDNA_tRNA_arc_chroms = c('JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1')
Ovis_mtDNA_tRNA_arc_start <- c(1,1026,2666,3699,3765,3839,4950,5018,5088,5193,5261,6872,6948,7704,9398,9814,11551,11620,11681,14083,15299,15368)
Ovis_mtDNA_tRNA_arc_stop <- c(68,1092,2740,3767,3836,3907,5016,5086,5160,5260,5328,6942,7015,7771,9466,9882,11619,11679,11750,14151,15368,15433)

Ovis_mtDNA_rRNA_arc_labels_weitz = c('12S_rRNA', '16S_rRNA')
Ovis_mtDNA_rRNA_arc_labels = c('12S_rRNA', '16S_rRNA')
Ovis_mtDNA_rRNA_arc_chroms = c('JN181255.1', 'JN181255.1')
Ovis_mtDNA_rRNA_arc_start <- c(69, 1093)
Ovis_mtDNA_rRNA_arc_stop <- c(1025, 2265)

Ovis_mtDNA_intergenic_arc_labels_weitz = c('int1', 'int2', 'int3', 'int4', 'int5', 'int6', 'int7', 'int8', 'int9', 'int10', 'int11', 'int12', 'int13', 'int14', 'int15', 'int16')
Ovis_mtDNA_intergenic_arc_labels = c('int1', 'int2', 'int3', 'int4', 'int5', 'int6', 'int7', 'int8', 'int9', 'int10', 'int11', 'int12', 'int13', 'int14', 'int15', 'int16')
Ovis_mtDNA_intergenic_arc_chroms = c('JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1', 'JN181255.1')
Ovis_mtDNA_intergenic_arc_start <- c(2266, 2741, 3837, 5017, 5087, 5161, 5329, 6943, 7016, 7701, 7772, 9813, 11680, 14152, 15296, 15434)
Ovis_mtDNA_intergenic_arc_stop <- c(2665, 2742, 3838, 5017, 5087, 5192, 5329, 6947, 7016, 7703, 7772, 9813, 11680, 14155, 15298, 16463)


CDStrack = BioCircosArcTrack('Ovis_CDS', Ovis_mtDNA_CDS_arc_chroms, Ovis_mtDNA_CDS_arc_start, Ovis_mtDNA_CDS_arc_stop,
                             minRadius = 0.45, maxRadius = 1, opacities = 0.25, labels = Ovis_mtDNA_CDS_arc_labels, colors = c("#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#35B5A5","#35B5A5"))
tRNAtrack = BioCircosArcTrack('Ovis_tRNA', Ovis_mtDNA_tRNA_arc_chroms, Ovis_mtDNA_tRNA_arc_start, Ovis_mtDNA_tRNA_arc_stop,
                              minRadius = 0.45, maxRadius = 1, opacities = 0.8, labels = Ovis_mtDNA_tRNA_arc_labels, colors = "#FFFF92")
rRNAtrack = BioCircosArcTrack('Ovis_rRNA', Ovis_mtDNA_rRNA_arc_chroms, Ovis_mtDNA_rRNA_arc_start, Ovis_mtDNA_rRNA_arc_stop,
                              minRadius = 0.45, maxRadius = 1, opacities = 0.8, labels = Ovis_mtDNA_rRNA_arc_labels, colors = "#FD6C5C")
Intergenic_track = BioCircosArcTrack('Ovis_int', Ovis_mtDNA_intergenic_arc_chroms, Ovis_mtDNA_intergenic_arc_start, Ovis_mtDNA_intergenic_arc_stop,
                                     minRadius = 0.45, maxRadius = 1, opacities = 0.8, labels = Ovis_mtDNA_intergenic_arc_labels, colors = "#B3ACDC")

#tracklist = CDStrack

# create a link track to indicate gene names of SNPs
link1 <- rep('JN181255.1', sniplen)
link2 <- rep('JN181255.1', sniplen)
linkpos1 <- sniptest[which(sniptest$Genus == "Ovis" & sniptest$Strain_ID == "metasample"),]$Pos
linkpos2 <- sniptest[which(sniptest$Genus == "Ovis" & sniptest$Strain_ID == "metasample"),]$Pos + 1
linklabs <- sniptest[which(sniptest$Genus == "Ovis" & sniptest$Strain_ID == "metasample"),]$Gene

trackback = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 0.46,
                                     borderSize = 0, fillColors = "#EEFFEE") 


trackback2 = BioCircosBackgroundTrack("mySNPBackgroundTrack", minRadius = 0.45, maxRadius = 1,
                                      borderSize = 0.1, borderColors = "black", fillColors = "snow") 

snpColors <- brewer.pal(length(levels(as.factor(sniptest[which(sniptest$Genus == "Ovis" & sniptest$Strain_ID == "metasample"),]$Type))),"Dark2")
names(snpColors) <- levels(as.factor(sniptest[which(sniptest$Genus == "Ovis" & sniptest$Strain_ID == "metasample"),]$Type))
snpColors2 <- snpColors[sniptest[which(sniptest$Genus == "Ovis" & sniptest$Strain_ID == "metasample"),]$Type]

tracklist = BioCircosLinkTrack('snps', link1, link2,
                               linkpos1, link2, linkpos1, linkpos2,
                               maxRadius = 0.3, #labels = linklabs,
                               axisColor = 'white', labelColor = as.vector(snpColors2))

snplist = trackback2 + BioCircosSNPTrack('snptrack', link1, linkpos1, 
                                         sniptest[which(sniptest$Genus == "Ovis" & sniptest$Strain_ID == "metasample"),]$SNP_depth, 
                                         colors = "gray20", 
                                         labels = sniptest[which(sniptest$Genus == "Ovis" & sniptest$Strain_ID == "metasample"),]$Product, 
                                         minradius = 0.46, 
                                         maxRadius = 0.99, 
                                         shape = 'circle',
                                         size = 2)

tracklist <- tracklist + snplist
tracklist <- tracklist + BioCircosTextTrack('Genus', '\'Ovis spp.\'',
                                            x = -0.45, y = -1.32, color = 'black', weight = 'bold')

tracklist <- tracklist + CDStrack + tRNAtrack + rRNAtrack + Intergenic_track
#tracklist <- tracklist + BioCircosLineTrack('gcskew', gcskew_cov$Contig, gcskew_cov$Position,
#             gcskew_cov$GC_skew, color = '#40B9D4', width = 1, maxRadius = 1, minRadius = 0.9)

BioCircos(genome = Ovis, genomeFillColor = c("grey95", "grey95"), genomeTicksScale = 1e3, displayGenomeBorder = T,
          genomeBorderColor = 'black', zoom = T,  genomeBorderSize = 0.2,
          genomeLabelDx = 3.1, genomeLabelDy = 25, genomeLabelOrientation = -1,
          genomeTicksLen = 3.5, genomeTicksTextSize = "1.5em", tracklist)

#################################################################################################
## Start Combined

snippy <- read.table("Combined-k15-m90-mapped-MF614822.1-snps-tab-for-R-input-snp-proportions.txt", sep = "\t", stringsAsFactors = F, na.strings = "NA")

colnames(snippy) <- c("Strain_ID", "Genus", "Chrom", "Pos", "Type", "Ref", "Alt", "Evidence", "FType", "Strand", "NT_pos", "AA_pos", "Effect", "Locus_tag", "Gene", "Product")

snippy[is.na(snippy)] <- ''

sniptest <- snippy
newprotsamps <- c("KY426014","KY426015","MF181946","MF181947","metasample","MF614804","MF614805","MF614806","MF614807","MF614808","MF614809","MF614810","MF614811","MF614812","MF614813","MF614814","MF614815","MF614816","MF614817","MF614818","MF614819","MF614820","MF614821","MF614822","MF614823","MF614824","MF614825","MF614826","MF614827","MF614828", "OM181946", "MF181947")
sniptest$protocol <- ifelse(sniptest$Strain_ID %in% newprotsamps, 'new_protocol', 'existing_protocol')

# use a high coverage CLas strain KY426014 as the example to plot

sniplen <- length(sniptest[which(sniptest$Genus == "Combined" & sniptest$Strain_ID == "metasample"),]$Pos)

# create a new column in sniptest to hold the read depth reported in the "Evidence"
# column by snippy
Evidence_depth <- c()
for(i in sniptest$Evidence){
  Evidence_depth <- c(Evidence_depth, as.numeric(strsplit(strsplit(strsplit(strsplit(i, " ")[[1]], " ")[[1]],":")[[1]], " ")[[2]]))
  print(as.numeric(strsplit(strsplit(strsplit(strsplit(i, " ")[[1]], " ")[[1]],":")[[1]], " ")[[2]]))
}

sniptest$SNP_depth <- Evidence_depth
# create a variable to change shape of SNP value if SNP inside a coding sequence
# or intergenic region
#sniptest$shape <- ifelse(sniptest$FType == "CDS", 'circle', 'rect')

# BioCircos plots

# genome track
Combined <- list('MF614822.1' = 15020)

Combined_mtDNA_CDS_arc_labels_weitz = c('nad2','cox1','cox2','atp8','atp6','cox3','nad3','nad5','nad4','nad4L','nad6','cob','nad1')
Combined_mtDNA_CDS_arc_labels = c('nad2','cox1','cox2','atp8','atp6','cox3','nad3','nad5','nad4','nad4L','nad6','cob','nad1')
Combined_mtDNA_CDS_arc_chroms = c('MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1')
Combined_mtDNA_CDS_arc_start <- c(203,1360,2954,3747,3893,4567,5423,6127,7804,9042,9461,9937,11162)
Combined_mtDNA_CDS_arc_stop <- c(1174,2889,3617,3899,4567,5349,5773,7747,9048,9317,9943,11079,12073)

Combined_mtDNA_tRNA_arc_labels_weitz = c('tRNA-Ile','tRNA-Gln','tRNA-Met','tRNA-Trp','tRNA-Cys','tRNA-Tyr','tRNA-Leu','tRNA-Lys','tRNA-Asp','tRNA-Gly','tRNA-Ala','tRNA-Arg','tRNA-Asn','tRNA-Ser','tRNA-Glu','tRNA-Phe','tRNA-His','tRNA-Thr','tRNA-Pro','tRNA-Ser','tRNA-Leu','tRNA-Val')
Combined_mtDNA_tRNA_arc_labels = c('tRNA-Ile','tRNA-Gln','tRNA-Met','tRNA-Trp','tRNA-Cys','tRNA-Tyr','tRNA-Leu','tRNA-Lys','tRNA-Asp','tRNA-Gly','tRNA-Ala','tRNA-Arg','tRNA-Asn','tRNA-Ser','tRNA-Glu','tRNA-Phe','tRNA-His','tRNA-Thr','tRNA-Pro','tRNA-Ser','tRNA-Leu','tRNA-Val')
Combined_mtDNA_tRNA_arc_chroms = c('MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1')
Combined_mtDNA_tRNA_arc_start <- c(1,71,138,1173,1234,1295,2893,3618,3686,5362,5772,5832,5893,5959,6015,6065,7745,9331,9393,11078,12074,13272)
Combined_mtDNA_tRNA_arc_stop <- c(65,137,202,1233,1295,1356,2953,3687,3746,5424,5831,5892,5959,6014,6075,6126,7804,9394,9458,11140,12137,13335)

Combined_mtDNA_rRNA_arc_labels_weitz = c('12S_rRNA', '16S_rRNA')
Combined_mtDNA_rRNA_arc_labels = c('12S_rRNA', '16S_rRNA')
Combined_mtDNA_rRNA_arc_chroms = c('MF614822.1', 'MF614822.1')
Combined_mtDNA_rRNA_arc_start <- c(13336, 12138)
Combined_mtDNA_rRNA_arc_stop <- c(14093, 13271)

Combined_mtDNA_intergenic_arc_labels_weitz = c('int1', 'int2', 'int3', 'int4', 'int5', 'int6', 'int7', 'int8')
Combined_mtDNA_intergenic_arc_labels = c('int1', 'int2', 'int3', 'int4', 'int5', 'int6', 'int7', 'int8')
Combined_mtDNA_intergenic_arc_chroms = c('MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1', 'MF614822.1')
Combined_mtDNA_intergenic_arc_start <- c(66, 1357, 2890, 5350, 9318, 9459, 11141, 14094)
Combined_mtDNA_intergenic_arc_stop <- c(70, 1359, 2892, 5361, 9330, 9460, 11161, 15020)


CDStrack = BioCircosArcTrack('Combined_CDS', Combined_mtDNA_CDS_arc_chroms, Combined_mtDNA_CDS_arc_start, Combined_mtDNA_CDS_arc_stop,
                             minRadius = 0.45, maxRadius = 1, opacities = 0.25, labels = Combined_mtDNA_CDS_arc_labels, colors = c("#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#6DCBBA","#35B5A5","#35B5A5"))
tRNAtrack = BioCircosArcTrack('Combined_tRNA', Combined_mtDNA_tRNA_arc_chroms, Combined_mtDNA_tRNA_arc_start, Combined_mtDNA_tRNA_arc_stop,
                              minRadius = 0.45, maxRadius = 1, opacities = 0.8, labels = Combined_mtDNA_tRNA_arc_labels, colors = "#FFFF92")
rRNAtrack = BioCircosArcTrack('Combined_rRNA', Combined_mtDNA_rRNA_arc_chroms, Combined_mtDNA_rRNA_arc_start, Combined_mtDNA_rRNA_arc_stop,
                              minRadius = 0.45, maxRadius = 1, opacities = 0.8, labels = Combined_mtDNA_rRNA_arc_labels, colors = "#FD6C5C")
Intergenic_track = BioCircosArcTrack('Combined_int', Combined_mtDNA_intergenic_arc_chroms, Combined_mtDNA_intergenic_arc_start, Combined_mtDNA_intergenic_arc_stop,
                                     minRadius = 0.45, maxRadius = 1, opacities = 0.8, labels = Combined_mtDNA_intergenic_arc_labels, colors = "#B3ACDC")

#tracklist = CDStrack

# create a link track to indicate gene names of SNPs
link1 <- rep('MF614822.1', sniplen)
link2 <- rep('MF614822.1', sniplen)
linkpos1 <- sniptest[which(sniptest$Genus == "Combined" & sniptest$Strain_ID == "metasample"),]$Pos
linkpos2 <- sniptest[which(sniptest$Genus == "Combined" & sniptest$Strain_ID == "metasample"),]$Pos + 1
linklabs <- sniptest[which(sniptest$Genus == "Combined" & sniptest$Strain_ID == "metasample"),]$Gene

trackback = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 0.46,
                                     borderSize = 0, fillColors = "#EEFFEE") 


trackback2 = BioCircosBackgroundTrack("mySNPBackgroundTrack", minRadius = 0.45, maxRadius = 1,
                                      borderSize = 0.1, borderColors = "black", fillColors = "snow") 

snpColors <- brewer.pal(length(levels(as.factor(sniptest[which(sniptest$Genus == "Combined" & sniptest$Strain_ID == "metasample"),]$Type))),"Dark2")
names(snpColors) <- levels(as.factor(sniptest[which(sniptest$Genus == "Combined" & sniptest$Strain_ID == "metasample"),]$Type))
snpColors2 <- snpColors[sniptest[which(sniptest$Genus == "Combined" & sniptest$Strain_ID == "metasample"),]$Type]

tracklist = BioCircosLinkTrack('snps', link1, link2,
                               linkpos1, link2, linkpos1, linkpos2,
                               maxRadius = 0.3, #labels = linklabs,
                               axisColor = 'white', labelColor = as.vector(snpColors2))

snplist = trackback2 + BioCircosSNPTrack('snptrack', link1, linkpos1, 
                                         sniptest[which(sniptest$Genus == "Combined" & sniptest$Strain_ID == "metasample"),]$SNP_depth, 
                                         colors = "gray20", 
                                         #labels = sniptest[which(sniptest$Genus == "Combined" & sniptest$Strain_ID == "metasample"),]$Product, 
                                         minradius = 0.46, 
                                         maxRadius = 0.99, 
                                         shape = 'circle',
                                         size = 2)

tracklist <- tracklist + snplist
tracklist <- tracklist + BioCircosTextTrack('Genus', '\'Diaphorina citri\'',
                                            x = -0.45, y = -1.32, color = 'black', weight = 'bold')

tracklist <- tracklist + CDStrack + tRNAtrack + rRNAtrack + Intergenic_track
#tracklist <- tracklist + BioCircosLineTrack('gcskew', gcskew_cov$Contig, gcskew_cov$Position,
#             gcskew_cov$GC_skew, color = '#40B9D4', width = 1, maxRadius = 1, minRadius = 0.9)

BioCircos(genome = Combined, genomeFillColor = c("grey95", "grey95"), genomeTicksScale = 1e3, displayGenomeBorder = T,
          genomeBorderColor = 'black', zoom = T,  genomeBorderSize = 0.2,
          genomeLabelDx = 3.1, genomeLabelDy = 25, genomeLabelOrientation = -1,
          genomeTicksLen = 3.5, genomeTicksTextSize = "1.4em", tracklist)

