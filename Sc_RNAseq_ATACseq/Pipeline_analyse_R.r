##############################
########### R ################
##############################

#En premier lieu nous allons effectuer le controle qualité.
workdir = "/work/user/ocuvier/scMultiome"
mesdatas = "/work/user/ocuvier/scMultiome/Beaf32/outs/filtered_feature_bc_matrix.h5"

library(Seurat)
library(Signac)
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(ggplot2)
library(cowplot)
library(magrittr)
library(tidyverse)
library(rtracklayer)


counts = Read10X_h5("/work/user/ocuvier/scMultiome/Beaf32/outs/filtered_feature_bc_matrix.h5")

#On obtient une liste 
head(names(counts))
#on obtient Gene Expression" "Peaks"

all(colnames(counts[[1]])==colnames(counts[[2]]))
# True ça nous montre qu'on a bien les mêmes cellules dans le RNAseq et dans l'ATACseq

#On crée l'objet Seurat pour faire plein de choses dessus

fragpath <- "/work/user/ocuvier/scMultiome/Beaf32/outs/atac_fragments.tsv.gz"
annotation <- import(paste0("/work/user/ocuvier/scMultiome/genes/genes.gtf.gz"))

Singlecell = CreateSeuratObject(
    counts = counts$'Gene Expression',
    assay = "RNAseq"
)

Singlecell[["ATAC"]] = CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":","-"),
    fragments = fragpath,
    annotation = annotation
)
rm(counts) #Remove counts
str(Singlecell)

#Singlecell@assays$RNAseq@counts  contient les counts

colnames(Singlecell@meta.data)

#Mesure du pourcentage du RNA mitochondrial.

DefaultAssay(Singlecell) = "RNAseq"
mito.genes = grep(pattern = "^mt:",
    x = rownames(x = Singlecell@assays$RNAseq@counts), value = TRUE)
percent.mito = Matrix::colSums(Singlecell@assays$RNAseq@counts[mito.genes, ])/
    Matrix::colSums(Singlecell@assays$RNAseq@counts)
Singlecell = AddMetaData(object = Singlecell, metadata = percent.mito, col.name = "percent.mt")


#On regarde les sites TSS = des emplacements spécifiques sur un brin d'ADN où la transcription 
#(le processus de copie de l'ADN en ARN) commence. Les sites de début de transcription sont cruciaux 
#dans l'expression génique, car ils marquent le début des instructions pour créer une protéine spécifique 
#ou une molécule d'ARN.

DefaultAssay(Singlecell) = "ATAC"
Singlecell = NucleosomeSignal(object = Singlecell)
Singlecell = TSSEnrichment(object=Singlecell, fast = TRUE)

#saveRDS(Singlecell,"/home/ocuvier/work/conda/MData/save_data/Singlecell_preprocessed.RDATA")
Singlecell = readRDS("/home/ocuvier/work/conda/MData/save_data/Singlecell_preprocessed.RDATA")

pdf("/home/ocuvier/work/conda/MData/figures/figure1.pdf")
DensityScatter(Singlecell, x = 'nCount_ATAC', y ='TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

Singlecell$high.tss = ifelse(Singlecell$TSS.enrichment > 1 , 'High', 'Low')
Singlecell$high.tss


### A VOIR ###ss
# summary(Singlecell$nucleosome_signal)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.1111  0.8957  1.0211  0.9983  1.1248  2.7143 
 
Singlecell$nucleosome_group = ifelse(Singlecell$nucleosome_signal > 1.3 , 'NS > 1.3', 'NS < 1.3')
pdf("/home/ocuvier/work/conda/MData/figures/figure2bis.pdf")
FragmentHistogram(object = Singlecell, group.by = 'nucleosome_group',region="2L-1-20000")
dev.off()

#saveRDS(Singlecell,"/home/ocuvier/work/conda/MData/save_data/Singlecell_qualitymeasurements.RDATA")
Singlecell = readRDS("/home/ocuvier/work/conda/MData/save_data/Singlecell_qualitymeasurements.RDATA")


# Comme la moyenne et mediane est à 1 il y a moins de nucléosomes que ce qu'on attend. Le génome droso est plus petit et plus dense en gène 
# Pour le TSS => ce qu'on compte comme le background est peut etre le signal voisin. 

# Normalement il y 14000 pic atac. et peut etre qu'on a moins de pic mais que ce sont de bons pics. et qu'il y a eu peu de transposase ajoutée. Une concentration plus faible que ce qui se met d'habitude = pas de saturation des cellules. 
# Les métriques de qualités peuvent conduire à d'autres supposition. 

# je peux exporter le fichier dans python pour langer les commandes import et tout... et la sortie de ces commandes c'est un vecteurs ou un tableau à importer dans R.


#SUR PYTHON 
import scrublet as scr
import scipy.io
import numpy as np
import os
import pandas as pd

input_dir = "/work/user/ocuvier/scMultiome/Beaf32/outs/"
counts_matrix = scipy.io.mmread(input_dir + "filtered_feature_bc_matrix/matrix.mtx.gz").T.tocsc()

scrub = scr.Scrublet(counts_matrix,
expected_doublet_rate=0.1,
sim_doublet_ratio=2,
n_neighbors = 8)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=1, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=25)

# Preprocessing...
# Simulating doublets...
# Embedding transcriptomes using PCA...
# Calculating doublet scores...
# Automatically set threshold at doublet score = 0.35 ## Sur la ligne 2 du goem vline
# Detected doublet rate = 4.7%
# Estimated detectable doublet fraction = 13.9%
# Overall doublet rate:
# 	Expected   = 10.0%
# 	Estimated  = 34.2%
# Elapsed time: 45.5 seconds

#Il faut exporter sur csv 

# Créer un DataFrame à partir de doublet_scores et predicted_doublets
df = pd.DataFrame({'DoubletScore': doublet_scores, 'PredictedDoublet': predicted_doublets})

# Définir le chemin d'accès et le nom du fichier CSV
output_file = "/work/user/ocuvier/scMultiome/Beaf32/Doublet_Scores.csv"

# Enregistrer le DataFrame dans un fichier CSV
df.to_csv(output_file, index=False)

print("Les résultats du Scrublet ont été enregistrés dans le fichier CSV:", output_file)

#Les résultats du Scrublet ont été enregistrés dans le fichier CSV: /work/user/ocuvier/scMultiome/Beaf32/Doublet_Scores.csv



### SUR R ###

py = read.csv("/work/user/ocuvier/scMultiome/Beaf32/Doublet_Scores.csv")


Singlecell@meta.data$DoubletScore = py$DoubletScore

library(ggplot2)

pdf("/home/ocuvier/work/conda/MData/figures/figure3.pdf")
ggplot(Singlecell@meta.data, aes(x = DoubletScore, after_stat(ndensity))) +
  geom_histogram(bins = 200, colour ="lightgrey")+
  geom_vline(xintercept = 0.35, colour = "red", linetype = 2)+
  geom_vline(xintercept = 0.3, colour = "green", linetype = 2) 
dev.off()


#Pour détecter le nombre de doublet, je laisse le seuil à 0.35 (j'ai également essayé 0.3 mais ça donne la meme chose)
Singlecell@meta.data$Predicted_doublets = ifelse(py$DoubletScore > 0.3, "Doublet", "Singlet")
table(Singlecell@meta.data$Predicted_doublets)

# Doublet Singlet 
#    1822   36551 = #4.9% sont considéré comme des doublets. 


# #Je ne sais pas si on prend ce qui se trouve sur la ligne verte donc j'ai aussi essayé avec un seuil à 0.29
# Singlecell@meta.data$Predicted_doublets = ifelse(py$DoubletScore > 0.29, "Doublet", "Singlet")
# table(Singlecell@meta.data$Predicted_doublets)

# # Doublet Singlet 
# #    4279   34094 #12.6% sont considérés comme doublet.

library(cowplot)
library(Seurat)

Singlecell$percent.mt = PercentageFeatureSet(Singlecell,pattern = "^mt:")
ribo.genes <- grep(pattern = "rRNA", x = rownames(x = Singlecell@assays$RNA@data), value = TRUE)

pdf("/home/ocuvier/work/conda/MData/figures/figure4.pdf")
VlnPlot(Singlecell,features = c("nFeature_RNAseq", "nCount_RNAseq", "percent.mt","nFeature_ATAC","nCount_ATAC"), group.by = "Predicted_doublets",ncol = 3)
dev.off()


# Singlecell[["percent.mt"]] <- PercentageFeatureSet(object = Singlecell, pattern = "mt:")
# Singlecell[["percent.ribo"]] <- PercentageFeatureSet(object = Singlecell, pattern = "rRNA:")


plot1 <- FeatureScatter(Singlecell, feature1 = "nCount_RNAseq", feature2 = "percent.mt")
plot2 <- FeatureScatter(Singlecell, feature1 = "nCount_RNAseq", feature2 = "nFeature_RNAseq")
pdf("/home/ocuvier/work/conda/MData/figures/figure5.pdf",height = 8,width = 14)
plot1 + plot2
dev.off()

# saveRDS(Singlecell,"/home/ocuvier/work/conda/MData/save_data/Singlecell_doublet.RDATA")
Singlecell = readRDS("/home/ocuvier/work/conda/MData/save_data/Singlecell_doublet.RDATA")

Singlecell$sample[which(grepl(colnames(Singlecell),pattern="-1"))] = "Luc-KD_1"
Singlecell$sample[which(grepl(colnames(Singlecell),pattern="-2"))] = "Luc-KD_2"
Singlecell$sample[which(grepl(colnames(Singlecell),pattern="-3"))] = "Beaf32-KD_1"
Singlecell$sample[which(grepl(colnames(Singlecell),pattern="-4"))] = "Beaf32-KD_2"

table(Singlecell$sample[which(Singlecell$percent.mt<30)])
counts_RNA <- Singlecell@assays$RNA@counts

seqlevelsStyle(TAD_Regions) = "ensembl" #Pour que la version UCSC soit au 

# sous laforme uc sc on change caavec seqlevel au format ensembl pour faire la 
#recherche des sites beafs32 au niv des tads"
head(TAD_Regions)
DefaultAssay(Singlecell) = "ATAC"

#Version Ensembl :
beaf32_peaks_bed = rtracklayer::import.bed("/work/user/nschickele/ANALYSES_SINGLECELL/DATA/CHIPSEQ/BED/Beaf32_trimmed_filt_sort_summits.bed")
head(beaf32_peaks_bed)

beaf32_peaks_bed = GenomeInfoDb::keepStandardChromosomes(beaf32_peaks_bed, 
species = "Drosophila_melanogaster","coarse") #Enleve les chromosomes 
annexes, les régions pas sur les chromosomes de 2L 3L 4X 2R 3R Y.
# Je dois Comparer les TAD_Regions aux peaks de beaf32 pour trouver les gènes 
#dans un rayon de 5 à 10kbpour trouver les gènes qui sont en bordure TAD.
#d'abord on resize les regions tad, on les veut chevauchantes. Et on va séparer 
#les starts et les ends.
#resize(x, width, fix="start", use.names = TRUE, ignore.strand = FALSE) avec x = notre objet GRange et width la taille.

TAD_Regions_Start = resize(TAD_Regions, 1,fix="start", use.names = TRUE, ignore.strand = FALSE)+500
TAD_Regions_End = resize(TAD_Regions, 1,fix="end", use.names = TRUE, ignore.strand = FALSE)+500 #ça fait 500 de chaque coté donc 1000 en tout.
#rezise = je veux pus ou moins autour de ce tad

Liste_GR_TAD_Regions = c(TAD_Regions_Start, TAD_Regions_End)
Liste_TAD_Borders= GenomicRanges::reduce(Liste_GR_TAD_Regions)

#finOverlaps, il va trouver les overlaps entre beaf et tad.
#findOverlaps(query, subject) = le query on va garder ces informations et tout, donc notre query c'est beaf32

Overlaps_Beaf_TAD = subsetByOverlaps(beaf32_peaks_bed,Liste_TAD_Borders)
Overlaps_TAD_Beaf = subsetByOverlaps(Liste_TAD_Borders,beaf32_peaks_bed)

# pareil mais inverse
DefaultAssay(Singlecell) = "ATAC"
#On établi la liste des gènes
Liste_gene= Annotation(Singlecell)
Liste_gene= Liste_gene[Liste_gene$type == "gene"]
# toute l'annotation du genome (17807 genes qui font partie du genomevia flybase)
DefaultAssay(Singlecell) = "SCT"
table(Liste_gene$gene_name %in% rownames(Singlecell))
FALSE TRUE
17
7791 10016
# 10016 gene au moins exprimé dans une cellule
Gene_Singlecell = Liste_gene[Liste_gene$gene_name %in% rownames(Singlecell)] 

#On sélectionne ceux qui sont présent dans au moins une cellule.
#resize les TAD overlapés avec des sites Beaf32 pour chercher les gènes avoisinnants les TAD.

TAD_Beaf_Resize = resize(Overlaps_TAD_Beaf, 1,fix="center", use.names = TRUE, ignore.strand = FALSE)+5000
# 10000 pb autour des bordures tads
# reduce c'est pour enlever les regions qui s'overlap entrelles (qui se surperposent ) 

TAD_Beaf_Resize = GenomicRanges::reduce(TAD_Beaf_Resize)
#GRanges object with 1436 ranges and 0 metadata columns:
# on recup promo ceux auxquelson s'attend que le pic atac interagie

Promoteur_gene = promoters(Gene_Singlecell,upstream=500, downstream= 500) #On trouve les promoteurs des gènes à peu près
# On cherche là où se recoupe la position des gènes detectés avec les beaf32 qui ont été autour de régions TAD : ça nous permettra de voir les gènes 
#autour des régons TAD.

Overlaps_Gene_TAD = subsetByOverlaps(Promoteur_gene,TAD_Beaf_Resize)
head(Overlaps_Gene_TAD$gene_name)
#[1] "DhpD" "CG14641" "abs" "Gel" "Vps24" "MP1"
# Peaks to gene frontiere TAD 
#----------------------------------------------------------------------------------------

Singlecell <- LinkPeaks(
object = Singlecell,
peak.assay = "ATAC",
expression.assay = "SCT",
genes.use = c("DhpD" ,"CG14641" ,"abs" ,"Gel" ,"Vps24" ,"MP1" ))
feature_bordureTAD = c("DhpD" , "CG14641", "abs" , "Gel" , "Vps24" , "MP1")

p1 <- CoveragePlot(
object = Singlecell,
region = "DhpD",
features = "DhpD",
expression.assay = "SCT",
group.by = "sample",
extend.upstream = 1000,
extend.downstream = 2000
)

## + Coverage plot pour tous les genes, ça donne rien ...
## On test autre chose
18
# Peaks to gene frontiere TAD - Trier par z score 
#---------------------------------------------------------------------

LinkPeak_Singlecell_Gene_Present_TAD = readRDS("/home/ocuvier/work/conda/MData/save_data/LinkPeak_Singlecell_Gene_Present_TAD.RDATA")
# trier par z score
top_n(as.data.frame(Links(LinkPeak_Singlecell_Gene_Present_TAD,assay="ATAC")),n=20,wt=zscore)
DefaultAssay(Singlecell) = "ATAC"
Annotation(Singlecell)$tx_id = Annotation(Singlecell)$transcript_id
# genes : "CG3036", "nAChRalpha6","CG30456","Drs","Ids","CG10863","ImpL2","form3","CG32354"

p1 <- CoveragePlot(
object = Singlecell,
region = "CG3036",
features = "CG3036",
expression.assay = "SCT",
group.by = "sample",
extend.upstream = 5000,
extend.downstream = 10000
)

## + Coverage plot pour tous les genes, ça donne rien ...
##############################################################################################
# Verification de la depletion de BEAF32 - BULK
##############################################################################################
# EST CE QUE la depletion de Beaf32 a bien fonctionné ?
#'BEAF-32'%in%rownames(Singlecell@assays$SCT@data)

pdf("/home/ocuvier/work/conda/eData/figures/figure27.pdf", width = 20, height = 5)
FeaturePlot(Singlecell,
features = "BEAF-32",
split.by = "sample",
reduction = "umap")
dev.off()
pdf("/home/ocuvier/work/conda/eData/figures/figure28.pdf", width = 10, height = 5)
RidgePlot(Singlecell, features = "BEAF-32", group.by = "sample", ncol = 1)
dev.off()
pdf("/home/ocuvier/work/conda/eData/figures/figure29.pdf")

DotPlot(Singlecell, features = "BEAF-32", group.by = "sample") + RotatedAxis()
dev.off()

# LA REPONSE EST OUI
# la deplation de beaf32 a bien marché on voit bien sur la fig 29 que beaf kd1 et d2
# ne ssont pas exprimé
# remarque :
# le gene meme si il est bien exprimé on peut ne pas le detecter
#de base les niveaux d'expression de beaf32 n'est pasenorme comparé aux 
#features Singlecell. Apres depletion (les transcrit sont comme en 
#"competition", du aunombrede read limité; partagé avec les fragments atac et 
#les transcrit (les genes biens exprimé sont facilements sequencé et detecté. 
#donc en moyenne le nfeatures qu'on voit estime a env 10000 gene exprimé par 
#cellule alors qu'on en detecte que 3 a 4 milles, du au fait que les genes de 
#menage soient ultra abondant (c 'est pour ca que pour faire ca bien on 
#sequence bcp bcp de cellule, ou alors on reduit le nombre de cellule et on 
#augmente le nombre de reads')))

##############################################################################################
# Expression Differentielle - PSEUDOBULK
##############################################################################################


Singlecell$condition = str_remove_all(Singlecell$sample,pattern ="_[1|2]")
# Ordonner avec log2FC décroissant 
#--------------------------------------------------------------------------

# RNA, Identification des gènes différentiellement exprimés
DefaultAssay(Singlecell) = "SCT"
pseudo_rna_differential = FindMarkers(Singlecell, group.by = "condition", ident.1="Beaf32-KD", ident.2="Luc-KD")
pseudo_rna_differential_avglog2FC = pseudo_rna_differential %>% dplyr::filter(avg_log2FC> 0.3 & p_val_adj <0.05) %>% dplyr::arrange(desc(avg_log2FC))
pseudo_rna_differential_avglog2FC_names = rownames(pseudo_rna_differential %>% dplyr::filter(avg_log2FC> 0.3 & p_val_adj <0.05)%>% dplyr::arrange(desc(avg_log2FC)))
head(pseudo_rna_differential_avglog2FC_names)

#[1] "CG7299" "CG16713" "CG10702" "CG32037" "CG42807" "Cyp12a4"
length(unique(pseudo_rna_differential_avglog2FC_names))

# [1] 94
#ATAC, Évaluation de la variabilité du signal ATAC
#20

DefaultAssay(Singlecell) = "ATAC"
pseudo_atac_differential = FindMarkers(Singlecell, group.by = "condition", ident.1="Beaf32-KD", ident.2="Luc-KD")
pseudo_atac_differential_avglog2FC = pseudo_atac_differential %>% dplyr::filter(avg_log2FC> 0.3 & p_val_adj <0.05) %>% dplyr::arrange(desc(avg_log2FC))
pseudo_atac_differential_avglog2FC_names = rownames(pseudo_atac_differential %>% dplyr::filter(avg_log2FC> 0.3 & p_val_adj <0.05)%>% dplyr::arrange(desc(avg_log2FC)))
head(pseudo_atac_differential_avglog2FC_names)
#[1] "2R-14771639-14772572" "X-13713303-13714139" "2L-19069642-19070402"
#[4] "X-8554059-8554940" "X-4129135-4129967" "3R-15461777-15462655"

length(unique(pseudo_atac_differential_avglog2FC_names))
# [1] 32

# Voir Tableau :
pseudo_atac_differential %>% dplyr::filter(avg_log2FC> 0.3 & p_val_adj <0.05) %>% dplyr::arrange(desc(avg_log2FC))
# seuil plus severe:
# pseudo_rna_differential %>% dplyr::filter(avg_log2FC> 0.3 & p_val_adj <0.05 & pct.1 > 0.2)%>% dplyr::arrange(desc(avg_log2FC))
#on garde pas c'est trop severe
# LinkPeaks 
#-----------------------------------------------------------------------------------------
# on refait linkpeaks avec les genes trouvés precedement

Singlecell <- LinkPeaks(
object = Singlecell,
peak.assay = "ATAC",
expression.assay = "SCT",
genes.use = pseudo_rna_differential_avglog2FC_names
)
# Testing 94 genes and 14727 peaks
# Évaluation de la corrélation pour examiner la relation entre le signal ATAC et l'expression génique
unique(Singlecell@assays$ATAC@links$gene[which(Singlecell@assays$ATAC@links$peak%in%rownames(p & Singlecell@assays$ATAC@links$gene%in%rownames(pseudo_rna_differential_avglog2FC))])
unique(Singlecell@assays$ATAC@links$gene[which(Singlecell@assays$ATAC@links$peak%in%rownames & Singlecell@assays$ATAC@links$gene%in%rownames(pseudo_rna_differential_avglog2FC))])

# [1] "CG10702" "AdamTS-A" "CG43759"
21
# "CG10702" "AdamTS-A" "CG43759" = gene differencielements exprimés significatifs avec des pics significatifs

# Visualisation 
#---------------------------------------------------------------------------------
# [1] "CG10702" "AdamTS-A" "CG43759"
# Expression differentielle par rapport aux cluster trouvés precedement

DefaultAssay(Singlecell) = "SCT"
featureRNA <- c("CG10702" , "AdamTS-A" ,"CG43759")
pdf("/home/ocuvier/work/conda/eData/figures/figureexpressionn_featureRNA.pdf", width = 20, height = 15)

FeaturePlot(Singlecell,
features = featureRNA,
reduction = "umap")
dev.off()

# Ca donne rien
featureRNA <- c("CG10702" , "AdamTS-A" ,"CG43759")
pdf("/home/ocuvier/work/conda/eData/figures/figuregene_cluster.pdf", width = 15, height = 10)
DotPlot(Singlecell, features = featureRNA) + RotatedAxis()
dev.off()

# "CG43759" semble plus s'exprimer dans le cluster 2
# CoveragePlot - "CG10702" "AdamTS-A" "CG43759"
DefaultAssay(Singlecell) = "ATAC"
p1 <- CoveragePlot(
object = Singlecell,
region = "CG10702" ,
features = "CG10702" ,
expression.assay = "SCT",
group.by = "condition",
extend.upstream = 1000,
extend.downstream = 1000
)

p2 <- CoveragePlot(
object = Singlecell,
region = "AdamTS-A" ,
features = "AdamTS-A" ,
expression.assay = "SCT",
group.by = "condition",
extend.upstream = 10000,
extend.downstream = 10000
)

p3 <- CoveragePlot(
object = Singlecell,
region = "CG43759" ,
features = "CG43759" ,
expression.assay = "SCT",
group.by = "condition",
extend.upstream = 10000,
extend.downstream = 10000
)

pdf("/home/ocuvier/work/conda/eData/figures/figurecoveragePlot.pdf", width = 10, height = 30)
patchwork::wrap_plots(p1,p2,p3, ncol = 1)
dev.off()

# Expression plot- "CG10702" "AdamTS-A" "CG43759"
p1bis <- ExpressionPlot(
object = Singlecell,
features = "CG10702",
assay = "SCT",
group.by = "condition"
)

p2bis <- ExpressionPlot(
object = Singlecell,
features = "AdamTS-A",
assay = "SCT",
group.by = "condition"
)

p3bis <- ExpressionPlot(
object = Singlecell,
features = "CG43759",
assay = "SCT",
group.by = "condition"
)

pdf("/home/ocuvier/work/conda/eData/figures/ExpressionPlot.pdf", width = 15, height = 5)
patchwork::wrap_plots(p1bis,p2bis,p3bis, ncol = 3)
dev.off()
# Comptage du nombre de cellule par condition beaf32KD & Luc KD 
#------------------------------------

# matrice vide
cluster_means <- matrix(nrow = 14, ncol = 2)
for (i in 1:14) {
cluster = table(Singlecell$condition[Singlecell@meta.data$seurat_clusters == i-1])

mean_beaf32 = cluster["Beaf32-KD"]
mean_luc = cluster["Luc-KD"]
somme = mean_beaf32+ mean_luc
cluster_means[i, 1] <- paste0(round(mean_beaf32 * 100 / somme), " %")
cluster_means[i, 2] <- paste0(round(mean_luc * 100 / somme), " %")
}

rownames(cluster_means) = c("0", "1" ,"2", "3", "4" ,"5" ,"6" ,"7", "8", "9", "10", "11", "12" ,"13")
colnames(cluster_means) = c("Beaf32-KD","Luc-KD")
saveRDS(Singlecell,"/home/ocuvier/work/conda/eData/save_data/Singlecell_dataTest.RDATA")
Singlecell = readRDS("/home/ocuvier/work/conda/eData/save_data/Singlecell_dataTest.RDATA")
save.image("17052024_sc.RDATA")
load("17052024_sc.RDATA")

##############################################################################################
# Visualisation de la depletion de BEAF32 - PSEUDOBULK
##############################################################################################

pdf("/home/ocuvier/work/conda/eData/figures/figureUMAP_condition3.pdf",width = 11, height = 5)
FeaturePlot(Singlecell,
features = "BEAF-32",
split.by = "condition",
reduction = "umap")
dev.off()

p3 <- FeaturePlot(Singlecell,
features = "BEAF-32",
split.by = "condition",
reduction = "umap.rna")

p4 <- FeaturePlot(Singlecell,
features = "BEAF-32",
split.by = "condition",
reduction = "umap.atac")

pdf("/home/ocuvier/work/conda/eData/figures/figureUMAP_condition.pdf", width = 11, height = 5)
patchwork::wrap_plots(p3,ncol = 1)
patchwork::wrap_plots(p4,ncol = 1)
dev.off()

pdf("/home/ocuvier/work/conda/eData/figures/figureUMAP_condition2.pdf", width = 22, height = 5)
patchwork::wrap_plots(p3,p4,ncol = 2)
dev.off()

### On peut faire les limites et conclusion de nos analyses.

