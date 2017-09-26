# All Aanalyses performed with R 3.3.3

# Load libraries 
library(limma)
library(DESeq2)
library(edgeR)
library(splines)
library(RColorBrewer)
library(ggplot2)

## Set working directory
setwd("exposure")

################################################################################################
#
#           1. Bring in files for field and lab data and trim counts files  
#
#
################################################################################################

# Load lab counts file generated from stringTie
lab_all <- read.csv("lab_gene_count_matrix.csv",header=T, row.names=1)

# Count number of transcripts and number of individuals
dim(lab_all)

# Subset data with outliers removed (BONGC6, PSOGC1)
lab_all = lab_all[,-which(colnames(lab_all) %in% c("BONGC6","PSOGC1"))]

# Trim data: Each row (i.e. transcript) must have greater than 2 counts per million (cpm) and at least 3 of the individuals must have counts data
keep <- rowSums(cpm(lab_all)>2) >= 3     

# Keep data that meet the filtering criteria
lab <- lab_all[keep,]

# Count number of transcripts left after trimming
dim(lab)

# Load field counts file generated from stringTie
field_all <- read.csv("field_gene_count_matrix.csv",header=T, row.names=1)

# Count number of transcripts and number of individuals
dim(field_all)

# Subset data wiht outlier removed (MX71)
field_all = field_all[,-which(colnames(field_all) %in% c("MX71"))]

# Trim data: Each row (i.e. transcript) must have greater than 2 counts per million (cpm) and at least 3 of the individuals must have counts data
keep <- rowSums(cpm(field_all)>2) >= 3     

# Keep data that meet the filtering criteria
field <- field_all[keep,]

# Count number of transcripts left after trimming
dim(field)


#####################################################################################################
#
#   2. Identify candidate gene lists for the two drainages (Lab and Field separately) 
#           
#
#####################################################################################################
#### Lab ####
# Group lab by exposure treatment and population
group = c(rep("TacoFreshControl",5), rep("TacoFreshHigh",5), rep("TacoFreshLow",6), rep("TacoFreshMedium",6),
          rep("PuyaFreshControl",6), rep("PuyaFreshHigh",3), rep("PuyaFreshLow",6), rep("PuyaFreshMedium",5),
          rep("PuyaSulfidicControl",6), rep("PuyaSulfidicHigh",3), rep("PuyaSulfidicLow",6), rep("PuyaSulfidicMedium",6),
          rep("TacoSulfidicControl",5), rep("TacoSulfidicHigh",5), rep("TacoSulfidicLow",6), rep("TacoSulfidicMedium",6))

# Generate dataframe with lab sample names and group 
samples = data.frame(cbind(colnames(lab)), as.character(group))
colnames(samples) = c("samples", "group")

# Create a DGEList object to hold the dataset 
m = DGEList(counts = lab, group = samples$group)

# Calculate normalized factors based on raw library sizes
m = calcNormFactors(m)

# Create a design marix 
design <- model.matrix(~0+group)

# Add column names based on sample names in the lab group
colnames(design) <- levels(factor(samples$group))

# Estimate common dispersion and tagwise dispersion 
m <- estimateDisp(m, design)

# Given tagwise dispersion and a design matrix, glmFIT fits the negative binomial GLM for each tag. Will then produce a DGEGLM object with new components. 
fit <- glmFit(m,design)

# Construct contrast matrix of comparisons (compare each exposure treatment to the control treatment for that drainage)
my.contrasts <- makeContrasts(
  Puya_low    = PuyaSulfidicLow - PuyaFreshControl,
  Puya_medium = PuyaSulfidicMedium - PuyaFreshControl,
  Puya_high   = PuyaSulfidicHigh - PuyaFreshControl,
  Taco_low    = TacoSulfidicLow - TacoFreshControl,
  Taco_medium = TacoSulfidicMedium - TacoFreshControl,
  Taco_high   = TacoSulfidicHigh - TacoFreshControl,
  levels = design
)

# Use the DGEGLM object to perform likelihood ratio test based on the contrast matrix of comparisons
PL_lrt = glmLRT(fit, contrast = my.contrasts[,"Puya_low"])
PM_lrt = glmLRT(fit, contrast = my.contrasts[,"Puya_medium"])
PH_lrt = glmLRT(fit, contrast = my.contrasts[,"Puya_high"])
TL_lrt = glmLRT(fit, contrast = my.contrasts[,"Taco_low"])
TM_lrt = glmLRT(fit, contrast = my.contrasts[,"Taco_medium"])
TH_lrt = glmLRT(fit, contrast = my.contrasts[,"Taco_high"])

# Summary of differential expression that was up, down or not significant in each comparison 
summary(decideTestsDGE(PL_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(PM_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(PH_lrt, adjust.method="BH", p.value = 0.05))

summary(decideTestsDGE(TL_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(TM_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(TH_lrt, adjust.method="BH", p.value = 0.05))

# Extract top expressed genes (based on significant upregulated and downregulated genes)
PL = topTags(PL_lrt, n = 473)
PM = topTags(PM_lrt, n = 284)
PH = topTags(PH_lrt, n = 1298)
TL = topTags(TL_lrt, n = 138)
TM = topTags(TM_lrt, n = 1058)
TH = topTags(TH_lrt, n = 71)

# Write logFC, logCPM, LR, Pvalue and FDR statistics for each significant gene to a .csv file 
write.csv(PL$table, file = "PL.csv")
write.csv(PM$table, file = "PM.csv")
write.csv(PH$table, file = "PH.csv")
write.csv(TL$table, file = "TL.csv")
write.csv(TM$table, file = "TM.csv")
write.csv(TH$table, file = "TH.csv")


# Field
# Group field by population
field.group = c(rep("PuyaFieldSulfidic", 5), rep("PuyaFieldFresh",6), rep("TacoFieldFresh",6), rep("TacoFieldSulfidic", 5))

# Generate dataframe with lab sample names and group 
field.samples = data.frame(cbind(colnames(field)), as.character(field.group))
colnames(field.samples) = c("samples", "group")

# Create a DGEList object to hold the dataset 
field.m = DGEList(counts = field, group = field.samples$group)

# Calculate normalized factors based on raw library sizes
field.m = calcNormFactors(field.m)

# Create a design marix 
field.design <- model.matrix(~0+field.group)

# Add column names based on sample names in the lab group
colnames(field.design) <- levels(factor(field.samples$group))

# Estimate common dispersion and tagwise dispersion 
field.m <- estimateDisp(field.m, field.design)

# Given tagwise dispersion and a design matrix, glmFIT fits the negative binomial GLM for each tag. Will then produce a DGEGLM object with new components. 
field.fit <- glmFit(field.m,field.design)

# Construct contrast matrix of comparisons (compare sulfide adapted to nonsulfide adapted populations in each drainage)
my.contrasts <- makeContrasts(
  Field_Puya = PuyaFieldSulfidic - PuyaFieldFresh, 
  Field_Taco = TacoFieldSulfidic - TacoFieldFresh,
  levels = field.design
)

# Use the DGEGLM object to perform likelihood ratio test based on the contrast matrix of comparisons
FieldP_lrt = glmLRT(field.fit, contrast = my.contrasts[,"Field_Puya"])
FieldT_lrt = glmLRT(field.fit, contrast = my.contrasts[,"Field_Taco"])

# Summary of differential expression that was up, down or not significant in each comparison 
summary(decideTestsDGE(FieldP_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(FieldT_lrt, adjust.method="BH", p.value = 0.05))

# Extract top expressed genes (based on significant upregulated and downregulated genes)
FieldP = topTags(FieldP_lrt, n = 1793)
FieldT = topTags(FieldT_lrt, n = 2462)

# Write logFC, logCPM, LR, Pvalue and FDR statistics for each significant gene to a .csv file 
write.csv(FieldP$table, file = "FieldP.csv")
write.csv(FieldT$table, file = "FieldT.csv")

#####################################################################################################
#
#   2A. Extract up- and downregulated genes for each constrasting comparison in lab and field separately 
#           
#
#####################################################################################################
# Puyacatengo Lab
PL.up   = rownames(PL[PL$table$logFC > 0,])
PL.down = rownames(PL[PL$table$logFC < 0,])
PM.up   = rownames(PM[PM$table$logFC > 0,])
PM.down = rownames(PM[PM$table$logFC < 0,])
PH.up   = rownames(PH[PH$table$logFC > 0,])
PH.down = rownames(PH[PH$table$logFC < 0,])

# Tacotalpa Lab
TL.up   = rownames(TL[TL$table$logFC > 0,])
TL.down = rownames(TL[TL$table$logFC < 0,])
TM.up   = rownames(TM[TM$table$logFC > 0,])
TM.down = rownames(TM[TM$table$logFC < 0,])
TH.up   = rownames(TH[TH$table$logFC > 0,])
TH.down = rownames(TH[TH$table$logFC < 0,])

# Puyacatengo Field
FieldP.up    = rownames(FieldP[FieldP$table$logFC > 0,])
FieldP.down  = rownames(FieldP[FieldP$table$logFC < 0,])

# Puyacatengo Lab 
FieldT.up    = rownames(FieldT[FieldT$table$logFC > 0,])
FieldT.down  = rownames(FieldT[FieldT$table$logFC < 0,])


###################################VENN DIAGRAMS#############################################################
###################################FOUR WAY VENN########################################################

# Create unique variable, and construct venns 
# Upregulated
universe.UP <- unique(c(PL.up, PM.up, PH.up, TL.up, TM.up, TH.up, FieldP.up, FieldT.up))
GroupA <- universe.UP %in% c(TL.up, TM.up, TH.up)
GroupB <- universe.UP %in% c(PL.up, PM.up, PH.up)
GroupC <- universe.UP %in% FieldT.up
GroupD <- universe.UP %in% FieldP.up
input.df <- data.frame(TacoLab=GroupA, PuyaLab=GroupB, TacoField=GroupC, PuyaField=GroupD)
colnames(input.df)=c("TacoL","PuyaL", "TacoF", "PuyaF")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)
Q1all.up <- universe.UP[which(input.df["TacoL"] == T & input.df["PuyaL"] == T & input.df["TacoF"] == T & input.df["PuyaF"] == T)]

#Downregulated
universe.DOWN <- unique(c(PL.down, PM.down, PH.down, TL.down, TM.down, TH.down, FieldP.down, FieldT.down))
GroupA <- universe.DOWN %in% c(TL.down, TM.down, TH.down)
GroupB <- universe.DOWN %in% c(PL.down, PM.down, PH.down)
GroupC <- universe.DOWN %in% FieldT.down
GroupD <- universe.DOWN %in% FieldP.down
input.df <- data.frame(TacoLab=GroupA, PuyaLab=GroupB, TacoField=GroupC, PuyaField=GroupD)
colnames(input.df)=c("TacoL","PuyaL", "TacoF", "PuyaF")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)
Q1all.down <- universe.DOWN[which(input.df["TacoL"] == T & input.df["PuyaL"] == T & input.df["TacoF"] == T & input.df["PuyaF"] == T)]

############################## TWO-WAY VENN DIAGRAMS##########################################
# Create unique variable, and construct venns 
# Tacotalpa Upregulated
PSOBONall.up = unique(c(TL.up, TM.up, TH.up))
universe.UP <- unique(c(PSOBONall.up, FieldT.up))
GroupA <- universe.UP %in% PSOBONall.up
GroupB <- universe.UP %in% FieldT.up
input.df <- data.frame(LT=GroupA, FT=GroupB)
colnames(input.df)=c("Lab-Taco","Field-Taco")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)
#create variable for upregulated candidate genes in the Tacotalpa (shared between lab and field)
Q1alltaco.UP <- universe.UP[which(input.df["Lab-Taco"] == T & input.df["Field-Taco"] == T)]

#Puyacatengo Upregulated
LLUIXAall.up = unique(c(PL.up, PM.up, PH.up))
universe.UP <- unique(c(LLUIXAall.up, FieldP.up))
GroupC <- universe.UP %in% LLUIXAall.up
GroupD <- universe.UP %in% FieldP.up
input.df <- data.frame(LP=GroupC, FP=GroupD)
colnames(input.df)=c("Lab-Puya","Field-Puya")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)
# create variable for upregulated candidate genes in the Puyacatengo (shared between lab and field)
Q1allpuya.UP <- universe.UP[which(input.df["Lab-Puya"] == T & input.df["Field-Puya"] == T)]

# Tacotalpa Downregulated
PSOBONall.down = unique(c(TL.down, TM.down, TH.down))
universe.DOWN <- unique(c(PSOBONall.down, FieldT.down))
GroupA <- universe.DOWN %in% PSOBONall.down
GroupB <- universe.DOWN %in% FieldT.down
input.df <- data.frame(LT=GroupA, FT=GroupB)
colnames(input.df)=c("Lab-Taco","Field-Taco")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)
# create variable for downregulated candidate genes in the Tacotalpa (shared between lab and field)
Q1alltaco.DOWN <- universe.DOWN[which(input.df["Lab-Taco"] == T & input.df["Field-Taco"] == T)]

# Venn Puyacatengo DOWN
LLUIXAall.down = unique(c(PL.down, PM.down, PH.down))
universe.DOWN <- unique(c(LLUIXAall.down, FieldP.down))
GroupC <- universe.DOWN %in% LLUIXAall.down
GroupD <- universe.DOWN %in% FieldP.down
input.df <- data.frame(LP=GroupC, FP=GroupD)
colnames(input.df)=c("Lab-Puya","Field-Puya")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)
# create variable for downregulated candidate genes in the Puyacatengo (shared between lab and field)
Q1allpuya.DOWN <- universe.DOWN[which(input.df["Lab-Puya"] == T & input.df["Field-Puya"] == T)]

#####################################################################################################
#
#           2B. MDS plots 
#           
#
#####################################################################################################
################################ ALL LAB AND FIELD ###################################
# Note this corresponds to Figure S1 in Passow et al. In revision
# Combine lab and field data
total <- cbind(lab_all,field_all)

# Group replicates
group <-c(rep("lab",85),rep("field",22))

# Trim out any transcripts with no counts
total <- total[rowSums(total) > 0, ]

# Count number of transcripts and individuals 
dim(total)

# Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=total, group=group)

# Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

# Construct MDS plot on top 10,000 expressed transcripts 
mds <- plotMDS(dg, top=10000, gene.selection = "common")

# Convert "names" into symbols and change color to represent the different organs
# Blue corresponds to lab data and green corresponds to field data
plot(mds, col=c(rep("blue",85),rep("green",22)), pch=16)

# Extract residuals
datamdsall <- (mds$cmdscale.out)

# Save residuals into .csv file to clean up MDS plots in Excel
write.csv(datamdsall, file="all_data_mds.csv")

################################ ALL LAB ###################################
# Note this corresponds to Figure S2A-B in Passow et al. In revision

# Organize the samples into their respective groups 
d <- subset(lab_all, select = c(BONGC1, BONGC2, BONGC3, BONGC4, BONGC5,
                                BONGL1, BONGL2, BONGL3, BONGL4, BONGL5, BONGL6,
                                BONGM1, BONGM2, BONGM3, BONGM4, BONGM5, BONGM6,
                                BONGH1, BONGH2, BONGH3, BONGH4, BONGH5,
                                PSOGC2, PSOGC3, PSOGC4, PSOGC5, PSOGC6,
                                PSOGL1, PSOGL2, PSOGL3, PSOGL4, PSOGL5, PSOGL6,
                                PSOGM1, PSOGM2, PSOGM3, PSOGM4, PSOGM5, PSOGM6,
                                PSOGH1, PSOGH2, PSOGH3, PSOGH5, PSOGH6,
                                IXAGC1, IXAGC2, IXAGC3, IXAGC4, IXAGC5, IXAGC6,
                                IXAGL1, IXAGL2, IXAGL3, IXAGL4, IXAGL5, IXAGL6,
                                IXAGM1, IXAGM2, IXAGM3, IXAGM5, IXAGM6,
                                IXAGH1, IXAGH2, IXAGH3,
                                LLUGC1, LLUGC2, LLUGC3, LLUGC4, LLUGC5, LLUGC6,
                                LLUGL1, LLUGL2, LLUGL3, LLUGL4, LLUGL5, LLUGL6,
                                LLUGM1, LLUGM2, LLUGM3, LLUGM4, LLUGM5, LLUGM6,
                                LLUGH2, LLUGH5, LLUGH6))

# Group repliates
group <-c(rep("TNonsulfidic",22),rep("Tsulfidic",22),rep("PNonsulfidic",20),rep("Psulfidic",21))

# Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

# Count number of transcripts and individuals
dim(d)

# Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

# Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

# Construct MDS plot on top 10,000 expressed transcripts 
mds <- plotMDS(dg, top=10000, gene.selection = "common")

# Convert "names" into symbols and change color to represent the different habitats
# Blue corresponds to control, yellow is 0.5mM, orange is 3mM and red is 6mM
# Circles are Tacotalpa nonsulfidic, squates are Tacotalpa sulfidic, triangles are Puyacatengo nonsulfidic and diamonds are Puyacatengo sulfidic
# Note that due to the difficultly of reading both drainages on one figure, we made two separate figures for each drainage based on the mds output
plot(mds, col=c(rep("blue",5),rep("yellow",6),rep("orange",6),rep("red",5),
                rep("blue",5),rep("yellow",6),rep("orange",6),rep("red",5),
                rep("blue",6),rep("yellow",6),rep("orange",5),rep("red",3),
                rep("blue",6),rep("yellow",6),rep("orange",6),rep("red",3)), pch=c(rep(16,22),rep(15,22),rep(17,20),rep(18,21)))

# Extract residuals
Lab <- (mds$cmdscale.out)

# Save residuals into .csv file to clean up MDS plots in Excel
write.csv(Lab, file="lab_mds.csv")

################################ALL FIELD###################################
# Note this corresponds to Figure S3 in Passow et al. In revision

# Organize field by population
d <- subset(field_all, select = c(MX57, MX59, MX60, MX61, MX62, MX63,
                                  MX73, MX74, MX75, MX76, MX77,
                                  MX50, MX51, MX52, MX53, MX54, MX55,
                                  MX44, MX45, MX46, MX48, MX49))

# Group repliates
group <-c(rep("BON",6),rep("PSO",5),rep("IXA",6),rep("LLU",5))

# Trim out any transcripts with no counts
d <- d[rowSums(d) > 0, ]

# Count number of transcripts and individuals
dim(d)

# Create a DGEList object to hold the dataset to be analysed in edgeR 
dgeg <- DGEList(counts=d, group=group)

# Calculate normalized factors based on raw library sizes
dg <- calcNormFactors(dgeg)

# Construct MDS plot on top 10,000 genes with the largest standard deviations between samples  
mds <- plotMDS(dg, top=10000, gene.selection = "common")

# Convert "names" into symbols and change color to represent the different habitats
# Circles correspond to the Tacotalpa drainage and Squares correspond to the Puyacatnego drainage
# Blue correspond to Tacotalpa nonsulfidic, yellow are Tacotalpa sulfidic, light blue are Puyacatengo nonsulfidic and goldenrod are Puyacatnego sulfidic
plot(mds, col=c(rep("blue",6),rep("yellow",5),rep("lightblue",6),rep("goldenrod",5)), 
     pch=c(rep(16,11),rep(15,11)))

# Extract residuals
Fieldall <- (mds$cmdscale.out)

# Save residuals into .csv file to clean up MDS plots in Excel
write.csv(Fieldall, file="Field_mds.csv")

################################################################################################
#
#           2C. Extract gene annotations for candidate genes from Table S2 BLAST output 
#           
#
################################################################################################
# Load lab BLAST outputfile (Table S2 was saved as a .csv)
BLAST <- read.csv("Table S2 - Poecilia mexicana Annotations.csv")

# subset gene IDs, subject sequence IDs and protein IDs
BLAST_annot <- subset(BLAST, select = c(gene.ID, Subject.sequence.ID, Protein.annotations, gene.name, E.value))


# All output below was used for gene annotations in Table 1 and Table S3
######################################TACOTALPA AND PUYACATENGO UPREGULATED#####################################
# make gene names from Tacotalpa candidate genes a data.frame
Q1all.UP <- data.frame(Q1all.up)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Q1all.UP) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
All_candidates.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Q1all.UP, by = "gene.ID")

######################################TACOTALPA AND PUYACATENGO DOWNREGULATED#####################################
# make gene names from Tacotalpa candidate genes a data.frame
Q1all.DOWN <- data.frame(Q1all.down)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Q1all.DOWN) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
All_candidates.down <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Q1all.DOWN, by = "gene.ID")

######################################TACOTALPA UPREGULATED#####################################
# make gene names from Tacotalpa candidate genes a data.frame
Q1alltaco1.UP <- data.frame(Q1alltaco.UP)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Q1alltaco1.UP) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Taco_candidates.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Q1alltaco1.UP, by = "gene.ID")

######################################TACOTALPA DOWNREGULATED#####################################
# make gene names from Tacotalpa candidate genes a data.frame
Q1alltaco1.DOWN <- data.frame(Q1alltaco.DOWN)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Q1alltaco1.DOWN) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Taco_candidates.down <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Q1alltaco1.DOWN, by = "gene.ID")

######################################PUYACATENGO UPREGULATED#####################################
# make gene names from Tacotalpa candidate genes a data.frame
Q1allpuya1.UP <- data.frame(Q1allpuya.UP)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Q1allpuya1.UP) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Puya_candidates.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Q1allpuya1.UP, by = "gene.ID")

######################################PUYACATENGO UPREGULATED#####################################
# make gene names from Tacotalpa candidate genes a data.frame
Q1allpuya1.DOWN <- data.frame(Q1allpuya.DOWN)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Q1allpuya1.DOWN) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Puya_candidates.down <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Q1allpuya1.DOWN, by = "gene.ID")

################################################################################################
#
#           3. Identify evidence of constitutive expression
#           
#
################################################################################################
# Set up contrast matrix for Tacotalpa lab in nonsulfidic control and sulfidic control
constiutive.contrasts <- makeContrasts( 
  Tconstit = TacoSulfidicControl - TacoFreshControl, 
  Pconstit = PuyaSulfidicControl - PuyaFreshControl, 
  levels = design
)

# Use the DGEGLM object to perform likelihood ratio test based on the contrast matrix of comparisons
TConstit_lrt = glmLRT(fit, contrast = constiutive.contrasts[,"Tconstit"])
PConstit_lrt = glmLRT(fit, contrast = constiutive.contrasts[,"Pconstit"])

# Summary of differential expression that was up, down or not significant in each comparison 
summary(decideTestsDGE(TConstit_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(PConstit_lrt, adjust.method="BH", p.value = 0.05))

# Extract up- (logFC > 0) and downregulated (logFC < 0) out separately for each drainage based on FDR cut off of 0.05
Tconst = topTags(TConstit_lrt, n = sum(summary(decideTestsDGE(TConstit_lrt, p=0.05, adjust="BH"))[c(1,3),]))
Pconst = topTags(PConstit_lrt, n = sum(summary(decideTestsDGE(PConstit_lrt, p=0.05, adjust="BH"))[c(1,3),]))

# Separate up- (logFC > 0) and downregulated (logFC < 0) out into two variables for each drainage
# Tacotalpa
Tconst.up   = rownames(Tconst[Tconst$table$logFC > 0,])
Tconst.down = rownames(Tconst[Tconst$table$logFC < 0,])

# Puyacatengo
Pconst.up   = rownames(Pconst[Pconst$table$logFC > 0,])
Pconst.down = rownames(Pconst[Pconst$table$logFC < 0,])

# Rename variables for Tacotalpa
Q2Tacolab.UP = Tconst.up
Q2Tacolab.DOWN = Tconst.down

### Intersect candidate genes (OBJECTIVE 1) with genes with evidence of consititutive expression in Tacotalpa
# UPREGULATED
Q1Q2taco.UP <- intersect(Q1alltaco.UP, Q2Tacolab.UP)
# DOWNREGULATED
Q1Q2taco.DOWN <- intersect(Q1alltaco.DOWN, Q2Tacolab.DOWN)

# Rename variables for Puyacatengo
Q2Puyalab.UP = Pconst.up
Q2Puyalab.DOWN = Pconst.down

### Intersect candidate genes (OBJECTIVE 1) with genes with evidence of consititutive expression in Puyacatengo
# UPREGULATED
Q1Q2Puya.UP <- intersect(Q1allpuya.UP, Q2Puyalab.UP)
# DOWNREGULATED
Q1Q2Puya.DOWN <- intersect(Q1allpuya.DOWN, Q2Puyalab.DOWN)

#################################VENN DIAGRAMS###############################################################
#################################TWO WAY VENNS###############################################################
# Create unique variable, and construct venns 
universe.UP <- unique(c(Q1Q2taco.UP, Q1Q2Puya.UP))
GroupA <- universe.UP %in% Q1Q2taco.UP
GroupB <- universe.UP %in% Q1Q2Puya.UP
input.df <- data.frame(LT=GroupA, FT=GroupB)
colnames(input.df)=c("Lab-Taco","Lab-Puya")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)
# create variable for genes with evidence of constiutituve expression that are upregulated in Tacotalpa and Puyacatengo
Evolvedtacopuya.up <- universe.UP[which(input.df["Lab-Taco"] == T & input.df["Lab-Puya"] == T)]
Evolvedtaco.up1 <- universe.UP[which(input.df["Lab-Taco"] == T & input.df["Lab-Puya"] == F)]
Evolvedpuya.up1 <- universe.UP[which(input.df["Lab-Taco"] == F & input.df["Lab-Puya"] == T)]

## Combine genes in with evidence of constitutive expression in the Tacotalpa that are shared and unique 
Evolvedtaco.up <- unique(c(Evolvedtacopuya.up, Evolvedtaco.up1))

## Combine genes in with evidence of constitutive expression in the Puyacatengo that are shared and unique 
Evolvedpuya.up <- unique(c(Evolvedtacopuya.up, Evolvedpuya.up1))

# Create unique variable, and construct venns 
universe.DOWN <- unique(c(Q1Q2taco.DOWN, Q1Q2Puya.DOWN))
GroupC <- universe.DOWN %in% Q1Q2taco.DOWN
GroupD <- universe.DOWN %in% Q1Q2Puya.DOWN
input.df <- data.frame(LP=GroupC, FP=GroupD)
colnames(input.df)=c("Lab-Taco","Lab-Puya")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)
# create variable for genes with evidence of constiutituve expression that are downregulated in Tacotalpa and Puyacatengo
Evolvedtacopuya.down <- universe.DOWN[which(input.df["Lab-Taco"] == T & input.df["Lab-Puya"] == T)]
Evolvedtaco.down1 <- universe.DOWN[which(input.df["Lab-Taco"] == T & input.df["Lab-Puya"] == F)]
Evolvedpuya.down1 <- universe.DOWN[which(input.df["Lab-Taco"] == F & input.df["Lab-Puya"] == T)]

## Combine genes in with evidence of constitutive expression in the Tacotalpa that are shared and unique 
Evolvedtaco.down <- unique(c(Evolvedtacopuya.down, Evolvedtaco.down1))

## Combine genes in with evidence of constitutive expression in the Puyacatengo that are shared and unique 
Evolvedpuya.down <- unique(c(Evolvedtacopuya.down, Evolvedpuya.down1))

################################################################################################
#
#           3A. Extract gene annotations for constitutive expression from Table S2 BLAST output 
#           
#
################################################################################################

# All output below can be found in Table S7
######################################TACOTALPA UPREGULATED#####################################
# make gene names from Tacotalpa candidate genes a data.frame
Evolvedtaco1.up <- data.frame(Evolvedtaco.up)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Evolvedtaco1.up) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Taco_constit.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Evolvedtaco1.up, by = "gene.ID")

######################################TACOTALPA DOWNREGULATED#####################################
# make gene names from Tacotalpa candidate genes a data.frame
Evolvedtaco1.down <- data.frame(Evolvedtaco.down)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Evolvedtaco1.down) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Taco_constit.down <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Evolvedtaco1.down, by = "gene.ID")

######################################PUYACATENGO UPREGULATED#####################################
# make gene names from Tacotalpa candidate genes a data.frame
Evolvedpuya1.up <- data.frame(Evolvedpuya.up)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Evolvedpuya1.up) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Puya_constit.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Evolvedpuya1.up, by = "gene.ID")

######################################PUYACATENGO DOWNREGULATED#####################################
# make gene names from Tacotalpa candidate genes a data.frame
Evolvedpuya1.down <- data.frame(Evolvedpuya.down)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Evolvedpuya1.down) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Puya_constit.down <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Evolvedpuya1.down, by = "gene.ID")

################################################################################################
#
#           4. Identify evidence of plasticity (ancestral, adaptive and canalization)
#           
#
################################################################################################
# Set up contrast matrix for each drainage in the lab with comparisons between exposure and control 
# Puyacatengo and Tacotalpa
plasticity.contrasts <- makeContrasts(
  PuyaFresh_controllow    = PuyaFreshLow - PuyaFreshControl,
  PuyaFresh_controlmedium = PuyaFreshMedium - PuyaFreshControl,
  PuyaFresh_controlhigh   = PuyaFreshHigh - PuyaFreshControl,
  
  PuyaSulfidic_controllow    = PuyaSulfidicLow - PuyaSulfidicControl,
  PuyaSulfidic_controlmedium = PuyaSulfidicMedium - PuyaSulfidicControl,
  PuyaSulfidic_controlhigh   = PuyaSulfidicHigh - PuyaSulfidicControl,
  
  TacoFresh_controllow    = TacoFreshLow - TacoFreshControl,
  TacoFresh_controlmedium = TacoFreshMedium - TacoFreshControl,
  TacoFresh_controlhigh   = TacoFreshHigh - TacoFreshControl,
  
  TacoSulfidic_controllow    = TacoSulfidicLow - TacoSulfidicControl,
  TacoSulfidic_controlmedium = TacoSulfidicMedium - TacoSulfidicControl,
  TacoSulfidic_controlhigh   = TacoSulfidicHigh - TacoSulfidicControl,
  levels = design
)

# Use the DGEGLM object to perform likelihood ratio test based on the contrast matrix of comparisons for each drainage
# Tacotalpa lab in sulfidic control and sulfidic exposed 
TCL_lrt = glmLRT(fit, contrast = plasticity.contrasts[,"TacoSulfidic_controllow"])
TCM_lrt = glmLRT(fit, contrast = plasticity.contrasts[,"TacoSulfidic_controlmedium"])
TCH_lrt = glmLRT(fit, contrast = plasticity.contrasts[,"TacoSulfidic_controlhigh"])

# Summary of differential expression that was up, down or not significant in each comparison 
summary(decideTestsDGE(TCL_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(TCM_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(TCH_lrt, adjust.method="BH", p.value = 0.05))

# Extract up- (logFC > 0) and downregulated (logFC < 0) out separately for each drainage based on FDR cut off of 0.05
# Note, for TCL_plast, you will get an error due to no evidence of differential expression (0 up and downregulated)
TCL_plast = topTags(TCL_lrt, sum(summary(decideTestsDGE(TCL_lrt, p=0.05, adjust="BH"))[c(1,3),]))
TCM_plast = topTags(TCM_lrt, sum(summary(decideTestsDGE(TCM_lrt, p=0.05, adjust="BH"))[c(1,3),]))
TCH_plast = topTags(TCH_lrt, sum(summary(decideTestsDGE(TCH_lrt, p=0.05, adjust="BH"))[c(1,3),]))

# Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05
# Tacotalpa sulfidic
# Note, for PSOCL.up and PSOCL.down, you will get an error due to no evidence of differential expression (0 up and downregulated)
PSOCL.up   = rownames(TCL_plast[TCL_plast$table$logFC > 0,])
PSOCL.down = rownames(TCL_plast[TCL_plast$table$logFC < 0,])

PSOCM.up   = rownames(TCM_plast[TCM_plast$table$logFC > 0,])
PSOCM.down = rownames(TCM_plast[TCM_plast$table$logFC < 0,])

PSOCH.up   = rownames(TCH_plast[TCH_plast$table$logFC > 0,])
PSOCH.down = rownames(TCH_plast[TCH_plast$table$logFC < 0,])

# Combine all three exposures into a list for Tacotalpa sulfidic
# UPREGULATED
# PSOCL.up is empty 
PSOall.UP <- unique(c(PSOCM.up, PSOCH.up))
# DOWNREGULATED
# PSOCL.down is empty
PSOall.DOWN <- unique(c(PSOCM.down, PSOCH.down))

###Intersect candidate genes (OBJECTIVE 1) for Tacotalpa with plastic genes for Tacotalpa sulfidic
#UPREGULATED
Q1Q3pso.UP <- intersect(Q1alltaco.UP, PSOall.UP)
#DOWNREGULATED
Q1Q3pso.DOWN <- intersect(Q1alltaco.DOWN, PSOall.DOWN)


#####
# Use the DGEGLM object to perform likelihood ratio test based on the contrast matrix of comparisons for each drainage
#Tacotalpa lab in nonsulfidic control and nonsulfidic exposed 
TFCL_lrt = glmLRT(fit, contrast = plasticity.contrasts[,"TacoFresh_controllow"])
TFCM_lrt = glmLRT(fit, contrast = plasticity.contrasts[,"TacoFresh_controlmedium"])
TFCH_lrt = glmLRT(fit, contrast = plasticity.contrasts[,"TacoFresh_controlhigh"])

# Summary of differential expression that was up, down or not significant in each comparison 
summary(decideTestsDGE(TFCL_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(TFCM_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(TFCH_lrt, adjust.method="BH", p.value = 0.05))

# Extract up- (logFC > 0) and downregulated (logFC < 0) out separately for each drainage based on FDR cut off of 0.05
TFCL_plast = topTags(TFCL_lrt, sum(summary(decideTestsDGE(TFCL_lrt, p=0.05, adjust="BH"))[c(1,3),]))
TFCM_plast = topTags(TFCM_lrt, sum(summary(decideTestsDGE(TFCM_lrt, p=0.05, adjust="BH"))[c(1,3),]))
TFCH_plast = topTags(TFCH_lrt, sum(summary(decideTestsDGE(TFCH_lrt, p=0.05, adjust="BH"))[c(1,3),]))

# Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05
# Tacotalpa nonsulfidic
BONCL.up   = rownames(TFCL_plast[TFCL_plast$table$logFC > 0,])
BONCL.down = rownames(TFCL_plast[TFCL_plast$table$logFC < 0,])

BONCM.up   = rownames(TFCM_plast[TFCM_plast$table$logFC > 0,])
BONCM.down = rownames(TFCM_plast[TFCM_plast$table$logFC < 0,])

BONCH.up   = rownames(TFCH_plast[TFCH_plast$table$logFC > 0,])
BONCH.down = rownames(TFCH_plast[TFCH_plast$table$logFC < 0,])

# Combine all three exposures into a list for Tacotalpa nonsulfidic
# UPREGULATED
BONall.UP <- unique(c(BONCL.up, BONCM.up, BONCH.up))
# DOWNREGULATED
BONall.DOWN <- unique(c(BONCL.down, BONCM.down, BONCH.down))

### Intersect candidate genes (OBJECTIVE 1) for Tacotalpa with plastic genes for Tacotalpa nonsulfidic
# UPREGULATED
Q1Q3bon.UP <- intersect(Q1alltaco.UP, BONall.UP)
# DOWNREGULATED
Q1Q3bon.DOWN <- intersect(Q1alltaco.DOWN, BONall.DOWN)


###
# Use the DGEGLM object to perform likelihood ratio test based on the contrast matrix of comparisons for each drainage
# Puyacatengo lab in nonsulfidic control and nonsulfidic exposed
PFCL_lrt = glmLRT(fit, contrast = plasticity.contrasts[,"PuyaFresh_controllow"])
PFCM_lrt = glmLRT(fit, contrast = plasticity.contrasts[,"PuyaFresh_controlmedium"])
PFCH_lrt = glmLRT(fit, contrast = plasticity.contrasts[,"PuyaFresh_controlhigh"])

# Summary of differential expression that was up, down or not significant in each comparison 
summary(decideTestsDGE(PFCL_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(PFCM_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(PFCH_lrt, adjust.method="BH", p.value = 0.05))

# Extract up- (logFC > 0) and downregulated (logFC < 0) out separately for each drainage based on FDR cut off of 0.05
PFCL_plast = topTags(PFCL_lrt, sum(summary(decideTestsDGE(PFCL_lrt, p=0.05, adjust="BH"))[c(1,3),]))
PFCM_plast = topTags(PFCM_lrt, sum(summary(decideTestsDGE(PFCM_lrt, p=0.05, adjust="BH"))[c(1,3),]))
PFCH_plast = topTags(PFCH_lrt, sum(summary(decideTestsDGE(PFCH_lrt, p=0.05, adjust="BH"))[c(1,3),]))

# Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05
# Puyacatengo nonsulfidic
IXACL.up   = rownames(PFCL_plast[PFCL_plast$table$logFC > 0,])
IXACL.down = rownames(PFCL_plast[PFCL_plast$table$logFC < 0,])

IXACM.up   = rownames(PFCM_plast[PFCM_plast$table$logFC > 0,])
IXACM.down = rownames(PFCM_plast[PFCM_plast$table$logFC < 0,])

IXACH.up   = rownames(PFCH_plast[PFCH_plast$table$logFC > 0,])
IXACH.down = rownames(PFCH_plast[PFCH_plast$table$logFC < 0,])

# Combine all three exposures into a list for Puyacatengo nonsulfidic
# UPREGULATED
IXAall.UP <- unique(c(IXACL.up, IXACM.up, IXACH.up))
# DOWNREGULATED
IXAall.DOWN <- unique(c(IXACL.down, IXACM.down, IXACH.down))

# Combine all three exposures into a list for Puyacatengo nonsulfidic
# UPREGULATED
Q1Q3ixa.UP <- intersect(Q1allpuya.UP, IXAall.UP)
# DOWNREGULATED
Q1Q3ixa.DOWN <- intersect(Q1allpuya.DOWN, IXAall.DOWN)

###
# Use the DGEGLM object to perform likelihood ratio test based on the contrast matrix of comparisons for each drainage
# subset counts for Puyacatengo lab in sulfidic control and sulfidic exposed 
PFL_lrt = glmLRT(fit, contrast = plasticity.contrasts[,"PuyaSulfidic_controllow"])
PFM_lrt = glmLRT(fit, contrast = plasticity.contrasts[,"PuyaSulfidic_controlmedium"])
PFH_lrt = glmLRT(fit, contrast = plasticity.contrasts[,"PuyaSulfidic_controlhigh"])

# Summary of differential expression that was up, down or not significant in each comparison 
summary(decideTestsDGE(PFL_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(PFM_lrt, adjust.method="BH", p.value = 0.05))
summary(decideTestsDGE(PFH_lrt, adjust.method="BH", p.value = 0.05))

# Extract up- (logFC > 0) and downregulated (logFC < 0) out separately for each drainage based on FDR cut off of 0.05
PFL_plast = topTags(PFL_lrt, sum(summary(decideTestsDGE(PFL_lrt, p=0.05, adjust="BH"))[c(1,3),]))
PFM_plast = topTags(PFM_lrt, sum(summary(decideTestsDGE(PFM_lrt, p=0.05, adjust="BH"))[c(1,3),]))
PFH_plast = topTags(PFH_lrt, sum(summary(decideTestsDGE(PFH_lrt, p=0.05, adjust="BH"))[c(1,3),]))

# Extract up- (logFC > 0) and downregulated (logFC < 0) out separately based on FDR cut off of 0.05
# Puyacatengo sulfidic
LLUCL.up   = rownames(PFL_plast[PFL_plast$table$logFC > 0,])
LLUCL.down = rownames(PFL_plast[PFL_plast$table$logFC < 0,])

LLUCM.up   = rownames(PFM_plast[PFM_plast$table$logFC > 0,])
LLUCM.down = rownames(PFM_plast[PFM_plast$table$logFC < 0,])

LLUCH.up   = rownames(PFH_plast[PFH_plast$table$logFC > 0,])
LLUCH.down = rownames(PFH_plast[PFH_plast$table$logFC < 0,])

# Combine all three exposures into a list for Puyacatengo sulfidic
# UPREGULATED
LLUall.UP <- unique(c(LLUCL.up, LLUCM.up, LLUCH.up))
# DOWNREGULATED
LLUall.DOWN <- unique(c(LLUCL.down, LLUCM.down, LLUCH.down))

# Combine all three exposures into a list for Puyacatengo sulfidic
# UPREGULATED
Q1Q3llu.UP <- intersect(Q1allpuya.UP, LLUall.UP)
# DOWNREGULATED
Q1Q3llu.DOWN <- intersect(Q1allpuya.DOWN, LLUall.DOWN)

##############################TWO WAY VENNS##################################################################
# Create unique variable, and construct venns 
universe.UP <- unique(c(Q1Q3pso.UP, Q1Q3bon.UP))
GroupA <- universe.UP %in% Q1Q3pso.UP
GroupB <- universe.UP %in% Q1Q3bon.UP
input.df <- data.frame(PSO=GroupA, BON=GroupB)
colnames(input.df)=c("PSO","BON")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)
# create variables with genes that have evidence of plasticity (ancestral, adaptive and canalization)
Tacoancestral.up <- universe.UP[which(input.df["PSO"] == T & input.df["BON"] == T)]
Tacoadaptive.up <- universe.UP[which(input.df["PSO"] == T & input.df["BON"] == F)]
Tacocanal.up <- universe.UP[which(input.df["PSO"] == F & input.df["BON"] == T)]

# Create unique variable, and construct venns
universe.UP <- unique(c(Q1Q3ixa.UP, Q1Q3llu.UP))
GroupC <- universe.UP %in% Q1Q3ixa.UP
GroupD <- universe.UP %in% Q1Q3llu.UP
input.df <- data.frame(IXA=GroupC, LLU=GroupD)
colnames(input.df)=c("IXA","LLU")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)
# create variables with genes that have evidence of plasticity (ancestral, adaptive and canalization)
Puyaancestral.up <- universe.UP[which(input.df["IXA"] == T & input.df["LLU"] == T)]
Puyaadaptive.up <- universe.UP[which(input.df["IXA"] == F & input.df["LLU"] == T)]
Puyacanal.up <- universe.UP[which(input.df["IXA"] == T & input.df["LLU"] == F)]

# Create unique variable, and construct venns
universe.DOWN <- unique(c(Q1Q3pso.DOWN, Q1Q3bon.DOWN))
GroupA <- universe.DOWN %in% Q1Q3pso.DOWN
GroupB <- universe.DOWN %in% Q1Q3bon.DOWN
input.df <- data.frame(PSO=GroupA, BON=GroupB)
colnames(input.df)=c("PSO","BON")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)
# create variables with genes that have evidence of plasticity (ancestral, adaptive and canalization)
Tacoancestral.down <- universe.DOWN[which(input.df["PSO"] == T & input.df["BON"] == T)]
Tacoadaptive.down <- universe.DOWN[which(input.df["PSO"] == T & input.df["BON"] == F)]
Tacocanal.down <- universe.DOWN[which(input.df["PSO"] == F & input.df["BON"] == T)]

# Create unique variable, and construct venns
universe.DOWN <- unique(c(Q1Q3ixa.DOWN, Q1Q3llu.DOWN))
GroupC <- universe.DOWN %in% Q1Q3ixa.DOWN
GroupD <- universe.DOWN %in% Q1Q3llu.DOWN
input.df <- data.frame(IXA=GroupC, LLU=GroupD)
colnames(input.df)=c("IXA","LLU")
head(input.df)
a <- vennCounts(input.df)
vennDiagram(a)
# create variables with genes that have evidence of plasticity (ancestral, adaptive and canalization)
Puyaancestral.down <- universe.DOWN[which(input.df["IXA"] == T & input.df["LLU"] == T)]
Puyacanal.down <- universe.DOWN[which(input.df["IXA"] == T & input.df["LLU"] == F)]
Puyadaptive.down <- universe.DOWN[which(input.df["IXA"] == F & input.df["LLU"] == T)]

################################################################################################
#
#           4A. Extract gene annotations for plastic expression from Table S2 BLAST output 
#           
#
################################################################################################
# All output below can be found in Table S7
######################################TACOTALPA UPREGULATED#####################################
##################ANCESTRAL####################
# make gene names from Tacotalpa candidate genes a data.frame
Tacoancestral1.up<- data.frame(Tacoancestral.up)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Tacoancestral1.up) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Taco_ancestral.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Tacoancestral1.up, by = "gene.ID")

##################ADAPTIVE####################
# make gene names from Tacotalpa candidate genes a data.frame
Tacoadaptive1.up <- data.frame(Tacoadaptive.up)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Tacoadaptive1.up) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Taco_adaptive.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Tacoadaptive1.up, by = "gene.ID")

##################CANALIZATION####################
# make gene names from Tacotalpa candidate genes a data.frame
Tacocanal1.up <- data.frame(Tacocanal.up)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Tacocanal1.up) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Taco_canal.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Tacocanal1.up, by = "gene.ID")

######################################TACOTALPA DOWNREGULATED#####################################
##################ANCESTRAL####################
#NO GENES HAD EVIDENCE OF DOWNREGULATION IN TACOTALPA

##################ADAPTIVE####################
# make gene names from Tacotalpa candidate genes a data.frame
Tacoadaptive1.down <- data.frame(Tacoadaptive.down)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Tacoadaptive1.down) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Taco_adaptive.down <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Tacoadaptive1.down, by = "gene.ID")

##################CANALIZATION####################
# make gene names from Tacotalpa candidate genes a data.frame
Tacocanal1.down <- data.frame(Tacocanal.down)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Tacocanal1.down) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Taco_canal.down <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Tacocanal1.down, by = "gene.ID")

######################################PUYACATENGO UPREGULATED#####################################
##################ANCESTRAL####################
# make gene names from Tacotalpa candidate genes a data.frame
Puyaancestral1.up <- data.frame(Puyaancestral.up)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Puyaancestral1.up) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Puya_ancestral.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Puyaancestral1.up, by = "gene.ID")

##################ADAPTIVE####################
# make gene names from Tacotalpa candidate genes a data.frame
Puyaadaptive1.up <- data.frame(Puyaadaptive.up)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Puyaadaptive1.up) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Puya_adaptive.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Puyaadaptive1.up, by = "gene.ID")

##################CANALIZATION####################
# make gene names from Tacotalpa candidate genes a data.frame
Puyacanal1.up <- data.frame(Puyacanal.up)

# add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Puyacanal1.up) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Puya_canal.up <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Puyacanal1.up, by = "gene.ID")

######################################PUYACATENGO DOWNREGULATED#####################################
##################ANCESTRAL####################
# Make gene names from Tacotalpa candidate genes a data.frame
Puyaancestral1.down <- data.frame(Puyaancestral.down)

# Add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Puyaancestral1.down) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Puya_ancestral.down <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Puyaancestral1.down, by = "gene.ID")

##################ADAPTIVE####################
# Make gene names from Tacotalpa candidate genes a data.frame
Puyaadaptive1.down <- data.frame(Puyadaptive.down)

# Add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Puyaadaptive1.down) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Puya_adaptive.down <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Puyaadaptive1.down, by = "gene.ID")

##################CANALIZATION####################
# Make gene names from Tacotalpa candidate genes a data.frame
Puyacanal1.down <- data.frame(Puyacanal.down)

# Add column heading that matches the blast output (Table S2 - Poecilia mexicana Annotations.csv)
col_headings <- c("gene.ID")
names(Puyacanal1.down) <- col_headings

## Extract swissprot accessions (subject sequence IDs) and Protein annotations from Blast output based on transcript IDs (query sequence IDs)
Puya_canal.down <- merge(BLAST_annot[, c("gene.ID", "Subject.sequence.ID", "Protein.annotations", "gene.name", "E.value")], Puyacanal1.down, by = "gene.ID")

###############################################################################################
#
#           5. Identify genes not assigned to a category
#           
#
################################################################################################
# Numbers reported in Table 2
# Combined genes with evidence of constitutive expression and/or plasticity for Tacotalpa. Note that numbers reported in table 2 are the difference between Tacotalpa candidates - Tacoall.UP/Tacoall.DOWN
# UPREGULATED
Tacoall.UP <- unique(c(Tacoancestral.up, Tacoadaptive.up, Tacocanal.up, Evolvedtaco.up))
# DOWNREGULATED
Tacoall.DOWN <- unique(c(Tacoancestral.down, Tacoadaptive.down, Tacocanal.down, Evolvedtaco.down))

# Combined genes with evidence of constitutive expression and/or plasticity for Puyacatengo. Note that numbers reported in table 2 are the difference between Puyacatengo candidates - Puyaall.UP/Puyaall.DOWN
# UPREGULATED
Puyaall.UP <- unique(c(Puyaancestral.up, Puyaadaptive.up, Puyacanal.up, Evolvedpuya.up))
# DOWNREGULATED
Puyaall.DOWN <- unique(c(Puyaancestral.down, Puyadaptive.down, Puyacanal.down, Evolvedpuya.down))
