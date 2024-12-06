#### Import data ####
## Import raw counts
library(readr)
raw_counts <- read_delim("", #  <- Add path to the featureCounts output file named 'raw_counts_table.txt"
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
class(raw_counts) <- "data.frame"
rownames(raw_counts) <- raw_counts$Geneid
raw_counts <- raw_counts[,-c(1:6)]
colnames(raw_counts) <- gsub("/home/OSUMC.EDU/lafe05/network/X/Labs/Sehgal/Alessandro/Nazia_CRC/bam_sorted/","", colnames(raw_counts))
colnames(raw_counts) <- gsub("_sorted.bam","", colnames(raw_counts))

## Make sample info table
sample_info <- data.frame(sample_id = colnames(raw_counts),
                          condition = substring(colnames(raw_counts),1,c(nchar(colnames(raw_counts))-2)))

## Import gene annotation
library(readxl)
gencode_v43_annotation <- read_excel("") # <- add path to excel file named 'gencode_v43_annotation.xlsx' containg the gene annotation (this file is also available on GitHub)

#### Filtering ####
## Calculate RPM
rpm_matrix <- raw_counts
total_mapped_reads <- apply(rpm_matrix, 2, sum)
rpm_matrix <- rpm_matrix*(10^6)

for (c in 1:ncol(rpm_matrix)){
  rpm_matrix[,c] <- rpm_matrix[,c]/total_mapped_reads[c]
}

## Filter low-expressed genes
cutoff <- 1
drop <- which(apply(rpm_matrix,1, mean) < cutoff)
d <- raw_counts[-drop,]

#### Normalization and Differential Expression Analysis ####
## Make design matrix
library(limma)
library(edgeR)

design_matrix <- model.matrix(~0+sample_info$condition)
sample_types <- levels(factor(sample_info$condition))
colnames(design_matrix) <- sample_types

## Normalization and differential expression analysis
voomTransformed <- voom(d, design_matrix, plot=TRUE)
voomed.fitted <- lmFit(voomTransformed, design=design_matrix)

contrasts <- makeContrasts(T-N,
                           UT-N,
                           levels = design_matrix)

fit_2 <- contrasts.fit(voomed.fitted, contrasts)
fit_2 <-eBayes(fit_2)

## PCA
library(ggfortify)
pca_res <- prcomp(t(voomTransformed$E), scale. = TRUE)
autoplot(pca_res, data = sample_info, colour = 'condition', frame = FALSE, label = FALSE, size = 3) +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"),
        legend.key.size = unit(1, 'cm'),legend.key.height = unit(1, 'cm'), legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=18), legend.text = element_text(size=14))

#### Extract results ####
library("writexl")

## Make average norm expr table
norm_expr <- rpm_matrix[-drop,]
norm_expr$ENSEMBL <- rownames(norm_expr)
norm_expr <- norm_expr[,c(ncol(norm_expr),1:c(ncol(norm_expr)-1))]

norm_expr$Avg_RPM_N <- as.vector(apply(norm_expr[,c(2:3)],1,mean))
norm_expr$SD_N <- as.vector(apply(norm_expr[,c(2:3)],1,sd))

norm_expr$Avg_RPM_T <- as.vector(apply(norm_expr[,c(4:5)],1,mean))
norm_expr$SD_T <- as.vector(apply(norm_expr[,c(4:5)],1,sd))

norm_expr$Avg_RPM_UT <- as.vector(apply(norm_expr[,c(6:7)],1,mean))
norm_expr$SD_UT <- as.vector(apply(norm_expr[,c(6:7)],1,sd))

norm_expr <- norm_expr[,c(1,8:13)]

norm_expr <- merge(x = norm_expr, y = gencode_v43_annotation, by = "ENSEMBL")

## T-N
T_vs_N_res <- topTable(fit_2, coef = 1, number = Inf, adjust.method = "BH")
T_vs_N_res$ENSEMBL <- rownames(T_vs_N_res)
T_vs_N_res <- merge(x = T_vs_N_res, y = norm_expr, by = "ENSEMBL")
T_vs_N_res$expr_diff <- T_vs_N_res$Avg_RPM_T - T_vs_N_res$Avg_RPM_N

i <- 1

T_vs_N_res$test <- NA

while(i <= nrow(T_vs_N_res)){
  if(sign(T_vs_N_res[i,2])==sign(T_vs_N_res[i,20])){T_vs_N_res[i,21]<-"keep"}else{T_vs_N_res[i,21]<-"remove"}
  i <- i + 1
}

T_vs_N_res <- T_vs_N_res[T_vs_N_res$test=="keep",,drop=FALSE]
T_vs_N_res <- T_vs_N_res[,c(1,14:19,8:11,2,5,6)]
write_xlsx(T_vs_N_res,
           "ADD_OUPUT_FOLDER/T_vs_N_res.xlsx") # <- add path for output excel file containg the results

## UT-N
UT_vs_N_res <- topTable(fit_2, coef = 2, number = Inf, adjust.method = "BH")
UT_vs_N_res$ENSEMBL <- rownames(UT_vs_N_res)
UT_vs_N_res <- merge(x = UT_vs_N_res, y = norm_expr, by = "ENSEMBL")
UT_vs_N_res$expr_diff <- UT_vs_N_res$Avg_RPM_T - UT_vs_N_res$Avg_RPM_N

i <- 1

UT_vs_N_res$test <- NA

while(i <= nrow(UT_vs_N_res)){
  if(sign(UT_vs_N_res[i,2])==sign(UT_vs_N_res[i,20])){UT_vs_N_res[i,21]<-"keep"}else{UT_vs_N_res[i,21]<-"remove"}
  i <- i + 1
}

UT_vs_N_res <- UT_vs_N_res[UT_vs_N_res$test=="keep",,drop=FALSE]
UT_vs_N_res <- UT_vs_N_res[,c(1,14:19,8,9,12,13,2,5,6)]
write_xlsx(UT_vs_N_res,
           "ADD_OUPUT_FOLDER/UT_vs_N_res.xlsx") # <- add path for output excel file containg the results


