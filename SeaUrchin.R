data <- read.table(file = "Pliv_aH2p.stages.counts.tsv", header = TRUE)
colnames(data)[7:28] <- c("X128", "X16", "X2", "X32", "X4", "X64", "X8", "A", "B", "C", 
                        "D", "EB", "eggs", "EG", "E", "ID16", "LG", "MB1", "MB2", 
                        "PHB", "P", "SB")


######## Calcul of FPKM and TPM #####


countToFpkm <- function(count){
  n <- as.numeric(sum(count))
  (count * 10**6 *10**3) / (n * data$Length)
}


descrp <- data[, 1:6]

countToFpkmAll <- descrp
for(i in data[, 7:28]){
  countToFpkmAll <- cbind(countToFpkmAll, countToFpkm(i))
}
colnames(countToFpkmAll)[7:28] <- colnames(data)[7:28]


countToTpm <- function(count) {
  x <- count / data$Length
  x*10**6 / sum(x)
}


countToTpmAll <- descrp
for(i in data[, 7:28]){
  countToTpmAll <- cbind(countToTpmAll, countToTpm(i))
}

colnames(countToTpmAll)[7:28] <- colnames(data)[7:28]


###### Heatmaps ####


library(pheatmap)


correct_name <- function(x){
  ord=c("eggs","X2","X4","X8","X16","X32","X64","X128","EB","PHB", "SB","MB1","MB2","EG", "LG","ID16","P","A","B","C","D","E")
  names=(c("eggs","2_cells","4_cells","8_cells","16_cells","32_cells","64_cells","128_cells","Early_blastula","Blastula",
           "Swimming_blastula","Mesenchyme_blastula_1","Mesenchyme_blastula_2","Early_gastrula", "Late_gastrula","Prism","Pluteus",
           "Ovaries","Gut","Bodywall","Tubefeet","Lantern"))
  x2 <- x[ord,ord]
  colnames(x2) <- names
  rownames(x2) <- names
  return(x2)
}


correct_name_col <- function(x){
  ord=c("eggs","X2","X4","X8","X16","X32","X64","X128","EB","PHB", "SB","MB1","MB2","EG", "LG","ID16","P","A","B","C","D","E")
  names=(c("eggs","2_cells","4_cells","8_cells","16_cells","32_cells","64_cells","128_cells","Early_blastula","Blastula",
           "Swimming_blastula","Mesenchyme_blastula_1","Mesenchyme_blastula_2","Early_gastrula", "Late_gastrula","Prism","Pluteus",
           "Ovaries","Gut","Bodywall","Tubefeet","Lantern"))
  x2 <- x[,ord]
  colnames(x2) <- names
  return(x2)
}


############ General ############


all_Fpkm_norma <- log10((countToFpkmAll[,7:28] + 1))
all_Fpkm_norma_eucli <-  as.matrix(dist(t(all_Fpkm_norma)),method = "euclidean")
all_Fpkm_norma_eucli_c <- correct_name(all_Fpkm_norma_eucli)
all_Fpkm_norma_eucli_c[all_Fpkm_norma_eucli_c == 0] <- NA
pheatmap(all_Fpkm_norma_eucli_c)

all_Tpm_norma <- log10((countToTpmAll[,7:28] + 1))
all_Tpm_norma_eucli <-  as.matrix(dist(t(all_Tpm_norma)),method = "euclidean")
all_Tpm_norma_eucli_c <- correct_name(all_Tpm_norma_eucli)
all_Tpm_norma_eucli_c[all_Tpm_norma_eucli_c == 0] <- NA
pheatmap(all_Tpm_norma_eucli_c)


######### Interresting genes ##########


vecgene <- c("PL34811", "PL09343", "PL08133", "PL10416", "PL27413", "PL17193", "PL09726", "PL09256", "PL09288", "PL23455",
             "PL10256", "PL22067", "PL01842", "PL05901", "PL39104", "PL14032", "PL05498", "PL37403", "PL23999", "PL06744", 
             "PL40150", "PL31322", "PL07173", "PL36093", "PL20521")


gene_name <- function(x){
  df <- data.frame()
  for(i in vecgene){
    df <- rbind(df, x[x$Geneid == i, ])
  }
  rownames(df) <- c("Nodalin", "Lefty", "bmp2/4", "goosecoid", "chordin", "Not", "NK2.2", "foxA", "foxG", "fgfr1",
                    "univin", "ADMP1", "wnt2", "gata1/2", "ese", "tbx2/3", "smad6", "msx", "glypican5", "Dlx", "irxA", 
                    "atbf1", "wnt5", "ADMP2", "gcm")
  mat <- as.data.frame(df[,])
  return(mat)
}


gene <- countToTpmAll[countToTpmAll$Geneid %in% vecgene,]
abondance_gene <- gene_name(gene)
abondance_gene_norma <- log10((abondance_gene[,7:28] + 1))
abondance_gene_norma_c <- correct_name_col(abondance_gene_norma)
pheatmap(abondance_gene_norma_c)

gene_norma <- log10((gene[,7:28] + 1))
gene_norma_dist <-  as.matrix(dist(t(gene_norma),method = "euclidean"))
gene_norma_dist_c <- correct_name(gene_norma_dist)
gene_norma_dist_c[gene_norma_dist_c == 0] <- NA
pheatmap(gene_norma_dist_c)

#################### PCA ##################


library(viridis)


Tpm_all_c <- correct_name_col(countToTpmAll[, 7:28])
Fpkm_all_c <- correct_name_col(countToFpkmAll[, 7:28])


pca_res <- prcomp(t(as.matrix(log10(Tpm_all_c[,1:17]+1))), scale = FALSE)
summary(pca_res)
pca_res$x[,1]
par(xpd = T, mar = par()$mar + c(0,0,0,11))
plot(pca_res$x[,1], pca_res$x[,3], col=viridis_pal(option = "D")(17), pch = 19, xlab = "PC1", ylab = "PC3")
legend(66, 30,legend = c("8_cells","16_cells","2_cells","32_cells","4_cells","64_cells","128_cells","eggs","Early_blastula",
                           "Blastula","Swimming_blastula","Mesenchyme_blastula_2","Mesenchyme_blastula_1","Early_gastrula",
                           "Late_gastrula" ,"Pluteus", "Prism"), col=viridis_pal(option = "D")(17), lty=1:2)
par(mar=c(5, 4, 4, 2) + 0.1)


######## Clustering ###########


library(Mfuzz)
help(standardise)

Tpm_stade <- Tpm_all_c[1:17]
Tpm_stade_log <- as.matrix(log10(Tpm_stade +1 ))
Tpm_stade_red <- as.matrix(Tpm_stade[,c(1,6,9,10,11,13,15,16,17)])
table(complete.cases(Tpm_stade_red))
Tpm_stade_log_red=as.matrix(Tpm_stade_log[,c(1,6,9,10,11,13,15,16,17)])


Tpm_set <- new("ExpressionSet",exprs=Tpm_stade_log_red)
Tpm_set <- filter.std(Tpm_set,min.std=0, visu = TRUE)
Tpm_set <- standardise(Tpm_set)
Tpm_m <- mestimate(Tpm_set)
stg.Rn <- c('Egg','32C','EBl','B','SBl','MBl','Ga','Pr','Pl')


Tpm_dmin = Dmin(Tpm_set,Tpm_m,crange=seq(4,48,4),repeats=3,visu=TRUE)
Tpm_clust <-mfuzz(Tpm_set,c=28,m=Tpm_m)
pdf("mfuzz.plot.pdf")
mfuzz.plot(Tpm_set,cl=Tpm_clust,mfrow=c(3,3),time.labels=stg.Rn, new.window = F)
dev.off()

Tpm_over <- overlap(Tpm_clust)
overlap.plot(Tpm_clust ,over=Tpm_over,thres=0.05)

