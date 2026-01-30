counts <- read.csv("hil-1_recovery_Starved_counts.csv",header = T)
WS273_geneNames <- read.csv("WS273_geneNames.csv",header = T)
counts_merge <- merge (counts, WS273_geneNames, by.x = "WB_ID", by.y = "WB_id") 
counts_byBiotype <- counts_merge[,-1]
counts_byBiotype <- counts_byBiotype[,-13:-15]
library(edgeR)
library (ggplot2)
library(reshape2)
counts_byBiotype_melt <- melt(counts_byBiotype)

byBiotype<- aggregate (counts_byBiotype_melt$value,list(variable=counts_byBiotype_melt$variable,type=counts_byBiotype_melt$type),sum)

proteinCoding<-subset(counts_merge,type=="protein_coding_gene") 
counts2<-proteinCoding[,1:13]
rownames(counts2)<-counts2$WB_ID
counts2 <- counts2[,-1]

counts_groups <- c("N2Starved", "N2Starved", "N2Starved", "N2Starved","hil1Starved", "hil1Starved", "hil1Starved", "hil1Starved","hil1Starved", "hil1Starved", "hil1Starved", "hil1Starved")

d2<-DGEList(counts=counts2,group=factor(counts_groups))
dim(d2) #check the dimensions of the d2 without any filtering
keep_filter<-rowSums(cpm(d2)>1)>=4 #decide how you want to filter the data
d2<-d2[keep_filter,] #filter d2 to only include genes that passed filtering
dim(d2) #check the new dimensions

cpm_d2<-cpm(d2,normalized.lib.sizes = TRUE) #make a counts per million object containing normalized CPM
cpm_d2<-as.data.frame(cpm_d2)
cpm_d2_melt<-melt(cpm_d2)

conditions<-counts_groups
cpm_d2_df<-data.frame(cpm_d2)
cpm_d2_df$mean<-rowMeans(cpm_d2_df)
cpm_d2_df2<-cpm_d2_df[,1:12]/cpm_d2_df$mean #mean normalize
cpm_d2_df2<-log2(cpm_d2_df2+1) #log2 transform 
pca = prcomp(t(cpm_d2_df2)) #principal component analysis (PCA) on the log2 mean normalized CPM values
summary(pca)

write.table(cpm_d2,"CPM Starved.txt",quote=F,sep="\t")

pca_genes<-pca$x
pca_genes_dataframe<-as.data.frame(pca_genes)
pca_genes_dataframe<-data.frame(conditions,pca_genes_dataframe)

replicates<-c("rep2","rep3","rep4","rep5","rep2","rep3","rep4","rep5","rep6","rep7","rep8", "rep9")
PCA1<-(ggplot(pca_genes_dataframe,aes(x=PC1,y=PC2,colour=conditions))+
         geom_point(size=5)+
         ggtitle("PCA of hil-1 starved, 95% CI")+
         labs(x="PC1 (20.62% of variance)",y="PC2 (17.8% of variance)")+
         stat_ellipse(level = 0.95)+theme_classic(base_size = 15)+
         theme(aspect.ratio = 1))
PCA1


PCA2<-(ggplot(pca_genes_dataframe,aes(x=PC1,y=PC2,colour=conditions,shape=replicates))+
         geom_point(size=5)+
         ggtitle("PCA of hil-1 starved, replicates")+
         labs(x="PC1 (20.62% of variance)",y="PC2 (17.8% of variance)")+
         theme_classic(base_size = 15)+
         theme(legend.key.size = unit(.2, "cm")) +
         theme(aspect.ratio = 1)) 
PCA2

count_cormatrix<-round(cor(log2(cpm_d2+1),use="all.obs",method="pearson"),digits=2)

reorder_cormat <- function(count_cormatrix){
  # Use correlation between variables as distance
  dd <- as.dist((1-count_cormatrix)/2)
  hc <- hclust(dd)
  cormat <-count_cormatrix[hc$order, hc$order]
}

#use the count_cormatrix that is not reordered. 
count_cormatrix1<-reorder_cormat(count_cormatrix)
melted_cormatrix1<-melt(count_cormatrix1)
melted_cormatrix<-melt(count_cormatrix)
correlation_matrix<-(ggplot(data=melted_cormatrix,aes(x=Var1,y=Var2,fill=value))+
                       geom_tile(colour="white")+
                       scale_fill_gradient2(low="#900C3F",high="#6B33FF",mid="white",midpoint=0.95,limit=c(0.9,1))+
                       theme_classic(base_size = 10)+
                       theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
                       geom_text(aes(Var2, Var1, label = value), color = "black", size = 3)+
                       ggtitle("Corr. matrix of hil-1 starved")+theme(aspect.ratio = 1))
correlation_matrix

d2<-calcNormFactors(d2) #calculate normalization factors
d2<-estimateCommonDisp(d2) #calculate dispersion 
d2<-estimateTagwiseDisp(d2)

#hil-1/N2
de.tag<-exactTest(d2,d2$tagwise.dispersion,pair=c("N2Starved","hil1Starved")) #make sure to have these in the right order. The control, or denominator, should come first.
de.tag_sort<-topTags(de.tag,n=nrow(de.tag$table))$table
de.tag_top<-rownames(de.tag_sort)[de.tag_sort$FDR<=0.1]
de.tag_merge1<-merge(WS273_geneNames,de.tag_sort,by.x="WB_id",by.y = 0)
plotSmear(de.tag,de.tags = de.tag_top,main="hil-1starved/N2starved", cex = 1.5,cex.axis = 2, cex.lab = 2, cex.main = 2)
write.table(de.tag_merge1,"hil1pooledvsN2Starved.txt",quote=F,sep="\t")


#CelEsT plot
celest <- read.table(
  "celest_result_1.28.26.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
celest$Padj <- p.adjust(celest$Pval, method = "BH")
celest <- celest %>%
  mutate(log10P = -log10(Pval))

sig_thresh <- 0.05

label_df <- celest %>%
  filter(Pval < 0.01 & abs(Score) > 2)
sig_thresh <- 0.05

ggplot(celest, aes(x = Score, y = log10P)) +
  geom_point(
    aes(
      size  = log10P,
      color = Padj < sig_thresh
    ),
    alpha = 0.8
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = Gene_name),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  scale_size(range = c(1, 5)) +  
  scale_color_manual(
    values = c("grey70", "red"),
    guide = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    x = "TF activity",
    y = expression(-log[10](p-value))
  ) +
  theme_classic() + coord_fixed()

write.csv(celest, file = "celest_padj.csv")

#Comparisons to other datasets
nipicebp<-read.csv("hil1_nipi3_cebp.csv", header = T)
nipicebp_lists <- nipicebp[, c(1:8, 16)] 

hil1down <- as.list(nipicebp_lists$down_hil1)
hil1down <- hil1down[hil1down != ""]

hil1up <- as.list(nipicebp_lists$up_hil1)
hil1up <- hil1up[hil1up != ""]

nipi3down <- as.list(nipicebp_lists$down_nipi3)
nipi3down <- nipi3down[nipi3down != ""]

nipi3up <- as.list(nipicebp_lists$up_nipi3)
nipi3up <- nipi3up[nipi3up != ""]

cebp1 <- as.list(nipicebp_lists$cebp1_targets)
cebp1 <- cebp1[cebp1 != ""]

pmk1down <- as.list(nipicebp_lists$pmk1_dep_up)
pmk1down <- pmk1down[pmk1down != ""]

pmk1up <- as.list(nipicebp_lists$pmk1_dep_down)
pmk1up <- pmk1up[pmk1up != ""]

downpa14 <- as.list(nipicebp_lists$down_fletcher)
downpa14 <- downpa14[downpa14 != ""]

uppa14 <- as.list(nipicebp_lists$up_fletcher)
uppa14 <- uppa14[uppa14 != ""]


#Figure 4c: Comparisons to PA14 (Fletcher 2019)

fit <- euler(c(
  "Up in hil-1" = 177,
  "Down in hil-1" = 192,
  "Down in hil-1&Up on PA14" = 40,
  "Down on PA14" = 772,
  "Up on PA14" = 825,
  "Up in hil-1&Up on PA14" = 24,
  "Up in hil-1&Down on PA14" = 20,
  "Down in hil-1&Down on PA14" = 11
))

plot(fit,
     fills = list(fill = c("skyblue", "salmon", "lightgreen"), alpha = 0.6),
     labels = list(cex = 1.2),
     quantities = TRUE,
     main = "Venn Diagram"
)



#Stats
phyper(10, 803, 12265, 242,lower.tail = FALSE,log.p = FALSE) #down in both, p=0.89
phyper(39, 890, 12178, 242,lower.tail = FALSE,log.p = FALSE) #up in PA14 down in hil-1, p=1.39e-07
phyper(19, 803, 12265, 221,lower.tail = FALSE,log.p = FALSE) #up hil-1, down PA14
phyper(23, 890, 12178, 221,lower.tail = FALSE,log.p = FALSE) #up in both
#Bonferroni corrections applied

