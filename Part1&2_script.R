#install packages
install.packages("tidyr")
install.packages("dplR")
install.packages("pheatmap")
install.packages("reshape")
install.packages("gridExtra")
install.packages("ggplot2")
library(ggplot2)
library(tidyr)
library(dplR)
library(pheatmap)
library(reshape)
library(grid)
library(gridExtra)
library(png)


#Part 1 – Gene Expression Analysis
# ==============================================================================
# Task 1a : Heatmap plot of HBR and UHR
# ==============================================================================
#Loading the dataset
file_onlinez <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv"
Gene_counts <- read.csv(file_onlinez, header = T )


colnames(Gene_counts[,2:7])
plot.1a <- pheatmap(mat = Gene_counts[,2:7],
                    border_color = "Black",
                    legend = T,
                    fontsize_row = 7,
                    fontsize_col = 7,
                    labels_row = Gene_counts[,1],  #check this line
                    cluster_rows = T,
                    cluster_cols = T,
                    color = blues9)
# ==============================================================================
# Task 1b: volcano Plot log2FoldChange vs log10(Padj) from the DEG results
# ==============================================================================
#Load dataset
file1 <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv"
DEG <- read.csv(file1, header = T)

#plot
plot(x = DEG$log2FoldChange,
     y = DEG$X.log10PAdj,
     col = as.factor(DEG$significance),
     main = "volcano plot: log2FoldChange vs log10(Padj)",
     xlab = "log2Foldchange",
     ylab = "log10PAdj",
     pch = 19,
     cex = 0.7,
     abline(v = c(-2,2), h = c(1, 0), lty = 2,)
     
)
legend("topright",
       legend = c("Downregulated", "Not Significant", "Upregulated"),
       col = c(1, 2, 3),
       pch = 19)


# Part 2 – Breast Cancer Data Exploration
# ==============================================================================
# Task2c : scatter plot radius vs texture
# ==============================================================================
#lOADING THE DATASET
file_online2 <- "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
breast_cancer <- read.csv(file_online2, header = T)

ggplot(breast_cancer,
       aes(x = radius_mean, y = texture_mean, color = as.factor(diagnosis)))+
  geom_point(size = 2, alpha = 1)+
  coord_cartesian(xlim = c(0,25), ylim = c(0,40))+
  labs( x = "radius_mean",
        y = "texture_mean",
        color = "diagnosis",
        title = "Texture_mean Vs Radius_mean")+
  theme_classic()

# ==============================================================================
# Task2d : correlation heatmap 
# ==============================================================================
#Reshape the data
heatmap_cor_df <- breast_cancer[, c("radius_mean", "texture_mean", "perimeter_mean", 
                                    "area_mean", "smoothness_mean", "compactness_mean")]

#compute correlation matrix
cor_matrix <- cor(heatmap_cor_df)
cor_df <- melt(cor_matrix) 

#Rename column names manually

colnames(cor_df) <- c(variable.name = c("variable1","variable2"),
                      value.name = "value")
#plot
ggplot(cor_df,
       aes(x = variable2, y = variable1, fill = value ))+
  geom_tile(col = "Black")+
  geom_text(aes(label = round(value, 2)), size = 3)+
  scale_fill_gradient(low = "white", high = "blue")+
  labs(x = "",
       y = "",
       title = "Correlation Heatmap")+
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8,hjust = 1 ,angle = 90))


#Task 2e: scatter plot smoothness_mean vs compactness_mean
ggplot(breast_cancer,
                  aes(x = smoothness_mean, y = compactness_mean, colour = as.factor(diagnosis)))+
  geom_point(size = 2, alpha = 0.8)+
  labs(x = "smoothness_mean",
       y = "compactness_mean",
       title = "Smoothness_mean vs Compactness_mean",
       color = "diagnosis")+
  theme_classic()

# ==============================================================================
# Task2f : density plot area of distribution
# ==============================================================================
ggplot(breast_cancer,
                  aes(x = area_mean, fill = as.factor(diagnosis)))+
  geom_density(alpha = 0.5)+
  labs(x = "area_mean",
       y = "Density",
       fill = "Diagnosis",
       title = "Area Distribution")+
  theme_classic()
