#install packages
install.packages("tidyr")
install.packages("dplR")
install.packages("pheatmap")
install.packages("reshape")
install.packages("gridExtra")
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

ggsave("panel_volcano.png", )

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
plot.2e <- ggplot(breast_cancer,
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
plot.2f <- ggplot(breast_cancer,
       aes(x = area_mean, fill = as.factor(diagnosis)))+
  geom_density(alpha = 0.5)+
  labs(x = "area_mean",
       y = "Density",
       fill = "Diagnosis",
       title = "Area Distribution")+
  theme_classic()


#Part 3
#Install the following packages

install.packages("readxl")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("igraph")
install.packages("dplyr")
install.packages("tidyr")
library(readxl)
library(ggplot2)
library(dplyr)
library(igraph)
library(tidyr)
library(stringr)

#Helper function for transparent colors
transparent_color <- function(color, percent = 50, name = NULL){
  # color = color name
  # percent = % transparency
  # name = name of color
  
  ## Get rbg values for color
  
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and apha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 /100,
               names = name)
  
  ## save the color
  invisible(t.col)
}

#custom color palette
hb_pal <- c("#4e79a7", 
            "#8cd17d", 
            "#e15759", 
            "#fabfd2", 
            "#a0cbe8", 
            "#59a14f", 
            "#b07aa1", 
            "#ff9d9a", 
            "#f28e2b", 
            "#f1ce63",
            "#79706e",
            "#d4a6c8",
            "#e9e9e9",
            "#ffbe7d",
            "#bab0ac",
            "#9d7660",
            "#d37295",
            "#86bcb6",
            "#362a39",
            "#cd9942")


# test the color pallete
plot(1:length(hb_pal), 1:length(hb_pal),
     col = hb_pal, pch = 19)

# ==============================================================================
# Task 1 : Reproduce panel 2a: Cell-type ratio distributions
# ==============================================================================
## Read sheet a 
a.sheet <- read_excel("~/Downloads/hb_stage_2.xlsx", sheet = "a")

## Produce a Box plot New_ratio grouped by Cell_type

plot.3a <- ggplot(a.sheet,
       aes(x = cell_type , y = new_ratio ))+
  geom_boxplot(fill = hb_pal[c(1:10)], alpha = 1)+
  labs(x = "Cell_type",y = "Ratio",title = "Cell-Type Ratio Distribution",
  )+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





# ==============================================================================
#Task 2 : Reproduce panel 2b: Half-life vs alpha-life scatter
# ==============================================================================
##Read sheet b 
b.sheet <- read_excel("~/Downloads/hb_stage_2.xlsx", sheet = "b")

##create a subset based on threshold
b.sheet <- b.sheet %>% mutate(
  regime = case_when(
    log2(alpha) <= -3.5 & log2(half_life) >= 2.5 ~ "Regime1", #low alpha, high_life
    log2(alpha) <= -3.5 & log2(half_life) < 2.5 ~ "Regime2", #low alpha, low high_life
    log2(alpha) > -3.5 & log2(half_life) >= 2.5 ~ "Regime3", #high alpha, high high_life
    log2(alpha) > - 3.5 & log2(half_life) < 2.5 ~ "Regime4") #high alpha, low high_life
)

#create color
my_palette <- c("Regime1" = "#a0cbe8","Regime2" = "#79706e","Regime4" = "#8cd17d",
                "Regime3" =  "#e15759")

#Plot
plot.3b <- ggplot(b.sheet, aes(log2(half_life), log2(alpha))) +
  geom_point(aes(colour = regime), size = 2) +
  scale_color_manual(values = my_palette, name = "Regime") +
  geom_hline(yintercept = -3.5, linetype = "dashed") +
  geom_vline(xintercept = 2.5, linetype = "dashed") +
  
  #label and title
  labs(
    x = "log2(half_life)",
    y = "log2(alpha)",
    title = "Half-life vs alpha-life scatter"
  )+
  geom_text(data = subset(b.sheet, cell %in% c("Ccr2", "Camp")),
            aes(label = cell))+
  theme_classic()+
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold.italic"),
    axis.title = element_text(face = "bold.italic"))


#Why log2?
#Log2 helps normalize skewed distribution making the data normally distributed,
#Each unit chnage in log2 represent a doubline or halving in the original scale.

#What does the four quadrants mean?
#Top left quadrants represent low alpha and long half life,meaning these genes have low degradation rate.
#Bottom left quadrant represent low alpha and short half life,its biologically unsual low alpha should correlate with high half life.
#Toop right quadrant represent high alpha and a long half life, another biological unusual combination.Might represent genes with complex degradation.
#Bottom right represent unstable transcrip degradation ( high alpha and short half life). likely housekeeping genes.


# ==============================================================================
# Task 3 : Reproduce Heatmap across cell types and time
# ==============================================================================
#Read sheet c 
c.sheet <- read_excel("~/Downloads/hb_stage_2.xlsx", sheet = "c")

#convert to matrix
c.sheet <- read_excel("~/Downloads/hb_stage_2.xlsx", sheet = "c")
c.sheet_matrx <- as.matrix(c.sheet)   

#colnames and rownames
mat_c <- c.sheet_matrx[,2:50]
mat_c <- apply(mat_c,2, as.numeric)

col_names <- colnames(mat_c)
row_names <- c.sheet[1]

#extract cell type and time from a character vector
cell_type <- sub("[0-9]+h$", "", col_names)
time <- str_extract(col_names, "[0-9]+h")


#create column annotation dataframe
annotation_col <- data.frame(cellType = cell_type, Time = time, 
                             row.names = col_names)
#Define heatmap color
heat_col <- colorRampPalette(c("#d73027","white","#4e79a7"))(100)
max_abs <- max(abs(mat_c), na.rm = T)
breaks <- seq(-max_abs, max_abs, length.out = 130)

#plot
pheatmap(mat = mat_c,
         border_color = "grey",
         color = heat_col, breaks = breaks,
         annotation_col = annotation_col,
         legend = T,
         labels_row = c.sheet_matrx[,1],
         cluster_rows = T,
         cluster_cols = F,
         main = "Heatmap across Cell type and Time",
         show_rownames = F, show_colnames = F)

##Conceptual check
#Why clustering gene but not time 
## Clustering gene reveal co-regulated modules.
## Certain genes share the same biological functions
## Time follows a natural order 


# ==============================================================================
# Task 4 : Pathway enrichment map
# ==============================================================================
#read sheet d_1
d.sheet <- read_excel("~/Downloads/hb_stage_2.xlsx", sheet = "d_1")

#restructure the data
d.sheet_df1 <- d.sheet%>%
  select(pathway,2:8)%>%
  pivot_longer(
    cols = -pathway,
    names_to = "variable",
    values_to = "values"
  )

#control row order
d.sheet_df1$pathway <- factor(d.sheet_df1$pathway,
                              levels = unique(d.sheet_df1$pathway))

#plot
ggplot(d.sheet_df1,
       aes(x = variable, y = pathway, fill = values))+
  geom_tile(color = "grey60")+
  scale_fill_gradient2(low = "#e15759" , mid = "white",high = "royalblue")+
  scale_y_discrete(limit = rev)+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, 
                               colour = "Black"),
    axis.text.y = element_text(size = 9, face = 'plain', colour = "Black", 
                               hjust = 1),
    axis.title = element_blank())

##No clustering here because clustering would reorder pathway, breking biological logic
##Diverging palette clearly seperate positive vs negative responses.

# ==============================================================================
# Task 5 : Bubble plot of Kinetic regimes
# ==============================================================================
sheet.e <- read_excel("~/Downloads/hb_stage_2.xlsx", sheet = "e") #read sheet e

#plot
ggplot(sheet.e,
       aes(x = half_life, y = alpha, color = stage, size = count))+
  geom_point(alpha = 0.7)+
  coord_cartesian(xlim = c(0,50))+
  labs( x = "Half Life", y = "Alpha Life", color = "stage", size = "count")+
  theme_classic()+
  theme(axis.title = element_text(face = "bold"),
        legend.box = "vertical",
        legend.position = "right")



# ==============================================================================
# Task 6 : Reproduce panel 2f: Stacked proportions
# ==============================================================================
f.sheet <- read_excel("~/Downloads/hb_stage_2.xlsx", sheet = "f") #read sheet f
sub.df <- f.sheet[c(1,2,7,8),] ##subset B and Plasma 

#plot
ggplot(sub.df, aes(x = stage, y = proportion, fill = cell_type))+
  geom_bar(stat = "identity", position = "stack", width = 0.7)+
  scale_fill_manual(values = c(Plasma = hb_pal[1], B = hb_pal[4]),
                    breaks = c("Plasma", "B"))+
  coord_cartesian(ylim = c(0,0.3))+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.title = element_blank())

##Stacked Barplot is most useful here for composition data
##Easy to compare total proportion at each timepoint

# ===================================================================
# Task 7: Direct Cell-Cell Interaction Network
# ===================================================================
g.sheet <- read_excel("~/Downloads/hb_stage_2.xlsx", sheet = "g")
g.df <- as.data.frame(g.sheet)
row.names(g.df) <- g.df[[1]]

#remove first column 

g.df <- g.df[-1]

#convert to numeric matrix
mat_g <- as.matrix(g.df)
mode(mat_g) <- "numeric"

#remove diagonal self loops
diag(mat_g) <- 0

#build directed graph
g_graph <- graph_from_adjacency_matrix(mat_g,
                                       mode = "directed",
                                       weighted = T,
                                       diag = F)

#remove zero weighted edges
g_graph <- delete_edges(g_graph, E(g_graph)[E(g_graph)$weight == 0])

# visual characteristcs
E(g_graph)$color = "grey70"
w <- E(g_graph)$weight
E(g_graph)$width = 0.8 +2.0
E(g_graph)$arrow.size = 1

V(g_graph)$color = "pink"
V(g_graph)$frame.color = "grey70"
V(g_graph)$size = 20
V(g_graph)$label.font = 2
V(g_graph)$label.cex = 1.5

V(g_graph)$label.color <- "darkblue"
V(g_graph)$vertex.shape <- "square"

#plot
plot(g_graph,
     main = "Directed cell–cell interaction network")

#W#hy directed? cel singlalling flows in one direction, so directed edge capture
#this natural information flow rather than assuming all interction are reciprocal
##What does edge weight encode biologically?
#Edge weight represent how frequently two cell interacts, high edge weeight 
#intense signalling, stronger molecular binding between cells. 

# ============================================================================
# Task 8 : Final assembly
# ============================================================================
# Arrange all panels into a single figure







