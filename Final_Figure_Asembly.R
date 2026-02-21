# ============================================================================
# Task 8 : Final assembly
# ============================================================================
# Arrange all panels into a single figure
#capture base R and igrapgh as grobs 
#ggplot need to be saved as grob so grid.arrange won't run into error

hm_grob <- grid.grabExpr(
  pheatmap(mat = Gene_counts[,2:7],
           border_color = "Black",
           legend = T,
           fontsize_row = 7,
           fontsize_col = 7,
           labels_row = Gene_counts[,1],  #check this line
           cluster_rows = T,
           cluster_cols = T,
           color = blues9)
)

###############################################################################
#recreate ggplot as objects
vol_plot <- ggplot(DEG,
      aes(x = log2FoldChange, y = X.log10PAdj, 
          color = as.factor(significance)))+
  geom_point(alpha = 0.8, size = 2)+
  geom_vline(xintercept = c(-2, 1.8), linetype = "dashed")+
  geom_hline(yintercept = 0.5, linetype = "dashed")+
  scale_color_manual(values = vol_colors)+
  theme_classic()+
  labs(
    title = "Volcano plot: log2FoldChange vs log10(Padj)",
    x = "log2Foldchange",
    y = "log10PAdj",
    color = "Significance"
  )+
  theme(axis.title = element_text(family = "arial"))
vol_grob <- ggplotGrob(vol_plot) #capture as grob

###############################################################################scatter_p <- ggplot(breast_cancer, #assign as objects
       aes(x = radius_mean, y = texture_mean, color = as.factor(diagnosis)))+
  geom_point(size = 2, alpha = 1)+
  coord_cartesian(xlim = c(0,25), ylim = c(0,40))+
  labs( x = "radius_mean",
        y = "texture_mean",
        color = "diagnosis",
        title = "Texture_mean Vs Radius_mean")+
  theme_classic()

scatter_p_grob <- ggplotGrob(scatter_p) #capture as grob

###############################################################################
cor_pheatmap <- ggplot(cor_df,
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
cor_pheatmap_grob <- ggplotGrob(cor_pheatmap) #capture as grob

###############################################################################
scatter_p2 <- ggplot(breast_cancer,
                     aes(x = smoothness_mean, y = compactness_mean, colour = as.factor(diagnosis)))+
  geom_point(size = 2, alpha = 0.8)+
  labs(x = "smoothness_mean",
       y = "compactness_mean",
       title = "Smoothness_mean vs Compactness_mean",
       color = "diagnosis")+
  theme_classic()
scatter_p2_grob <- ggplotGrob(scatter_p2) #capture as grob

###############################################################################
density_p <- ggplot(breast_cancer,
                    aes(x = area_mean, fill = as.factor(diagnosis)))+
  geom_density(alpha = 0.5)+
  labs(x = "area_mean",
       y = "Density",
       fill = "Diagnosis",
       title = "Area Distribution")+
  theme_classic()
densit_grob <- ggplotGrob(density_p) #capture as grob

###############################################################################
## final aseembly

part1_panels <- grid.arrange(
  scatter_p2_grob,
  hm_grob,
  vol_grob,
  cor_pheatmap_grob,
  scatter_p_grob,
  densit_grob,
  ncol = 2
)

## save final figure 
ggsave(
  "final_multi_figure_part1&2.png",
  part1_panels,
  width = 10,
  height = 10,
  dpi = 300
)

#Assembly part 3
panel_3a <- ggplot(a.sheet,
                       aes(x = cell_type , y = new_ratio ))+
  geom_boxplot(fill = hb_pal[c(1:10)], alpha = 1)+
  labs(x = "Cell_type",y = "Ratio",title = "Cell-Type Ratio Distribution",
  )+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

panel_3a_grob <- ggplotGrob(panel_3a)

###############################################################################

panel_3b <- ggplot(b.sheet, aes(log2(half_life), log2(alpha))) +
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

panel_3b_grob <- ggplotGrob(panel_3b)

###############################################################################

panel_3c_grob <- grid.grabExpr(
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
)

###############################################################################

panel_3d <- ggplot(d.sheet_df1,
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

panel_3d_grob <- ggplotGrob(panel_3d)

###############################################################################

panel_3e <- ggplot(sheet.e,
                   aes(x = half_life, y = alpha, color = stage, size = count))+
  geom_point(alpha = 0.7)+
  coord_cartesian(xlim = c(0,50))+
  labs( x = "Half Life", y = "Alpha Life", color = "stage", size = "count")+
  theme_classic()+
  theme(axis.title = element_text(face = "bold"),
        legend.box = "vertical",
        legend.position = "right")

panel_3e_grob <- ggplotGrob(panel_3e)

###############################################################################

panel_3f <- ggplot(sub.df, aes(x = stage, y = proportion, fill = cell_type))+
  geom_bar(stat = "identity", position = "stack", width = 0.7)+
  scale_fill_manual(values = c(Plasma = hb_pal[1], B = hb_pal[4]),
                    breaks = c("Plasma", "B"))+
  coord_cartesian(ylim = c(0,0.3))+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.title = element_blank())

panel_3f_grob <- ggplotGrob(panel_3f)

###############################################################################

grid.newpage()
plot(g_graph,
     main = "Directed cell–cell interaction network")
grid.echo()
panel_3g_grob <- grid.grab()
class(panel_3g_grob)
#arrnage 

part3_panels <- grid.arrange(
  panel_3a_grob,
  panel_3b_grob,
  panel_3e_grob,
  panel_3f_grob,
  ncol = 4
)

#save output
ggsave("final_multi_panel3.png",
       part3_panels,
       height = 10,
       width = 15,
       dpi = 300
       )

part3cdg <- grid.arrange(
  panel_3d_grob,
  panel_3g_grob,
  panel_3c_grob,
  ncol = 3
)

ggsave("final_multi_panel3cdg.png",
       part3cdg,
       height = 10,
       width = 20,
       dpi = 300)

