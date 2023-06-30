library(garnett)
library(monocle3)

## Step 0: read your data
matdata <- read_csv('path_to_scRNAseq_count_matrix')
mat_cell_metadata <- read_csv('path_to_cell_metadata')
mat_genes_metadata <- read_csv('path_to_genes_metadata')

cds <- new_cell_data_set(matdat,
                         cell_metadata = mat_cell_metadata,
                         gene_metadata = mat_genes_metadata
)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 50, build_nn_index = TRUE, verbose = TRUE,
                      scaling = TRUE, norm_method = "log")
plot_pc_variance_explained(cds)

## Step 2: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds, reduction_method = "UMAP", 
                        cores = 16, max_components = 2, 
                        umap.n_neighbors = 15, preprocess_method = "PCA", 
                        build_nn_index = TRUE
                        )
## Step 3: Cluster the cells
cds <- cluster_cells(cds, 
                     reduction_method = 'UMAP', partition_qval = 0.01, cluster_method = "louvain")
## Step 4: Learn a graph
cds <- learn_graph(cds, 
                   verbose = TRUE, use_partition = FALSE
)

## Step 6: Order cells
cds <- order_cells(cds, reduction_method = "UMAP")

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="CD56bright"){
  cell_ids <- which(colData(cds)[, "Clusters.ident"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

## subset cells by branch
cds_sub <- choose_graph_segments(cds)

##work in 3D
cds_3d <- reduce_dimension(cds, max_components = 3, preprocess_method = "PCA", 
                           reduction_method = "UMAP")
cds_3d <- cluster_cells(cds_3d, cluster_method = "louvain")
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="Clusters")


gene_fits <- fit_models(cds, model_formula_str = "~Clusters", 
                        cores = 15)
fit_coefs <- coefficient_table(gene_fits)
emb_time_terms <- fit_coefs %>% filter(term == "Clusters")
#emb_time_terms <- emb_time_terms %>% mutate(q_value = p.adjust(p_value))
#sig_genes <- emb_time_terms %>% filter (q_value < 0.05) %>% pull(gene_short_name)

# With graph autocorrelation:
plot_cells(cds, genes=sig_genes[1:4])

marker_test_res <- top_markers(cds, group_cells_by="Clusters", 
                               reference_cells=10000, cores=16)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="Clusters",
                    ordering_type="maximal_on_diag",
                    max.size=3)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.05) %>%
  group_by(cell_group) %>%
  top_n(3, wt=specificity)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="Clusters",
                    ordering_type="cluster_row_col",
                    max.size=3)

AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                       colData(cds)$cell.type %in% c("AFD")]
plot_genes_in_pseudotime(cds,
                         color_cells_by="Clusters",
                         min_expr=0.5)
