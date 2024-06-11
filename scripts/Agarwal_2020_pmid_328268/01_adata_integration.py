# %% Set up ----

# Import required libraries
import os
import numpy as np
import pandas as pd
import scanpy as sc

# Set working directory
os.chdir("/scratch/ben/rnaseq/")

# Set path to output figure directory
sc.settings.figdir="/scratch/ben/rnaseq/billy/data/output/Agarwal_2020_pmid_328268/adata_integration/"
output_path = str(sc.settings.figdir)+"/"

# %% Load filtered adata ----

adata = sc.read_h5ad(filename="/scratch/ben/rnaseq/asap/data/input/agarwal_astrocytes/annot_filtered_adata.h5ad") 

# Reorder cell type categories
adata.obs["cell_type"] = adata.obs["cell_type"].cat.reorder_categories(
  ["DaN", "excitatory neuron", "inhibitory neuron", "GABA", "astrocyte", "OPC", "oligodendrocyte", "microglia", "endothelial"]
  )

# Reorder tissue - cell type categories
adata.obs["tissue_cell_type"] = adata.obs["tissue_cell_type"].cat.reorder_categories(
  ['SN - DaN', 'SN - GABA',
   'CTX - excitatory neuron', 'CTX - inhibitory neuron',
   'SN - astrocyte', 'CTX - astrocyte', 
   'SN - OPC', 'CTX - OPC',
   'SN - oligodendrocyte', 'CTX - oligodendrocyte',
   'SN - microglia', 'SN - endothelial']
   )

# Define broad cell type and set as category
adata.obs["broad_cell_type"] = np.where(
  adata.obs["cell_type"].isin(['DaN', 'GABA', 'inhibitory neuron', 'excitatory neuron']), 
  "neuron", 
  adata.obs["cell_type"]
  )

adata.obs["broad_cell_type"] = adata.obs["broad_cell_type"].astype("category")
adata.obs["broad_cell_type"] = adata.obs["broad_cell_type"].cat.reorder_categories(
  ['neuron', 'astrocyte', 'OPC', 'oligodendrocyte', 'microglia', 'endothelial']
  )

# Define tissue - broad cell type category
adata.obs["tissue_broad_cell_type"] = adata.obs["tissue"].astype(str) + " - " + adata.obs["broad_cell_type"].astype(str)
adata.obs["tissue_broad_cell_type"] = adata.obs["tissue_broad_cell_type"].astype("category")
adata.obs["tissue_broad_cell_type"] = adata.obs["tissue_broad_cell_type"].cat.reorder_categories(
  ["SN - astrocyte", "CTX - astrocyte", 
   "SN - OPC", "CTX - OPC",
   "SN - oligodendrocyte", "CTX - oligodendrocyte",
   "SN - neuron", "CTX - neuron", 
   "SN - microglia", "SN - endothelial"]
   )

# %% Plot cell clustering results ----

# Expression level of marker genes on UMAP
leiden_key = "leiden scVI "
sc.pl.umap(adata, 
           neighbors_key="neighbors_scvi",
           color=[leiden_key+"0.4", "tissue_cell_type", 
                  "GFAP", "AQP4", 
                  "VCAN", "OLIG1",
                  "MOG", "MOBP", 
                  "TH", "SLC6A3", 
                  "GAD1", "GAD2",
                  "SATB2", "SLC17A7",
                  "RGS5", "OLR1"],
           title=[leiden_key+"0.4", "tissue - cell type", 
                  "GFAP", "AQP4", 
                  "VCAN", "OLIG1",
                  "MOG", "MOBP", 
                  "TH", "SLC6A3", 
                  "GAD1", "GAD2",
                  "SATB2", "SLC17A7",
                  "RGS5", "OLR1"],
           layer="log_norm_counts", sort_order=True, frameon=False, edges=True, cmap="Reds", ncols=4,
           gene_symbols="gene_name",  vmin=0, vmax="p99",
           legend_fontsize="medium", save="_marker_genes_long.pdf")
# %% 
sc.pl.umap(adata, 
           neighbors_key="neighbors_scvi",
           color=["sample_id", "tissue", "cell_type", "tissue_cell_type", 
                  "GFAP", "VCAN", "MOG", "SNAP25"],
           title=["sample ID", "tissue", "cell type", "tissue - cell type", 
                  "GFAP", "VCAN", "MOG", "SNAP25"],
           layer="log_norm_counts", sort_order=True, frameon=False, edges=True, cmap="Reds", ncols=4,
           gene_symbols="gene_name",  vmin=0, vmax="p99",
           legend_fontsize="medium", save="_marker_genes.pdf", wspace=0.2)

# %% 
# Expression level of TFEB, TFE3, TFEC and MITF on UMAP
# => TFEB and TFE3 only
sc.pl.umap(adata, 
           neighbors_key="neighbors_scvi",
           color=["tissue_cell_type", "tissue_broad_cell_type", "TFEB", "TFE3"],
           title=["tissue - cell type", "tissue - broad cell type", "TFEB", "TFE3"],
           layer="log_norm_counts", sort_order=False, frameon=False, edges=True, cmap="Reds", ncols=5,
           gene_symbols="gene_name",  vmin=0, vmax="p99",
           legend_fontsize="medium", save="_tfeb_tfe3.pdf")

# %% 
# => TFEB, TFE3, TFEC and MITF 
sc.pl.umap(adata, 
           neighbors_key="neighbors_scvi",
           color=["cell_type", "broad_cell_type", "TFEB", "TFE3", "TFEC", "MITF"],
           title=["cell type", "broad cell type", "TFEB", "TFE3", "TFEC", "MITF"],
           layer="log_norm_counts", sort_order=False, frameon=False, edges=True, cmap="Reds", ncols=3,
           gene_symbols="gene_name",  vmin=0, vmax="p99",
           legend_fontsize="medium", save="_tfeb_tfe3_tfec_mitf.pdf")

# Expression level of marker genes on dotplot
# => TFEB and TFE3 only
sc.pl.dotplot(
  adata=adata,
  var_names=["TFEB", "TFE3"],
  groupby="cell_type",
  gene_symbols="gene_name",
  layer="log_norm_counts",
  var_group_rotation=0,
  mean_only_expressed=True,
  save="tfeb_tfe3.pdf"
)

sc.pl.dotplot(
  adata=adata,
  var_names=["TFEB", "TFE3"],
  groupby="broad_cell_type",
  gene_symbols="gene_name",
  layer="log_norm_counts",
  var_group_rotation=0,
  mean_only_expressed=True,
  save="tfeb_tfe3_broad.pdf"
)

# => TFEB, TFE3, TFEC and MITF 
sc.pl.dotplot(
  adata=adata,
  var_names=["TFEB", "TFE3", "TFEC", "MITF"],
  groupby="cell_type",
  gene_symbols="gene_name",
  layer="log_norm_counts",
  var_group_rotation=0,
  mean_only_expressed=True,
  save="tfeb_tfe3_tfec_mitf.pdf"
)

sc.pl.dotplot(
  adata=adata,
  var_names=["TFEB", "TFE3", "TFEC", "MITF"],
  groupby="broad_cell_type",
  gene_symbols="gene_name",
  layer="log_norm_counts",
  var_group_rotation=0,
  mean_only_expressed=True,
  save="tfeb_tfe3_tfec_mitf_broad.pdf"
)

# %%
