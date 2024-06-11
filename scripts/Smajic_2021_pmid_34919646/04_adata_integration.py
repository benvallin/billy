# %% Set up ----

# Import required libraries
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch

# Set working directory
os.chdir("/scratch/ben/rnaseq/")

# Set path to output figure directory
sc.settings.figdir="/scratch/ben/rnaseq/billy/data/output/Smajic_2021_pmid_34919646/adata_integration/"
output_path = str(sc.settings.figdir)+"/"

torch.set_float32_matmul_precision('high')

# %% Load filtered adata ----

# Define batch variable
batch_key = "donor_id"

# Define method for HVG detection
hvg_detection_method = "seurat_v3"
# hvg_detection_method = "seurat"

# Load filtered adatasn
adata = sc.read_h5ad(filename="billy/data/input/Smajic_2021_pmid_34919646/raw_adata.h5ad")

# Store raw counts 
adata.layers["raw_counts"] = adata.X.copy()

# %% Normalize and log transform gene counts ----

# Perform normalization (to 10,000 reads per cell)
sc.pp.normalize_total(adata=adata,
                      target_sum=1e4,
                      inplace=True)

# Perform log transformation (using natural logarithm)
sc.pp.log1p(adata)

# Store normalized and log transformed counts matrix in norm_log_counts
adata.layers["log_norm_counts"] = adata.X.copy()

# %% scVI integration and clustering ----

# Set count matrix to log normalised counts 
adata.X = adata.layers["log_norm_counts"].copy()

# Find HVG on each batch separately
if hvg_detection_method == "seurat_v3":
  sc.pp.highly_variable_genes(adata=adata, flavor="seurat_v3", layer="raw_counts", n_top_genes=2000,
                              batch_key=batch_key)
elif hvg_detection_method == "seurat":
  sc.pp.highly_variable_genes(adata=adata, batch_key=batch_key)

# Prepare adata for scVI integration
# => Create a copy of adata containing only the HVG
adata_hvg = adata[:, adata.var["highly_variable"]].copy()

# => Set up the HVG-only adata
scvi.model.SCVI.setup_anndata(adata=adata_hvg, 
                              layer="raw_counts", 
                              batch_key=batch_key)

# scvi.model.SCVI.setup_anndata(adata=adata_hvg, 
#                               layer="raw_counts", 
#                               batch_key=batch_key,
#                               continuous_covariate_keys=["total_counts", "pct_counts_ribo"])

# # Build scVI model
# model_scvi = scvi.model.SCVI(adata=adata_hvg)
# model_scvi.view_anndata_setup()

# # Train scVI model
# model_scvi.train()

# # Write ref_model_scvi to disk
# model_scvi.save(dir_path="billy/data/input/Smajic_2021_pmid_34919646/scvi_models/scvi_model", overwrite=True)

# Load ref_model_scvi from disk
model_scvi = scvi.model.SCVI.load(dir_path="billy/data/input/Smajic_2021_pmid_34919646/scvi_models/scvi_model", 
                                  adata=adata_hvg)

# Pass the scVI-normalised counts to the HVG-only adata
adata_hvg.layers["norm_counts_scvi"] = model_scvi.get_normalized_expression(library_size=1e4)

# Pass the scVI-corrected embedding to the HVG-only adata
adata_hvg.obsm["X_scvi"] = model_scvi.get_latent_representation()
adata_hvg.obsm["X_scvi_mde"] = scvi.model.utils.mde(adata_hvg.obsm["X_scvi"])

# Transfer the scVI-corrected embedding to the full adata
for i in ["X_scvi", "X_scvi_mde"]:
  adata.obsm[i] = adata_hvg.obsm[i]

# Compute neighbors using the scVI-corrected embedding
sc.pp.neighbors(adata, use_rep = "X_scvi", key_added="neighbors_scvi")

# Compute UMAP
sc.tl.umap(adata, neighbors_key="neighbors_scvi")

# Find clusters
leiden_key = "leiden scVI "
for i in np.round(np.arange(0.1, 1, 0.1), 1):
  sc.tl.leiden(adata, 
               neighbors_key="neighbors_scvi",
               key_added=leiden_key+str(i),
               resolution=i)

# Write annotated filtered adata to disk
adata.copy().write_h5ad(filename="billy/data/input/Smajic_2021_pmid_34919646/scvi_adata.h5ad") 

# %% Plot cell clustering results ----

# Expression level of marker genes on UMAP
sc.pl.umap(adata, 
           neighbors_key="neighbors_scvi",
           color=["cell_type", "donor_id", "MAP2", "TH", "SLC17A6", 
                  "GAD2", "GRIK1", "CADPS2", "VCAN", "MOBP", 
                   "AQP4", "FOXJ1", "CLDN5", "PDGFRB", "CD74"],
           title=["cell type", "donor id", "MAP2", "TH", "SLC17A6", 
                  "GAD2", "GRIK1", "CADPS2", "VCAN", "MOBP", 
                   "AQP4", "FOXJ1", "CLDN5", "PDGFRB", "CD74"],
           layer="log_norm_counts", sort_order=False, frameon=False, edges=True, cmap="Reds", ncols=5,
           gene_symbols="gene_name",  vmin=0, vmax="p99",
           legend_fontsize="medium", save="_marker_genes.pdf")

# Expression level of TFEB, TFE3, TFEC and MITF on UMAP
# => TFEB and TFE3 only
sc.pl.umap(adata, 
           neighbors_key="neighbors_scvi",
           color=["cell_type", "broad_cell_type", "TFEB", "TFE3"],
           title=["cell type", "broad cell type", "TFEB", "TFE3"],
           layer="log_norm_counts", sort_order=False, frameon=False, edges=True, cmap="Reds", ncols=5,
           gene_symbols="gene_name",  vmin=0, vmax="p99",
           legend_fontsize="medium", save="_tfeb_tfe3.pdf")

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
