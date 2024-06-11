# %% Set up ----

# Import required libraries
import os
import pandas as pd
import scanpy as sc
import decoupler as dc

# Set working directory
os.chdir("/scratch/ben/rnaseq/")

# %% Load annotated filtered adata ----

# Load annotated filtered adata
adata = sc.read_h5ad(filename="/scratch/ben/rnaseq/billy/data/input/Smajic_2021_pmid_34919646/scvi_adata.h5ad")

# Set count matrix to raw counts
adata.X = adata.layers["raw_counts"].copy()

# %% Compute pseudobulks for combinations of donor ID and cell type  ----

# Compute pseudo-bulk profile 
adata_pb = dc.get_pseudobulk(
  adata=adata.copy(),
  sample_col="donor_id",
  groups_col="cell_type",
  layer="raw_counts",
  mode="sum",
  min_cells=0,
  min_counts=0
  )

# Write sample metadata to disk
adata_pb.obs.to_csv("billy/data/input/Smajic_2021_pmid_34919646/sample_metadata_pb.csv", index=False)

# Tidy count matrix for R import
adata_pb = pd.DataFrame(adata_pb.X, index=adata_pb.obs.index.tolist(), columns=adata_pb.var_names.tolist()).reset_index(names="donor_id_cell_type")

# Write count matrix to disk
adata_pb.to_csv("billy/data/input/Smajic_2021_pmid_34919646/raw_count_matrix_pb.csv", index=False)

# %% Compute pseudobulks for combinations of donor ID and broad cell type  ----

# Compute pseudo-bulk profile 
adata_pb = dc.get_pseudobulk(
  adata=adata.copy(),
  sample_col="donor_id",
  groups_col="broad_cell_type",
  layer="raw_counts",
  mode="sum",
  min_cells=0,
  min_counts=0
  )

# Write sample metadata to disk
adata_pb.obs.to_csv("billy/data/input/Smajic_2021_pmid_34919646/sample_metadata_pb_broad_cell_type.csv", index=False)

# Tidy count matrix for R import
adata_pb = pd.DataFrame(adata_pb.X, index=adata_pb.obs.index.tolist(), columns=adata_pb.var_names.tolist()).reset_index(names="donor_id_broad_cell_type")

# Write count matrix to disk
adata_pb.to_csv("billy/data/input/Smajic_2021_pmid_34919646/raw_count_matrix_pb_broad_cell_type.csv", index=False)

# %%
