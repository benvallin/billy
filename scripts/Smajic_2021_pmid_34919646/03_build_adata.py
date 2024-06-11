# %% Set up ----

# Import required libraries
import os
import numpy as np
import pandas as pd
import scanpy as sc

# Set working directory
os.chdir("/scratch/ben/rnaseq/")

# %% Build adata ----

adata = sc.read_10x_h5(filename="billy/data/input/Smajic_2021_pmid_34919646/counts.h5")

# %% Add sample/gene metadata ----

# Load sample metadata
sample_metadata = pd.read_csv("billy/data/input/Smajic_2021_pmid_34919646/sample_metadata.csv")
sample_metadata = sample_metadata.set_index(keys="barcode", drop=False)

# Reorder rows to match adata.obs_names
sample_metadata = sample_metadata.loc[adata.obs_names]

# Add sample metadata to adata
adata.obs = sample_metadata

# Set cell type as category and reorder levels
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
adata.obs["cell_type"] = adata.obs["cell_type"].cat.reorder_categories(
['DaNs', 'CADPS2+ neurons', 'Excitatory', 
 'Inhibitory', 'GABA', 'Astrocytes',
 'Ependymal', 'OPCs', 'Oligodendrocytes', 
 'Microglia','Endothelial cells', 'Pericytes']
 )

# Define broad cell type and set as category
adata.obs["broad_cell_type"] = np.where(
  adata.obs["cell_type"].isin(['DaNs', 'CADPS2+ neurons', 'Excitatory', 'Inhibitory', 'GABA']), 
  "Neurons", 
  adata.obs["cell_type"]
  )

adata.obs["broad_cell_type"] = adata.obs["broad_cell_type"].astype("category")
adata.obs["broad_cell_type"] = adata.obs["broad_cell_type"].cat.reorder_categories(
  ['Neurons', 'Astrocytes',
 'Ependymal', 'OPCs', 'Oligodendrocytes', 
 'Microglia','Endothelial cells', 'Pericytes']
  )

# Load gene metadata
gene_metadata = pd.read_csv("ref_seq/feature_metadata/gene_metadata_slim.csv")
gene_metadata = gene_metadata.set_index(keys="ensembl_gene_id", drop=False)

# List ensembl gene IDs from adata that are present/missing in gene metadata
known_ensembl_gene_ids = [i for i in adata.var_names if i in gene_metadata.index]
unknown_ensembl_gene_ids = [i for i in adata.var_names if i not in gene_metadata.index]

# Discard unknown ensembl gene IDs from adata
adata = adata[:, known_ensembl_gene_ids]

# Filter/reorder rows to match adata.var_names
gene_metadata = gene_metadata.loc[adata.var_names]

# Annotate mitochondrial and ribosomal genes
gene_metadata["mt"] = gene_metadata["gene_name"].str.startswith("MT-")
gene_metadata["ribo"] = gene_metadata["gene_name"].str.startswith((("RPS", "RPL")))

# Annotate genes located on chromosomes X, Y and M
gene_metadata["chr_y"] = gene_metadata["chr_name"] == "chrY"
gene_metadata["chr_x"] = gene_metadata["chr_name"] == "chrX"
gene_metadata["chr_m"] = gene_metadata["chr_name"] == "chrM"

# Add gene metadata to adata
adata.var = gene_metadata
adata.var_names = adata.var["ensembl_gene_id_version"]

# Compute QC metrics 
sc.pp.calculate_qc_metrics(adata=adata,
                           expr_type="counts",
                           var_type="genes",
                           qc_vars=["mt", "ribo", "chr_x", "chr_y", "chr_m"],
                           percent_top=None,
                           inplace=True,
                           log1p=False)

# Add XIST counts to observations metadata 
adata.obs["xist_count"] = adata.copy().X[:,adata.var["gene_name"].str.match("XIST")].toarray()

# Write raw adata to disk
adata.copy().write_h5ad(filename="billy/data/input/Smajic_2021_pmid_34919646/raw_adata.h5ad") 

# %%
