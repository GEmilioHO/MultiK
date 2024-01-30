# MultiK
## NOTE: This updated version fixes bugs and allows processing of SCT and LogNorm normalized single and integrated Seurat objects.
Single cell RNA-seq analysis tool that objectively selects multiple insightful numbers of clusters (*K*).


## Getting started


## Installing
**MultiK** relies on the following R packages: **Seurat**, **SigClust**, please install them before installing **MultiK**.

```{}
#install.packages("Seurat")
#install.packages("SigClust")
```

```{}
#install.packages("devtools")
#library(devtools)
#install_github("GEmilioHO/MultiK")
```

**MultiK( )** is the main function that implements the subsampling and application of the Seurat clustering over multiple resolution parameters. Details about this function can be found in the user manual:

```{}
?MultiK
```


The main function, **MultiK( )**, takes a Seurat object with the normalized expression matrix and other parameters set by default values if not specified. MultiK explores a range of resolution parameters (from 0.05 to 2.00 with a step size of 0.05) in Seurat clustering, and aggregates all the clustering runs that give rise to the same K groups regardless of the resolution parameter and computes a consensus matrix for each K.

Note: MultiK re-selects highly variable genes in each subsampling run. Also, MultiK, by default, uses 30 principal components and 20 K-nearest neighbors in Seurat clustering.  



## MultiK workflow

## Example
Prior to running MultiK, data should be normalised and scaled.
```{}
library(Seurat)
library(MultiK)

seurat.obj <- readRDS("path/to/seurat_object.rds")

# Perform LogNormalisation
seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize",
                            scale.factor = 10000)
seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
seurat.obj <- ScaleData(seurat.obj, features = rownames(seurat.obj))

# Alternatively, perform SCT normalisation
seurat.obj <- SCTransform(seurat.obj, vst.flavor = "v2", verbose = TRUE)

# Set DefaultAssay to the desired assay to use.
DefaultAssay(seurat.obj) <- "SCT" # or DefaultAssay(seurat.obj) <- "RNA"

# If running for an integrated object:
DefaultAssay(seurat.obj) <- "integrated"
```
### Step 1: Run **MultiK** main algorithm to determine optimal Ks

Run subsampling and consensusing clustering to generate output for evaluation (this step can take a long time). For demonstration purpose, we are running 10 reps here. For real data pratice, we recommend using at least 100 reps. Set nPC to the number of principal components used for running PCA and finding neighbours.
```{}
# For RNA and SCT assays:
multik <- MultiK(seurat.obj, nPC = 30, reps = 10)

# For integrated assays: 
# batch = variable for identifying different batches;
# integrated.assay.norm.method = normalisation method used for each batch
#                                prior to integration (LogNorm or SCT)
multik <- MultiK(seurat.obj, nPC = 30, reps = 10, batch = "orig.ident",
                 integrated.assay.norm.method = "LogNorm")
```

Make MultiK diagnostic plots: 
```{}
DiagMultiKPlot(multik$k, multik$consensus)
```

### Step 2: Assign _classes_ and _subclasses_

Get the clustering labels at optimal K level (using either the upper or lower K). Set nPC to the number of principal components used for running PCA and finding neighbours.
```{}
clusters <- getClusters(seurat.obj, optK = 3, nPC = 30)
```

Run SigClust at optimal K level:
```{}
pval <- CalcSigClust(seurat.obj, clusters$clusters)
```

Make diagnostic plots (this includes a dendrogram of cluster centroids with the pairwise SigClust _p_ values mapped on the nodes, and a heatmap of the pairwise SigClust _p_ values)
```{}
PlotSigClust(seurat.obj, clusters$clusters, pval)
```

Add cluster labels to Seurat Object.
```{}
seurat.obj$multiK.clusters <- clusters$clusters
```

## Contact
If you have any questions, please contact: **Siyao Liu** (<siyao@email.unc.edu>) (original developer) </br>
Alternatively, contact: **Gabriel Herrera-Oropeza** (<gabriel.herrera_oropeza@kcl.ac.uk>) (SCT and integration modifications)
