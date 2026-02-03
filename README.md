## SpaNiche: spatial niche analysis to explore colocalization patterns and cellular interactions in spatial transcriptomics data

We propose a computational framework for spatial niche analysis
(SpaNiche) in spatial transcriptomics data to uncover colocalization
patterns and infer potential ligand-receptor interactions. SpaNiche
leverages graph-regularized joint non-negative matrix factorization to
integrate information from cell abundance and ligand-receptor
expression, identifying spatial colocalization patterns among cell types
while providing insights into associated ligand-receptor interactions.
In addition, SpaNiche employs consensus clustering to define “ecotypes”,
enhancing its utility in multi-sample spatial transcriptomics datasets.

<figure>
<img src="./tutorial/example/figure_1_framework.png"
alt="Schematic overview of SpaNiche" />
<figcaption aria-hidden="true">Schematic overview of
SpaNiche</figcaption>
</figure>

## Installation

SpaNiche has been tested on both Windows (R version 4.2.1) and
Linux-based systems (R version 4.2.3), and should be compatible with
other R environments.

You can install SpaNiche from GitHub using:

    devtools::install_github("SiyuanHuang1/SpaNiche")

## Application examples

-   [1a. Graph-regularized integrated NMF in single-sample spatial
    transcriptomics data using SpaNiche](./tutorial/1a.singlesample_iNMF.md)

This analysis typically reveals colocalization patterns of different
cell type combinations within a spatial transcriptomics slide, along
with their corresponding ligand-receptor interaction features. Because
it jointly factorizes cell type abundance and ligand-receptor matrices,
it often yields better correspondence between patterns and interactions.
Moreover, the workflow is more streamlined compared to two-step
approaches (e.g., spatial clustering followed by CellphoneDB/cellchat
analysis).
