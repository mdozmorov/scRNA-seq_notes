# scRNA-seq data analysis tools and papers

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com)

Single-cell RNA-seq related tools and genomics data analysis resources. Tools are sorted by publication date, reviews and most recent publications on top. Unpublished tools are listed at the end of each section. Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes. See [scATAC-seq_notes](https://github.com/mdozmorov/scATAC-seq_notes) for scATAC-seq related resources.

# Table of content

<!-- START doctoc generated TOC please keep comment here to allow auto update -->

<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Awesome](#awesome)
- [Courses](#courses)
  - [Tutorials](#tutorials)
- [Preprocessing pipelines](#preprocessing-pipelines)
  - [Format conversion](#format-conversion)
  - [Visualization pipelines](#visualization-pipelines)
- [Quality control](#quality-control)
- [Normalization](#normalization)
- [Integration, Batch correction](#integration-batch-correction)
- [Imputation](#imputation)
- [Dimensionality reduction](#dimensionality-reduction)
- [Clustering](#clustering)
  - [Time, trajectory inference](#time-trajectory-inference)
  - [Networks](#networks)
  - [RNA velocity](#rna-velocity)
- [Differential expression](#differential-expression)
  - [Differential abundance](#differential-abundance)
- [Downstream analysis]($downstream-analysis)
- [CNV](#cnv)
- [Splicing](#splicing)
- [Annotation, subpopulation identification](#annotation-subpopulation-identification)
  - [Cell markers](#cell-markers)
    - [Immune markers](#immune-markers)
    - [Brain markers](#brain-markers)
- [Immuno-analysis](#immuno-analysis)
- [Cell-cell interactions](#cell-cell-interactions)
- [Phylogenetic inference](#phylogenetic-inference)
- [Simulation](#simulation)
  - [Power](#power)
  - [Benchmarking](#benchmarking)
- [Deep learning](#deep-learning)
- [Spatial transcriptomics](#spatial-transcriptomics)
- [Technology](#technology)
  - [10X Genomics](#10x-genomics)
    [10X QC](#10x-qc)
- [Data](#data)
  - [Human](#human)
    [Cancer](#cancer)
  - [Mouse](#mouse)
  - [Brain](#brain)
- [Links](#links)
  - [Papers](#papers)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Awesome

- Review of scRNA-seq cell type annotation methods. Table 1 - tools grouped by methods, each described in the text. Tale 2 - annotation databases, CellMaker and PangaloDB are the largest. Table 3 - tools and cell types used for benchmarking. <details>
    <summary>Paper</summary>
    Pasquini, Giovanni. “Automated Methods for Cell Type Annotation on scRNA-Seq Data.” Computational and Structural Biotechnology Journal, 2021, 9. https://doi.org/10.1016/j.csbj.2021.01.015
</details>

- Overview of various steps and tools for data analysis of single cell transcriptomics (scRNA-seq), chromatin accessibility (scATAC-seq), surface protein (CITE-seq), antigen immune receptor repertoire (AIRR, TCR and BCR profiling), and spatial transcriptomics. QC, doublet removal ([scDoubletFinder](https://bioconductor.org/packages/scDblFinder/)), normalization ([Scran](https://bioconductor.org/packages/scran/)), batch removal ([Harmony](https://portals.broadinstitute.org/harmony/articles/quickstart.html) and others), cell cycle removal ([Tricycle](https://bioconductor.org/packages/tricycle/)), cell annotation, trajectory analysis ([dynguidelines](https://zouter.shinyapps.io/server/)), differential expression, gene set enrichment analysis, cell-cell communication. Data integration methods. [Glossary](https://www.nature.com/articles/s41576-023-00586-w#glossary). Extended version in the [Single-cell best practices](https://www.sc-best-practices.org/preamble.html) book. <details>
    <summary>Paper</summary>
    Heumos, Lukas, Anna C. Schaar, Christopher Lance, Anastasia Litinetskaya, Felix Drost, Luke Zappia, Malte D. Lücken, et al. “Best Practices for Single-Cell Analysis across Modalities.” Nature Reviews Genetics, March 31, 2023. https://doi.org/10.1038/s41576-023-00586-w.
</details>

- [Community-curated list of software packages and data resources for single-cell, including RNA-seq, ATAC-seq, etc.](https://github.com/seandavi/awesome-single-cell) by Sean Davis

- [Notes on scATAC-seq analysis](https://github.com/crazyhottommy/scATACseq-analysis-notes) by Ming Tang

- [Dave Tand's blog, single cell posts](https://davetang.org/muse/category/single-cell-2/)

- [www.scrna-tools.org](https://www.scrna-tools.org/) - The scRNA-tools database <details>
    <summary>Paper</summary>
    Zappia, Luke, and Fabian J. Theis. “Over 1000 Tools Reveal Trends in the Single-Cell RNA-Seq Analysis Landscape.” Genome Biology 22, no. 1 (December 2021): 301. https://doi.org/10.1186/s13059-021-02519-4.
    Zappia, Luke, Belinda Phipson, and Alicia Oshlack. "[Exploring the single-cell RNA-seq analysis landscape with the scRNA-tools database" https://doi.org/10.1371/journal.pcbi.1006245   PLoS computational biology, June 25, 2018
</details>

- [awesome-10x-genomics](https://github.com/johandahlberg/awesome-10x-genomics) - List of tools and resources related to the 10x Genomics GEMCode/Chromium system

- [awesome-deep-learning-single-cell-papers](https://github.com/OmicsML/awesome-deep-learning-single-cell-papers) - categorized list of latest scRNA-seq papers using deep learning

- [awesome-vdj](https://github.com/slowkow/awesome-vdj) - Tools and databases for analyzing HLA and VDJ genes, by [slowkow](https://github.com/slowkow)

- [SingleCell Omics](https://docs.google.com/spreadsheets/d/1IPe2ozb1Mny8sLvJaSE57RJr3oruiBoSudAVhSH-O8M/edit#gid=11468010) - A Google Doc with a structured collection of scRNA-seq methods, software, and many other scRNA-seq information, by @albertvilella

- [R_packages_for_scRNA-seq.pdf](R_packages_for_scRNA-seq.pdf) - Bioconductor software packages for single-cell analysis. <details>
    <summary>Paper</summary>
    Amezquita, Robert A., Aaron T. L. Lun, Etienne Becht, Vince J. Carey, Lindsay N. Carpp, Ludwig Geistlinger, Federico Martini, et al. "Orchestrating Single-Cell Analysis with Bioconductor" https://doi.org/10.1038/s41592-019-0654-x  Nature Methods, December 2, 2019.
</details>

- [single-cell-pseudotime](https://github.com/agitter/single-cell-pseudotime) - an overview of single-cell RNA-seq pseudotime estimation algorithms, comprehensive collection of links to software and accompanying papers, by Anthony Gitter


## Courses

- [Single-cell best practices](https://www.sc-best-practices.org/) - scRNA-seq analysis with Python/scanpy in. From preprocessing to trajectory inference, network analysis, immune receptor profiling, and more. [GitHub](https://github.com/theislab/single-cell-best-practices) with Jupyter notebooks/code.

- Review of single-cell transcriptomics technologies and analysis steps and software. Sample preparation, scRNA-seq preprocessing, QC, normalization, batch correction, dimensionaliry reduction. Downstream analysis on cell level (clustering, trajectory inference), gene level (differential expression, functional enrichment, network analysis). Table 1 - preprocessing pipelines and tools, brief description. Table 2 - clustering algorithms. <details>
    <summary>Paper</summary>
    Nayak, Richa, and Yasha Hasija. “A Hitchhiker’s Guide to Single-Cell Transcriptomics and Data Analysis Pipelines.” Genomics 113, no. 2 (March 2021): 606–19. https://doi.org/10.1016/j.ygeno.2021.01.007.
</details>

- [Orchestrating Single-Cell Analysis with Bioconductor](http://bioconductor.org/books/release/OSCA/) - scRNA-seq analysis overview within Bioconductor ecosystem, bookdown. SingleCellexperiment, scran and scater examples. Table S1 - summary of packages for data input, infrastructure, QC, integration, dimensionality reduction, clustering, pseudotime, differential expression, functional enrichment, simulation, benchmarking data, and data packages. Types of feature selection. Associated GitHub repos: [OrchestratingSingleCellAnalysis](https://github.com/Bioconductor/OrchestratingSingleCellAnalysis) [OSCABase](https://github.com/Bioconductor/OSCABase), [OrchestratingSingleCellAnalysis](https://github.com/seandavi/OrchestratingSingleCellAnalysis). [simpleSingleCell](https://bioconductor.org/packages/simpleSingleCell/) R package, a step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor, by Aaron Lun et al., [rendered version](http://bioinformatics.age.mpg.de/presentations-tutorials/presentations/modules/single-cell//bioconductor_tutorial.html). <details>
    <summary>Paper</summary>
    Amezquita RA, Lun AT, Becht E, Carey VJ, Carpp LN, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L. "Orchestrating single-cell analysis with Bioconductor" https://doi.org/10.1038/s41592-019-0654-x Nature methods. 2020 Feb;17(2):137-45.
</details>

- [Analysis of single cell RNA-seq data, www.singlecellcourse.org](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html) - step-by-step scRNA-seq analysis course. R-based, with code examples, explanations, exercises. From alignment (STAR) and QC (FASTQC) to introduction to R, SingleCellExperiment class, `scater` object, data exploration (reads, UMI), filtering, normalization (`scran`), batch effect removal (`RUV`, `ComBat`, `mnnCorrect`, GLM, `Harmony`), clustering and marker gene identification (`SINCERA`, `SC3`, tSNE, `Seurat`), feature selection (`M3Drop::M3DropConvertData`, `BrenneckeGetVariableGenes`), pseudotime analysis (`TSCAN`, `Monocle`, diffusion maps, `SLICER`, `Ouija`, `destiny`), imputation (`scImpute`, `DrImpute`, `MAGIC`), differential expression (Kolmogorov-Smirnov, Wilcoxon, `edgeR`, `Monocle`, `MAST`), data integration (`scmap`, cell-to-cell mapping, `Metaneighbour`, `mnnCorrect`, `Seurat`'s canonical correllation analysis). Search for scRNA-seq data ([scfind](https://github.com/hemberg-lab/scfind) R package), as well as [Hemberg group’s public datasets](https://hemberg-lab.github.io/scRNA.seq.datasets/). [Seurat chapter](https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html). ["Ideal" scRNA-seq pipeline](https://scrnaseq-course.cog.sanger.ac.uk/website/ideal-scrnaseq-pipeline-as-of-oct-2017.html). [Video lectures](https://www.youtube.com/watch?list=PLEyKDyF1qdOYAhwU71qlrOXYsYHtyIu8n&v=56n77bpjiKo). <details>
    <summary>Paper</summary>
    Andrews, Tallulah S., Vladimir Yu Kiselev, Davis McCarthy, and Martin Hemberg. "Tutorial: Guidelines for the Computational Analysis of Single-Cell RNA Sequencing Data." https://doi.org/10.1038/s41596-020-00409-w Nature Protocols, December 7, 2020. 
</details>

- [ANALYSIS OF SINGLE CELL RNA-SEQ DATA](https://broadinstitute.github.io/2019_scWorkshop/index.html) course by Orr Ashenberg, Dana Silverbush, Kirk Gosik

### Tutorials

- [Single-cell RNA-seq data analysis workshop](https://github.com/hbctraining/scRNA-seq_online) by the Harvard Chan Bioinformatics Core. [Lessons](https://github.com/hbctraining/scRNA-seq_online/tree/master/lessons) - hands-on Introduction to single-cell RNA-seq analysis using Seurat/RStudio, starting from count matrices. 

- [Rockfeller University scRNAseq tutorial](https://rockefelleruniversity.github.io/scRNA-seq/)

- [ASAP: Full pipeline on a project imported from the Human Cell Atlas](https://asap.epfl.ch/home/tutorial?t=full_pipeline)

- [Introduction to single-cell RNA-seq technologies](https://figshare.com/articles/Introduction_to_single-cell_RNA-seq_technologies/7704659/1), presentation by Lior Pachter. Key figures, references, statistics. [Slides](https://figshare.com/articles/Introduction_to_single-cell_RNA-seq_technologies/7704659/1), and [notes](https://liorpachter.wordpress.com/2019/02/19/introduction-to-single-cell-rna-seq-technologies/)

- [Machine learning for single cell analysis workshop](https://www.krishnaswamylab.org/workshop). Presentations (Google Slides) and Jupyter notebooks run on Google Colab

- [Pbmc dataset Roche analysis](https://almutlue.github.io/batch_dataset/pbmc_roche.html#)

- [Preprocessing and normalization of single-cell RNA-seq droplet data](https://kkorthauer.org/fungeno2019/singlecell/vignettes/1.2-preprocess-droplet.html)


## Preprocessing pipelines

- Assessment of 9 preprocessing pipelines (Cell Ranger, Optimus, salmon alevin, kallisto bustools, dropSeqPipe, scPipe, zUMIs, celseq2 and scruff) on 10X and CEL-Seq2 datasets ([scmixology](https://github.com/LuyiTian/sc_mixology) and others, 9 datasets total). All pipelines coupled with performant post-processing (normalization, filtering, etc.) produce comparable data quality in terms of clustering/agreement with known cell types. Low-expressed genes are discordant. Details and specific results of each pipeline. [GitHub with pre-/postprocessing scripts](https://github.com/YOU-k/preprocess) <details>
    <summary>Paper</summary>
    You, Yue, Luyi Tian, Shian Su, Xueyi Dong, Jafar S Jabbari, Peter F Hickey, and Matthew E Ritchie. "Benchmarking UMI-Based Single Cell RNA-Sequencing Preprocessing Workflows" https://doi.org/10.1186/s13059-021-02552-3  Genome Biology. 14 December 2021
</details>

- [Single cell current best practices tutorial, GitHub](https://github.com/theislab/single-cell-tutorial). QC (count depth, number of genes, % mitochondrial), normalization (global, downsampling, nonlinear), data correction (batch, denoising, imputation), feature selection, dimensionality reduction (PCA, diffusion maps, tSNE, UMAP), visualization, clustering (k-means, graph/community detection), annotation, trajectory inference (PAGA, Monocle), differential analysis (DESeq2, EdgeR, MAST), gene regulatory networks. Description of the bigger picture at each step, latest tools, their brief description, references. R-based Scater as the full pipeline for QC and preprocessing, Seurat for downstream analysis, scanpy Python pipeline. Links and refs to other tutorials. <details>
    <summary>Paper</summary>
    Luecken, Malte D., and Fabian J. Theis. "Current Best Practices in Single-Cell RNA-Seq Analysis: A Tutorial" https://doi.org/10.15252/msb.20188746  Molecular Systems Biology 15, no. 6 (June 19, 2019)
</details>

- [Alevin](https://github.com/COMBINE-lab/salmon) - end-to-end droplet-based scRNA-seq (10X Genomics) processing pipeline performing cell barcode detection (two-step whitelisting procedure), read mapping, UMI deduplication (parsimonious UMI graphs, PUGs), resolving multimapped reads (EM method to resolve UMI collisions), gene count estimation. Intelligently handles UMI deduplication and multimapped reads, resulting in more accurate gene abundance estimation. Input - sample-demultiplexed FASTQ, output - gene-level UMI counts. Compared against the Cell Ranger, dropEst, STAR and featureCount-based pipelines, UMI-tools, alevin is more accurate and quantifies a greater proportion of sequenced data, especially on combined genomes. Approx. 21X faster than Cell Ranger, low memory requirements, 10-12 threads optimal. C++ implementation, part of [Salmon](https://github.com/COMBINE-lab/salmon). [Alevin documentation](https://salmon.readthedocs.io/en/latest/alevin.html), [Tutorials](https://combine-lab.github.io/alevin-tutorial/#blog) that include visualization options. <details>
    <summary>Paper</summary>
    Srivastava, Avi, Laraib Malik, Tom Smith, Ian Sudbery, and Rob Patro. "Alevin Efficiently Estimates Accurate Gene Abundances from DscRNA-Seq Data" https://doi.org/10.1186/s13059-019-1670-y  Genome Biology, (December 2019)
</details>

- [bigSCale](https://github.com/iaconogi/bigSCale) - scalable analytical framework to analyze large scRNA-seq datasets, UMIs or counts. Pre-clustering, convolution into iCells, final clustering, differential expression, biomarkers.Correlation metric for scRNA-seq data based on converting expression to Z-scores of differential expression. Robust to dropouts. Matlab implementation. [Data, 1847 human neuronal progenitor cells](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102934) <details>
    <summary>Paper</summary>
    Iacono, Giovanni, Elisabetta Mereu, Amy Guillaumet-Adkins, Roser Corominas, Ivon Cuscó, Gustavo Rodríguez-Esteban, Marta Gut, Luis Alberto Pérez-Jurado, Ivo Gut, and Holger Heyn. "BigSCale: An Analytical Framework for Big-Scale Single-Cell Data." https://doi.org/10.1101/gr.230771.117 Genome Research 28, no. 6 (June 2018): 878–90. 
</details>

- [CALISTA](https://github.com/CABSEL/CALISTA) - clustering, lineage reconstruction, transition gene identification, and cell pseudotime single cell transcriptional analysis. Analyses can be all or separate. Uses a likelihood-based approach based on probabilistic models of stochastic gene transcriptional bursts and random technical dropout events, so all analyses are compatible with each other. Input - a matrix of normalized, batch-removed log(RPKM) or log(TPM) or scaled UMIs. Methods detail statistical methodology. Matlab and R version <details>
    <summary>Paper</summary>
    Papili Gao N, Hartmann T, Fang T, Gunawan R. [CALISTA: Clustering and LINEAGE Inference in Single-Cell Transcriptional Analysis" https://doi.org/10.3389/fbioe.2020.00018 Frontiers in bioengineering and biotechnology. 2020 Feb 4;8:18.
</details>

- [demuxlet](https://github.com/statgen/demuxlet) - Introduces the ‘demuxlet’ algorithm, which enables genetic demultiplexing, doublet detection, and super-loading for droplet-based scRNA-seq. Recommended approach when samples have distinct genotypes <details>
    <summary>Paper</summary>
    Kang, Hyun Min, Meena Subramaniam, Sasha Targ, Michelle Nguyen, Lenka Maliskova, Elizabeth McCarthy, Eunice Wan, et al. "Multiplexed Droplet Single-Cell RNA-Sequencing Using Natural Genetic Variation." https://doi.org/10.1038/nbt.4042 Nature Biotechnology 36, no. 1 (January 2018): 89–94. 
</details>

- [dropEst](https://github.com/hms-dbmi/dropEst) - pipeline for pre-processing, mapping, QCing, filtering, and quantifying droplet-based scRNA-seq datasets. Input - FASTQ or BAM, output - an R-readable molecular count matrix. Written in C++ <details>
    <summary>Paper</summary>
    Petukhov, Viktor, Jimin Guo, Ninib Baryawno, Nicolas Severe, David T. Scadden, Maria G. Samsonova, and Peter V. Kharchenko. "DropEst: Pipeline for Accurate Estimation of Molecular Counts in Droplet-Based Single-Cell RNA-Seq Experiments." https://doi.org/10.1186/s13059-018-1449-6 Genome Biology 19, no. 1 (December 2018): 78. 
 </details>
 
 - [kallistobus](https://www.kallistobus.tools/) - fast pipeline for scRNA-seq processing. New BUS (Barcode, UMI, Set) format for storing and manipulating pseudoalignment results. Includes RNA velocity analysis. Python-based <details>
    <summary>Paper</summary>
    Melsted, Páll, A. Sina Booeshaghi, Fan Gao, Eduardo da Veiga Beltrame, Lambda Lu, Kristján Eldjárn Hjorleifsson, Jase Gehring, and Lior Pachter. "Modular and Efficient Pre-Processing of Single-Cell RNA-Seq." https://doi.org/10.1101/673285 Preprint. Bioinformatics, June 17, 2019. 
</details>

- [PyMINEr](https://www.sciencescott.com/pyminer) - Python-based scRNA-seq processing pipeline. Cell type identification, detection of cell type-enriched genes, pathway analysis, co-expression networks and graph theory approaches to interpreting gene expression. Notes on methods: modified K++ clustering, automatic detection of the number of cell types, co-expression and PPI networks. Input: .txt or .hdf5 files. Detailed analysis of several pancreatic datasets <details>
    <summary>Paper</summary>
    Tyler, Scott R., Pavana G. Rotti, Xingshen Sun, Yaling Yi, Weiliang Xie, Michael C. Winter, Miles J. Flamme-Wiese, et al. "PyMINEr Finds Gene and Autocrine-Paracrine Networks from Human Islet ScRNA-Seq." https://doi.org/10.1016/j.celrep.2019.01.063 Cell Reports 26, no. 7 (February 2019): 1951-1964.e8. 
</details>

- [Scanpy](https://github.com/theislab/scanpy) - Python-based pipeline for preprocessing, visualization, clustering, pseudotime and trajectory inference, differential expression and network simulation <details>
    <summary>Paper</summary>
    Wolf, F. Alexander, Philipp Angerer, and Fabian J. Theis. "SCANPY: Large-Scale Single-Cell Gene Expression Data Analysis." https://doi.org/10.1186/s13059-017-1382-0  Genome Biology 19, no. 1 (06 2018): 15. 
</details>

- [STQ](https://github.com/TheJacksonLaboratory/STQ) - PDX spatial transcriptomics Nextflow pipeline. 10x Genomics Visium with matched H&E images. Xenome read classification, spatial gene expression with Space Ranger, extraction of B-allele frequencies (can be used for CNV inference), splicing quantivication with Velocyto, image segmentation with Inception, StarDist, or HoVer-Net CNN. Can be applied to one-species samples. Optimized for SLURM HPC. <details>
    <summary>Paper</summary>
    Domanskyi, Sergii, Anuj Srivastava, Jessica Kaster, Haiyin Li, Meenhard Herlyn, Jill C. Rubinstein, and Jeffrey H. Chuang. “Nextflow Pipeline for Visium and H&E Data from Patient-Derived Xenograft Samples.” Preprint. Bioinformatics, July 30, 2023. https://doi.org/10.1101/2023.07.27.550727.
</details>

- [RAPIDS & Scanpy Single-Cell RNA-seq Workflow](https://github.com/clara-parabricks/rapids-single-cell-examples/blob/master/notebooks/hlca_lung_gpu_analysis.ipynb) - real-time analysis of scRNA-seq data on GPU. [Tweet](https://twitter.com/johnny_israeli/status/1265762506993135618?s=20)

- [SCRAT](https://github.com/zji90/SCRAT) - a Single-Cell Regulome Analysis Toolbox R package and a [Shiny web service](https://zhiji.shinyapps.io/scrat). scRNA-seq and scATAC-seq analyses. Input - BAM files. Summarizes regulatory activities on gene or transcription factor binding sites (by ENCODE clusters, motifs, DHSs, genes, MSigDb gene sets, or custom genomic features), clustering, cell annotation, differential gene/TF activity analysis. [Supplementary Table S1](https://academic.oup.com/bioinformatics/article/33/18/2930/3823309#supplementary-data) compares with other tools. Supplementary results demonstrate application to human and mouse ESC data. <details>
    <summary>Paper</summary>
    Ji, Zhicheng, Weiqiang Zhou, and Hongkai Ji. “Single-Cell Regulome Data Analysis by SCRAT.” Edited by Inanc Birol. Bioinformatics 33, no. 18 (September 15, 2017): 2930–32. https://doi.org/10.1093/bioinformatics/btx315.
</details>

- [scPipe](https://bioconductor.org/packages/release/bioc/html/scPipe.html) - A preprocessing pipeline for single cell RNA-seq data that starts from the fastq files and produces a gene count matrix with associated quality control information. It can process fastq data generated by CEL-seq, MARS-seq, Drop-seq, Chromium 10x and SMART-seq protocols. Modular, can swap tools like use different aligners. <details>
    <summary>Paper</summary>
    Tian et al. "[scPipe: A flexible R/Bioconductor preprocessing pipeline for single-cell RNA-sequencing data" https://doi.org/10.1371/journal.pcbi.1006361  PLOS Computational Biology, 2018. 
</details>


- [SEQC](https://github.com/ambrosejcarr/seqc) - Single-Cell Sequencing Quality Control and Processing Software, a general purpose method to build a count matrix from single cell sequencing reads, able to process data from inDrop, drop-seq, 10X, and Mars-Seq2 technologies. <details>
    <summary>Paper</summary>
    Azizi, Elham, Ambrose J. Carr, George Plitas, Andrew E. Cornish, Catherine Konopacki, Sandhya Prabhakaran, Juozas Nainys, et al. "Single-Cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment." https://doi.org/10.1016/j.cell.2018.05.060  Cell, June 2018. 
</details>

- [zUMIs](https://github.com/sdparekh/zUMIs) - scRNA-seq processing pipeline that handles barcodes and summarizes  UMIs using exonic or exonic + intronic mapped reads (improves clustering, DE detection). Adaptive downsampling of oversequenced libraries. STAR aligner, Rsubread::featureCounts counting UMIs in exons and introns. <details>
    <summary>Paper</summary>
    Parekh, Swati, Christoph Ziegenhain, Beate Vieth, Wolfgang Enard, and Ines Hellmann. "ZUMIs - A Fast and Flexible Pipeline to Process RNA Sequencing Data with UMIs." https://doi.org/10.1093/gigascience/giy059  GigaScience 7, no. 6 (01 2018). 
</details>

- [ramdaq](https://github.com/rikenbit/ramdaq) - pipeline to analyze data from full-length single-cell RNA sequencing (scRNA-seq) methods. [Documentation](https://github.com/rikenbit/ramdaq/tree/master/docs)

- STAR alignment parameters: `–outFilterType BySJout, –outFilterMultimapNmax 100, –limitOutSJcollapsed 2000000 –alignSJDBoverhangMin 8, –outFilterMismatchNoverLmax 0.04, –alignIntronMin 20, –alignIntronMax 1000000, –readFilesIn fastqrecords, –outSAMprimaryFlag AllBestScore, –outSAMtype BAM Unsorted`. From Azizi et al., “Single-Cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment.”

### Format conversion

- [sceasy](https://github.com/cellgeni/sceasy) - R package to convert different single-cell data formats to each other, supports Seurat, SingleCellExperiment, AnnData, Loom

- [scKirby](https://github.com/neurogenomics/scKirby) - R package for automated ingestion and conversion of various single-cell data formats (SingleCellExperiment, SummarizedExperiment, HDF5SummarizedExperiment, Seurat, H5Seurat, anndata, loom, loomR, CellDataSet/monocle, ExpressionSet, and more). 

- [zellkonverter](https://bioconductor.org/packages/zellkonverter/) - R package for conversion between scRNA-seq objects (the Bioconductor SingleCellExperiment data structure and the Python AnnData-based single-cell analysis environment). [Tweet](https://twitter.com/tangming2005/status/1466865990667542536?s=20)

### Visualization pipelines

- [ShinyCell](https://github.com/SGDDNB/ShinyCell) - R package to convert single-cell RNA-seq data into Shiny-based apps to share and visually explore the data. Input - h5ad, loom, SCE, Seurat object. Output - processed files and Shiny scripts. [Example 1](http://shinycell1.ddnetbio.com/), [Example 2](http://shinycell2.ddnetbio.com/). <details>
    <summary>Paper</summary>
    Ouyang, John F., Uma S. Kamaraj, Elaine Y. Cao, and Owen JL Rackham. "ShinyCell: simple and sharable visualization of single-cell gene expression data." Bioinformatics 37, no. 19 (1 October 2021): https://doi.org/10.1093/bioinformatics/btab209
</details>

- [Kana](https://www.jkanche.com/kana/) - single-cell analysis in the browser, by Jayaram Kancherla ([@jkanche](https://github.com/jkanche)), Aaron Lun ([@LTLA](https://github.com/LTLA)). Client-side computations using the WebAssembly framework. Input - 10X genomics CellRanger's output (Matrix Market format), csv matrix or .h5 files. Preprocessing (removal of low-quality cells, Normalization and log-transformation, Modelling of the mean-variance trend across genes), PCA, Clustering (t-SNE/UMAP), Marker detection, custom cluster definition and marker analysis. Works with scATAC-seq data. [GitHub](https://github.com/jkanche/kana), [Tweet](https://twitter.com/jayaram/status/1480599647039016962?s=20). <details>
    <summary>Paper</summary>
    Lun, Aaron, and Jayaram Kancherla. “Powering Single-Cell Analyses in the Browser with WebAssembly.” Preprint. Bioinformatics, March 4, 2022. https://doi.org/10.1101/2022.03.02.482701.
</details>

- [cellxgene](https://github.com/chanzuckerberg/cellxgene) - An interactive exploratory visualization tool for single-cell transcriptomics data, web and desktop versions. Input - matrix-form datasets, metadata, pre-computed embeddings/clustering. Compatible with Seurat, Scanpy, Bioconductor, scVI [GitHub](https://github.com/chanzuckerberg/cellxgene) <details>
    <summary>Paper</summary>
    Megill, Colin, Bruce Martin, Charlotte Weaver, Sidney Bell, Lia Prins, Seve Badajoz, Brian McCandless, et al. "Cellxgene: A Performant, Scalable Exploration Platform for High Dimensional Sparse Matrices" https://doi.org/10.1101/2021.04.05.438318  Preprint. Systems Biology, April 6, 2021.
</details>

- [iCellR](https://github.com/rezakj/iCellR) - Single (i) Cell R package (iCellR) is an interactive R package to work with high-throughput single cell sequencing technologies (i.e scRNA-seq, scVDJ-seq and CITE-seq). <details>
    <summary>Paper</summary>
    Khodadadi-Jamayran, Alireza, Joseph Pucella, Hua Zhou, Nicole Doudican, John Carucci, Adriana Heguy, Boris Reizis, and Aristotelis Tsirigos. "ICellR: Combined Coverage Correction and Principal Component Alignment for Batch Alignment in Single-Cell Sequencing Analysis" https://www.biorxiv.org/content/10.1101/2020.03.31.019109v1.full  BioRxiv, April 1, 2020
</details>

- [Cerebro](https://github.com/romanhaa/Cerebro) - interactive scRNA-seq visualization from a Seurat object (v2 or 3), dimensionality reduction, clustering, identification and visualization of marker genes, enriched pathways (EnrichR), signatures (MSigDb),  expression of individual genes. [cerebroPrepare](https://github.com/romanhaa/cerebroPrepare) R package saves the Seurat object, to be visualized with [cerebroApp Shiny app](https://github.com/romanhaa/cerebroApp). Standalone and Docker versions are available. [GitHub](https://github.com/romanhaa/Cerebro). <details>
    <summary>Paper</summary>
    Hillje, Roman, Pier Giuseppe Pelicci, and Lucilla Luzi. "Cerebro: Interactive Visualization of ScRNA-Seq Data" https://doi.org/10.1093/bioinformatics/btz877  Bioinformatics, 1 April 2020 
</details>

- [iS-CellR](https://github.com/immcore/iS-CellR) - a Shiny app for scRNA-seq analysis. Can be insalled locally, run from GitHub, Docker. Input - count matrix. Filtering, normalization, dimensionality reduction, clustering, differential expression, co-expression, reports. <details>
    <summary>Paper</summary>
    Patel, Mitulkumar V. "IS-CellR: A User-Friendly Tool for Analyzing and Visualizing Single-Cell RNA Sequencing Data" https://doi.org/10.1093/bioinformatics/bty517   Bioinformatics 34, no. 24 (December 15, 2018)
</details>

- [iSEE](https://github.com/kevinrue/iSEEWorkshop2019) - Shiny app for interactive visualization of SummarizedExperiment scRNA-seq objects. [GitHub](https://github.com/csoneson/iSEE), [RNA-seq blog post](https://www.rna-seqblog.com/isee-an-interactive-shiny-based-graphical-user-interface-for-exploring-data-stored-in-summarizedexperiment-objects/), [Workshop](https://github.com/kevinrue/iSEEWorkshop2019). <details>
    <summary>Paper</summary>
    -Rue-Albrecht, Kevin, Federico Marini, Charlotte Soneson, and Aaron T.L. Lun. "ISEE: Interactive SummarizedExperiment Explorer" https://doi.org/10.12688/f1000research.14966.1  F1000Research 7 (June 14, 2018)
</details>

- [SPRING](https://github.com/AllonKleinLab/SPRING_dev) - a pipeline for data filtering, normalization and visualization using force-directed layout of k-nearest-neighbor graph. [Web-based](https://kleintools.hms.harvard.edu/tools/spring.html) (10,000 cells max) and [GitHub](https://github.com/AllonKleinLab/SPRING_dev). <details>
    <summary>Paper</summary>
    Weinreb, Caleb, Samuel Wolock, and Allon M. Klein. "SPRING: A Kinetic Interface for Visualizing High Dimensional Single-Cell Expression Data" https://doi.org/10.1093/bioinformatics/btx792  Bioinformatics (Oxford, England) 34, no. 7 (April 1, 2018)
</details>

- [Granatum](http://garmiregroup.org/granatum/app) - web-based scRNA-seq analysis. list of modules, including plate merging and batch-effect removal, outlier-sample removal, gene-expression normalization, imputation, gene filtering, cell clustering, differential gene expression analysis, pathway/ontology enrichment analysis, protein network interaction visualization, and pseudo-time cell series reconstruction. [Twitter](https://twitter.com/GarmireGroup/status/1185269818015940609). <details>
    <summary>Paper</summary>
    Zhu, Xun, Thomas K. Wolfgruber, Austin Tasato, Cédric Arisdakessian, David G. Garmire, and Lana X. Garmire. "Granatum: A Graphical Single-Cell RNA-Seq Analysis Pipeline for Genomics Scientists" https://doi.org/10.1186/s13073-017-0492-3  Genome Medicine 9, no. 1 (December 2017). 
</details>

- [SCope](https://github.com/aertslab/SCope) - Fast visualization tool for large-scale and high dimensional single-cell data in `.loom` format. R and Python scripts for converting scRNA-seq data to `.loom` format. 

- [singleCellTK](https://bioconductor.org/packages/singleCellTK/) - R/Shiny package for an interactive scRNA-Seq analysis. Input, raw counts in SingleCellExperiment. Analysis: filtering raw results, clustering, batch correction, differential expression, pathway enrichment, and scRNA-Seq study design.

- [scDataviz](https://github.com/kevinblighe/scDataviz) - single cell data vizualization and downstream analyses, by Kevin Blighe

- [scOrange](https://singlecell.biolab.si/) - visual pipeline builder for an in-depth analysis and visualization of scRNA-seq data. Works with 10X data, tab-delimited. Filtering, preprocessiong, differential gene expression, marker analysis, enrichment analysis, batch removal, clustering, tSNE. [Screenshots](https://singlecell.biolab.si/screenshots/), [Short video tutorials](https://www.youtube.com/playlist?list=PLmNPvQr9Tf-a4MrEG5thq3qzlkrF5NFbC). Python-based, Conda-installable. [GitHub](https://github.com/biolab/orange3-single-cell)

- [scCustomize](https://github.com/samuel-marsh/scCustomize) - an R package, Collection of functions created and/or curated to aid in the visualization and analysis of single-cell data. Extends Seurat, Liger visualization, helper functions to enhance analysis of Seurat objects.

- [UCSC Single Cell Browser](https://github.com/maximilianh/cellBrowser) - Python pipeline and Javascript scatter plot library for single-cell datasets. Pre-process an expression matrix by filtering, PCA, nearest-neighbors, clustering, t-SNE and UMAP and formats them for cbBuild. [Demo that includes several landmark datasets](https://cells.ucsc.edu/)

## Quality control

- [miQC](https://bioconductor.org/packages/miQC/) - data-driven identification of cells with high mitochondrial content (likely, dead cells) from scRNA-seq data. Joint statistical model the proportion of reads mapping to mtDNA genes and the number of detected genes, EM for parameter estimation (flexmix). Tested on various datasets processed with CellRanged and salon alevin - improves removal of compromised cells as compared with hard threshold. Bioconductor R package, integrates with scater. <details>
    <summary>Paper</summary>
    Hippen, Ariel A., Matias M. Falco, Lukas M. Weber, Erdogan Pekcan Erkan, Kaiyang Zhang, Jennifer Anne Doherty, Anna Vähärautio, Casey S. Greene, and Stephanie C. Hicks. "MiQC: An Adaptive Probabilistic Framework for Quality Control of Single-Cell RNA-Sequencing Data" https://doi.org/10.1371/journal.pcbi.1009290  PLOS Computational Biology, (August 24, 2021)
</details>

- [doubletD](https://github.com/elkebir-group/doubletD) - doublet detection in single-cell DNA-seq data. doublets in scRNA-seq data have a characteristic variant allele frequency spectrum due to increased copy number and allelic dropout. A maximum likelihood approach with a closed-form solution - stats in Methods. Simulated and real data, outperforms SCG, Scrublet, robust to the presence of CNAs, mixture of two cell types. Python3 implementation. <details>
    <summary>Paper</summary>
    Weber, Leah L, Palash Sashittal, and Mohammed El-Kebir. "DoubletD: Detecting Doublets in Single-Cell DNA Sequencing Data" https://doi.org/10.1093/bioinformatics/btab266  Bioinformatics, (August 4, 2021)
</details>

- [DropletQC](https://github.com/powellgenomicslab/DropletQC) - empty droplet identification. A novel metric - nuclear fraction. Damaged cells due to the depletion of cytoplasmic RNA will have a higher nuclear fraction compared to intact cells. Compared with 10X Cell Ranger, CellBlender, EmptyNN, EmptyDrops. [Scripts](https://github.com/powellgenomicslab/dropletQC_paper). <details>
    <summary>Paper</summary>
    Muskovic, Walter. “DropletQC: Improved Identification of Empty Droplets and Damaged Cells in Single-Cell RNA-Seq Data,” 2021, 9.
</details>

- [DropletUtils](https://bioconductor.org/packages/DropletUtils/) - Provides a number of utility functions for handling single-cell (RNA-seq) data from droplet technologies such as 10X Genomics. This includes data loading, identification of cells from empty droplets, removal of barcode-swapped pseudo-cells, and downsampling of the count matrix. <details>
    <summary>Paper</summary>
    Lun ATL, Riesenfeld S, Andrews T, Dao T, Gomes T, participants in the 1st Human Cell Atlas Jamboree, Marioni JC (2019). "EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data" https://doi.org/10.1186/s13059-019-1662-y  Genome Biol.
</details>

- [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) - doublet detection using gene expression data. Simulates artificial doublets, incorporate them into existing scRNA-seq data. Integrates with Seurat (Figure 1). Three input parameters (the expected number of doublets, the number of artificial doublets pN, the neighborhood size pN), need to be tailored to data with different number of cell types and magnitudes of transcriptional heterogeneity. Bimodality Coefficient maximization to select pN. Benchmarked against ground-truth scRNA-seq datasets. Not optimal for homogeneous data. <details>
    <summary>Paper</summary>
    -McGinnis, Christopher S., Lyndsay M. Murrow, and Zev J. Gartner. "DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors" https://doi.org/10.1016/j.cels.2019.03.003  Cell Systems 8, no. 4 (April 2019)
</details>

- [scrublet](https://github.com/AllonKleinLab/scrublet) - Detect doublets in single-cell RNA-seq data. <details>
    <summary>Paper</summary>
    Wolock, Samuel L, Romain Lopez, and Allon M Klein. "Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data" https://doi.org/10.1101/357368  Preprint. Bioinformatics, July 9, 2018. 
</details>

- [scater](https://bioconductor.org/packages/scater/) - A collection of tools for doing various analyses of single-cell RNA-seq gene expression data, with a focus on quality control. <details>
    <summary>Paper</summary>
    McCarthy et al. "[Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R" https://doi.org/10.1093/bioinformatics/btw777   Bioinformatics, 2017. 
</details>

- [celloline](https://github.com/Teichlab/celloline) - A pipeline to remove low quality single cell files. Figure 2 - 20 biological and technical features used for filtering. High mitochondrial genes = broken cells. <details>
    <summary>Paper</summary>
    Ilicic, Tomislav, Jong Kyoung Kim, Aleksandra A. Kolodziejczyk, Frederik Otzen Bagger, Davis James McCarthy, John C. Marioni, and Sarah A. Teichmann. "Classification of Low Quality Cells from Single-Cell RNA-Seq Data" https://doi.org/10.1186/s13059-016-0888-1  Genome Biology 17 (February 17, 2016)
</details>

## Normalization

- [sctransform](https://cran.r-project.org/web/packages/sctransform/) - using the Pearson residuals from regularized negative binomial regression where sequencing depth is utilized as a covariate to remove technical artifacts. Interfaces with Seurat. <details>
    <summary>Paper</summary>
    Hafemeister, Christoph, and Rahul Satija. "Normalization and Variance Stabilization of Single-Cell RNA-Seq Data Using Regularized Negative Binomial Regression" https://doi.org/10.1101/576827  BioRxiv, March 14, 2019.
</details>

- [SCnorm](https://www.biostat.wisc.edu/~kendzior/SCNORM/) - normalization for single-cell data. Quantile regression to estimate the dependence of transcript expression on sequencing depth for every gene. Genes with similar dependence are then grouped, and a second quantile regression is used to estimate scale factors within each group. Within-group adjustment for sequencing depth is then performed using the estimated scale factors to provide normalized estimates of expression. Good statistical methods description. <details>
    <summary>Paper</summary>
    Bacher, Rhonda, Li-Fang Chu, Ning Leng, Audrey P Gasch, James A Thomson, Ron M Stewart, Michael Newton, and Christina Kendziorski. "SCnorm: Robust Normalization of Single-Cell RNA-Seq Data" https://doi.org/10.1038/nmeth.4263  Nature Methods 14, no. 6 (April 17, 2017)
</details>

## Integration, Batch correction

- [Evaluation of 10 single-cell data integration methods and 4 preprocessing combinations](https://github.com/theislab/scib) on 77 batches of gene expression, chromatin accessibility, and simulated data (Table 1) in 9 integration tasks using 14 evaluation metrics. BBKNN, Scanorama, scVI perform well on complex tasks, Seurat performs well on simpler tasks but may eliminate biological signal. the use of Seurat v3 and Harmony is appropriate for simple integration tasks with distinct batch and biological structure. Batch in ATAC-seq is the most difficult to remove. Jupyter notebooks for full reproducibility. <details>
    <summary>Paper</summary>
    Luecken, Md, M Büttner, K Chaichoompu, A Danese, M Interlandi, Mf Mueller, Dc Strobl, et al. “Benchmarking Atlas-Level Data Integration in Single-Cell Genomics.” Nature Methods, 23 December 2021 https://doi.org/10.1038/s41592-021-01336-8
</details>

- [Benchmark of 14 methods for scRNA-seq batch correction](https://github.com/JinmiaoChenLab/Batch-effect-removal-benchmarking). Using five scenarios: different technologies and same cells, non-identical cell types, multiple batches, big data, simulated data. Four benchmarking metrics. Harmony (fast), LIGER, and Seurat 3 perform well overall. For differential expression, ComBat, limma, MNN Correct perform well. Detailed description of 9 datasets and download links. [Data and scripts](https://github.com/JinmiaoChenLab/Batch-effect-removal-benchmarking). <details>
    <summary>Paper</summary>
    Tran, Hoa Thi Nhu, Kok Siong Ang, Marion Chevrier, Xiaomeng Zhang, Nicole Yee Shin Lee, Michelle Goh, and Jinmiao Chen. "A Benchmark of Batch-Effect Correction Methods for Single-Cell RNA Sequencing Data" https://doi.org/10.1186/s13059-019-1850-9  Genome Biology 21, no. 1 (December 2020)
</details>

- [CellANOVA](https://github.com/Janezjz/cellanova) (cell state space analysis of variance) - scRNA-seq data integration. Statistical model to separate unwanted and biological variation. Operates on top of existing integration methods (Harmony, Seurat). Requires a control-pool set of samples, a set of samples whereby variation beyond what is preserved by the existing integration are not of interest to the study. Used to estimate a latent linear space that captures cell- and gene-specific unwanted batch variations. Tested on 4 experimental designs (case-control, longitudinal, irregular block design, scRNA and snRNA integration). Outperforms Seurat, Harmony, LIGER, Symphony integration, as judged by local inverse Simpson’s index (LISI). Python implementation. <details>
    <summary>Paper</summary>
    Zhang, Zhaojun, Divij Mathew, Tristan Lim, Kaishu Mason, Clara Morral Martinez, Sijia Huang, E. John Wherry, et al. “Signal Recovery in Single Cell Batch Integration,” May 8, 2023. https://doi.org/10.1101/2023.05.05.539614.
</details>

- [SMAI](https://github.com/rongstat/SMAI) - spectral manifold alignment and inference framework for alignment single-cell sequencing data. Includes SMAI-test of (partial) alignability against the null hypothesis that two single-cell datasets are alignable up to some similarity transformation, that is, combinations of scaling, translation, and rotation. SMAI-align incorporates a high-dimensional shuffled Procrustes analysis, which iteratively searches for the sample correspondence and the best similarity transformation that minimizes the discrepancy between the intrinsic low-dimensional signal structures of the datasets. References that current integration methods distort the biology. Compared with Seurat. LIGER, Harmony, fastMNN, Scanorama, evaluated on integration of diverse tissues, technologies, 13 integration tasks. Assessment of false positives in differential expression, outperforms all methods. R implementation, works on Seurat objects, [tutorial](https://rongstat.github.io/SMAI_guide.io/SMAI-tutorial.html). <details>
    <summary>Paper</summary>
    Ma, Rong, Eric D Sun, David Donoho, and James Zou. “Principled and Interpretable Alignability Testing and Integration of Single-Cell Data” 121, no. 10 (2024). https://doi.org/10.1073/pnas.2313719121
</details>

- [MOJITOO](https://github.com/CostaLab/MOJITOO) (Multi-mOdal Joint IntegraTion of cOmpOnents) - a single cell multi-modal integration method (does not require shared features, e.g., genes), uses canonical correlation analysis. Introduction of two main frameworks: metric learning (WNN, Schema) and latent variable learning (MOFA, scAI, totalVI, LIGER). Benchmarked against them on bi- and trimodal data (RNA, protein, ATAC). R implementation, compatible with Seurat, Signac. <details>
    <summary>Paper</summary>
    Cheng, Mingbo, Zhijian Li, and Ivan Gesteira Costa Filho. “MOJITOO: A Fast and Universal Method for Integration of Multimodal Single Cell Data.” Preprint. Bioinformatics, January 21, 2022. https://doi.org/10.1101/2022.01.19.476907.
</details>

- [RPCI](https://github.com/bioinfoDZ/RISC) (Reference Principal Component Integration) - R package RISC for integration of scRNA-seq data using the gene eigenvectors from a reference dataset as a single reference space. Compared with CCA, shared cell type-based strategies. Tested on simulated and experimental datasets against 11 other integration approaches (Scanorama, Harmony, fastMNN, Anchor, among others) using four metrics (kBET scores, LISI, ARI, SW). Robust when using two and more datasets (e.g., timecourse). Scanorama generally ranks second. [Code to reproduce the paper](https://codeocean.com/capsule/9098032/tree/v1). <details>
    <summary>Paper</summary>
    Liu, Yang, Tao Wang, Bin Zhou, and Deyou Zheng. “Robust Integration of Multiple Single-Cell RNA Sequencing Datasets Using a Single Reference Space.” Nature Biotechnology, March 25, 2021. https://doi.org/10.1038/s41587-021-00859-x.
</details>

- [MOFA2](https://github.com/bioFAM/MOFA2) - Multi-Omics Factor Analysis v2 (MOFA+), a statistical framework for the integration of single-cell multi-modal data. Reconstructs a low-dimensional representation of the data using variational inference (a stochastic variant parallelizable on GPU, 20-fold speed increase). Supports sparsity constraints, allowing to jointly model variation across multiple sample groups and data modalities. Infers K latent factors with associated feature weight matrices (per data modality, Figure 1a) that can be used for clustering, trajectory inference, variance decomposition etc. Input - multiple datasets measuring non-overlapping modalities, cells grouped by experiments, batches, or conditions. Python and R implementation. <details>
    <summary>Paper</summary>
    Argelaguet, Ricard, Damien Arnol, Danila Bredikhin, Yonatan Deloro, Britta Velten, John C. Marioni, and Oliver Stegle. “MOFA+: A Statistical Framework for Comprehensive Integration of Multi-Modal Single-Cell Data.” Genome Biology 21, no. 1 (December 2020): 111. https://doi.org/10.1186/s13059-020-02015-1.
</details>

- [CarDEC](https://github.com/jlakkis/CarDEC) - Count adapted regularized Deep Embedded Clustering, a joint deep learning model that simultaneously clusters, denoises, corrects for multiple batch effects in gene expression space (Figure 1). Outperforms scVI, DCA, MNN, scDeepCluster. Separately treats highly and lowly variable genes. Improves integration of omics generated by multiple technologies, pseudotime reconstruction. <details>
    <summary>Paper</summary>
    Lakkis, Justin, David Wang, Yuanchao Zhang, Gang Hu, Kui Wang, Huize Pan, Lyle Ungar, Muredach P. Reilly, Xiangjie Li, and Mingyao Li. “A Joint Deep Learning Model for Simultaneous Batch Effect Correction, Denoising and Clustering in Single-Cell Transcriptomics.” Preprint. Bioinformatics, September 25, 2020. https://doi.org/10.1101/2020.09.23.310003
</details>

- [iCellR](https://cran.r-project.org/web/packages/iCellR/index.html) - batch correction in scRNA-seq data. Combined Coverage Correction Alignment (CCCA) and Combined Principal Component Alignment (CPCA). CCCA - PCA into 30 dimensions, for each cell, take k=10 nearest neighbors, average gene expression, thus imputing the adjusted matrix. CPCA skips imputation, instead PCs themselves get averaged. Similar performance. Tested on nine PBMC datasets provided by the Broad institute to test batch effect. Outperforms MAGIC. [Data in text and .rda formats](https://genome.med.nyu.edu/results/external/iCellR/data/). <details>
    <summary>Paper</summary>
    Khodadadi-Jamayran, Alireza, Joseph Pucella, Hua Zhou, Nicole Doudican, John Carucci, Adriana Heguy, Boris Reizis, and Aristotelis Tsirigos. "ICellR: Combined Coverage Correction and Principal Component Alignment for Batch Alignment in Single-Cell Sequencing Analysis" https://doi.org/10.1101/2020.03.31.019109  Preprint. Bioinformatics, April 1, 2020
</details>

- [scAlign](https://bioconductor.org/packages/scAlign/) - a deep learning method for alignment and integration of scRNA-seq datasets. Bidirectional mapping via a low-dimensional space. Can perform unsupervised, semi-supervised, and supervised (by cell type labels) integration. Outperforms scVI, MNN, scmap, MINT, scMERGE, Scanorama, Seurat (two latter perform well). [GitHub](https://github.com/quon-titative-biology/scAlign), [Bioconductor R package](https://bioconductor.org/packages/scAlign/). <details>
    <summary>Paper</summary>
    Johansen, Nelson, and Gerald Quon. “ScAlign: A Tool for Alignment, Integration, and Rare Cell Identification from ScRNA-Seq Data.” Genome Biology 20, no. 1 (December 2019): 166. https://doi.org/10.1186/s13059-019-1766-4
</details>

- [BERMUDA](https://github.com/txWang/BERMUDA) - Batch Effect ReMoval Using Deep Autoencoders, for scRNA-seq data. Requires batches to share at least one common cell type. Five step framework: 1) preprocessing, 2) clustering of cells in each batch individually, 3) identifying similar cell clusters across different batches, 4) removing batch effect by training an autoencoder, 5) further analysis of batch-corrected data. Tested on simulated (splatter) and experimental (10X genomics) data. <details>
    <summary>Paper</summary>
    Wang, Tongxin, Travis S. Johnson, Wei Shao, Zixiao Lu, Bryan R. Helm, Jie Zhang, and Kun Huang. "BERMUDA: A Novel Deep Transfer Learning Method for Single-Cell RNA Sequencing Batch Correction Reveals Hidden High-Resolution Cellular Subtypes" https://doi.org/10.1186/s13059-019-1764-6  Genome Biology 20, no. 1 (December 2019). 
</details>

- [Scanorama](http://cb.csail.mit.edu/cb/scanorama/) - Python tool, integrates scRNA-seq datasets, identifies the shared cell types among all pairs of datasets (mutual nearest-neighbors matching in low-dimensional (100 SVD components) space) and uses this info for batch correction and merging. Tested on 26 scRNA-seq datasets, 9 technologies, simulated data. Compared with Seurat's CCA, scran MNN. [Links to many public datasets](https://www.nature.com/articles/s41587-019-0113-3#data-availability). <details>
    <summary>Paper</summary>
    Hie, Brian, Bryan Bryson, and Bonnie Berger. “Efficient Integration of Heterogeneous Single-Cell Transcriptomes Using Scanorama.” Nature Biotechnology, May 6, 2019. https://doi.org/10.1038/s41587-019-0113-3.
</details>

- [BBKNN](https://github.com/Teichlab/bbknn) (batch balanced k nearest neighbours) - batch correction for scRNA-seq data. Neighborhood graphs, balanced across all batches of the data, separately for each batch, that are merged. Main assumption (as in mnnCorrect) - at least some cells of the same type exist across batches. Preserves data structure allowing subsequent embedding, trajectory reconstruction. Python, compatible with SCANPY, very fast. <details>
    <summary>Paper</summary>
    Polański, Krzysztof, Matthew D Young, Zhichao Miao, Kerstin B Meyer, Sarah A Teichmann, and Jong-Eun Park. "BBKNN: Fast Batch Alignment of Single Cell Transcriptomes" https://doi.org/10.1093/bioinformatics/btz625  Bioinformatics, August 10, 2019
</details>

- [conos](https://github.com/hms-dbmi/conos) - joint analysis of scRNA-seq datasets through inter-sample mapping (mutual nearest-neighbor mapping) and constructing a joint graph. [Analysis scripts](http://pklab.med.harvard.edu/peterk/conos/). <details>
    <summary>Paper</summary>
    Barkas, Nikolas, Viktor Petukhov, Daria Nikolaeva, Yaroslav Lozinsky, Samuel Demharter, Konstantin Khodosevich, and Peter V. Kharchenko. "Joint Analysis of Heterogeneous Single-Cell RNA-Seq Dataset Collections" https://doi.org/10.1038/s41592-019-0466-z  Nature Methods, July 15, 2019. 
</details>

- [LIGER](https://github.com/MacoskoLab/liger) - R package for integrating and analyzing multiple single-cell datasets, across conditions, technologies (scRNA-seq and methylation), or species (human and mouse). Integrative nonnegative matrix factorization (W and H matrices), dataset-specific and shared patterns (metagenes, matrix H). Graphs of factor loadings onto these patterns (shared factor neighborhood graph), then comparing patterns. Alignment and agreement metrics to assess performance, LIGER outperforms Seurat on agreement. Analysis of published blood cells, brain. [Human/mouse brain data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126836). <details>
    <summary>Paper</summary>
    Welch, Joshua D., Velina Kozareva, Ashley Ferreira, Charles Vanderburg, Carly Martin, and Evan Z. Macosko. "Single-Cell Multi-Omic Integration Compares and Contrasts Features of Brain Cell Identity" https://doi.org/10.1016/j.cell.2019.05.006  Cell 177, no. 7 (June 13, 2019)
</details>

- [scMerge](https://github.com/SydneyBioX/scMerge/) - R package for batch effect removal and normalizing of multipe scRNA-seq datasets. fastRUVIII batch removal method. Tested on 14 datasets, compared with scran, MNN, ComBat, Seurat, ZINB-WaVE using Silhouette, ARI - better separation of clusters, pseudotime reconstruction. <details>
    <summary>Paper</summary>
    Lin, Yingxin, Shila Ghazanfar, Kevin Wang, Johann A. Gagnon-Bartsch, Kitty K. Lo, Xianbin Su, Ze-Guang Han, et al. "ScMerge: Integration of Multiple Single-Cell Transcriptomics Datasets Leveraging Stable Expression and Pseudo-Replication" https://doi.org/10.1101/393280 September 12, 2018. 
</details>

- [cellHarmony](https://github.com/AltAnalyze/cellHarmony-Align) - a Python 2.7 package for integration and comparison of scRNA-seq datasets, a part of [AltAnalyze](http://www.altanalyze.org/) workflow for RNA-Seq gene, splicing and pathway analysis. Uses a community clustering to produce a network graph and define communities in both the reference and query datasets, and alignment (label projection) strategy. Table 1 - comparison with other joint alignment and label projection methods. Differential expression using empirical Bayes moderated t-test and FDR, global, local, and co-regulated comparisons. Tested on several datasets, improves similarity to the author-defined ground truth. Support for 10x Genomics data format. <details>
    <summary>Paper</summary>
    DePasquale, Erica AK, Phillip Dexheimer, Daniel Schnell, Kyle Ferchen, Stuart Hay, Íñigo Valiente-Alandí, Burns C. Blaxall, H. Leighton Grimes, and Nathan Salomonis. “CellHarmony: Cell-Level Matching and Holistic Comparison of Single-Cell Transcriptomes.” Preprint. Bioinformatics, September 8, 2018. https://doi.org/10.1101/412080.
</details>

- [MNN](https://bioconductor.org/packages/scran/) - mutual nearest neighbors method for single-cell batch correction. Assumptions: MNN exist between batches, batch is orthogonal to the biology. Cosine normalization, Euclidean distance, a pair-specific barch-correction vector as a vector difference between the expression profiles of the paired cells using selected genes of interest and hypervariable genes. Supplementary note 5 - algorithm. mnnCorrect function in the [scran](https://bioconductor.org/packages/scran/) package. [Code for paper](https://github.com/MarioniLab/MNN2017/). <details>
    <summary>Paper</summary>
    Haghverdi, Laleh, Aaron T L Lun, Michael D Morgan, and John C Marioni. "Batch Effects in Single-Cell RNA-Sequencing Data Are Corrected by Matching Mutual Nearest Neighbors" https://doi.org/10.1038/nbt.4091  Nature Biotechnology, April 2, 2018. 
</details>

- [batchelor](https://bioconductor.org/packages/batchelor/) - Single-Cell Batch Correction Methods, by Aaron Lun. <details>
    <summary>Paper</summary>   
    Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018). "Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors" https://doi.org10.1038/nbt.4091  Nat. Biotechnol.
</details>

- [scLVM](https://github.com/PMBio/scLVM) - a modelling framework for single-cell RNA-seq data that can be used to dissect the observed heterogeneity into different sources and remove the variation explained by latent variables. Can correct for the cell cycle effect. Applied to naive T cells differentiating into TH2 cells. <details>
    <summary>Paper</summary>
    Buettner, Florian, Kedar N Natarajan, F Paolo Casale, Valentina Proserpio, Antonio Scialdone, Fabian J Theis, Sarah A Teichmann, John C Marioni, and Oliver Stegle. "Computational Analysis of Cell-to-Cell Heterogeneity in Single-Cell RNA-Sequencing Data Reveals Hidden Subpopulations of Cells" https://doi.org/10.1038/nbt.3102  Nature Biotechnology 33, no. 2 (March 2015)
    
    Buettner, Florian, Naruemon Pratanwanich, Davis J. McCarthy, John C. Marioni, and Oliver Stegle. "F-ScLVM: Scalable and Versatile Factor Analysis for Single-Cell RNA-Seq" https://doi.org/10.1186/s13059-017-1334-8  Genome Biology 18, no. 1 (December 2017) - f-scLVM - factorial single-cell latent variable model guided by pathway annotations to infer interpretable factors behind heterogeneity. PCA components are annotated by correlated genes and their enrichment in pathways. Docomposition of the original gene expression matrix to a sum of annotated, unannotated, and confounding components. Applied to their own naive T to TH2 cells, mESCs, reanalyzed 3005 neuronal cells. Simulated data. https://github.com/bioFAM/slalom
</details>

## Imputation

[Assessment of 18 scRNA-seq imputation methods](https://github.com/Winnie09/imputationBenchmark) (model-based, smooth-based, deep learning, matrix decomposition). Similarity of scRNA- and bulk RNA-seq profiles (Spearman), differential expression (MAST and Wilcoxon), clustering (k-means, Louvain), trajectory reconstruction (Monocle 2, TSCAN), didn't test velocity. scran for normalization. Imputation methods improve correlation with bulk RNA-seq, but have minimal effect on downstream analyses. MAGIC, kNN-smoothing, SAVER perform well overall. [Plate- and droplet-derived scRNA-seq cell line data, Additional File 4](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-020-02132-x/MediaObjects/13059_2020_2132_MOESM4_ESM.xlsx)), [Summary table of the functionality of all imputation methods, Additional File 5](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-020-02132-x/MediaObjects/13059_2020_2132_MOESM4_ESM.xlsx). <details>
    <summary>Paper</summary>
    Hou, Wenpin, Zhicheng Ji, Hongkai Ji, and Stephanie C. Hicks. "A Systematic Evaluation of Single-Cell RNA-Sequencing Imputation Methods" https://doi.org/10.1186/s13059-020-02132-x  Genome Biology 21, no. 1 (December 2020)
</details>

- [Deepimpute](https://github.com/lanagarmire/DeepImpute) - scRNA-seq imputation using deep neural networks. Sub-networks, each processes up to 512 genes needed to be imputed. Four layers: Input - dense (ReLU activation) - 20% dropout - output. MSE as loss function. Outperforms MAGIC, DrImpute, ScImpute, SAVER, VIPER, and DCA on multiple metrics (PCC, several clustering metrics). Using 9 datasets. <details>
    <summary>Paper</summary>
    Arisdakessian, Cédric, Olivier Poirion, Breck Yunits, Xun Zhu, and Lana X. Garmire. "DeepImpute: An Accurate, Fast, and Scalable Deep Neural Network Method to Impute Single-Cell RNA-Seq Data" https://doi.org/10.1186/s13059-019-1837-6  Genome Biology 20, no. 1 (December 2019)
</details>

- [SCRABBLE](https://github.com/tanlabcode/SCRABBLE) - scRNA-seq imputation constraining on bulk RNA-seq data. Matrix regularzation optimizing a three-term objective function. Compared with DrImpute, scImpute, MAGIC, VIPER on simulated and real data. [Datasets](https://github.com/tanlabcode/SCRABBLE_PAPER). R and Matlab implementation. <details>
    <summary>Paper</summary>
    Peng, Tao, Qin Zhu, Penghang Yin, and Kai Tan. "SCRABBLE: Single-Cell RNA-Seq Imputation Constrained by Bulk RNA-Seq Data" https://doi.org/10.1186/s13059-019-1681-8  Genome Biology 20, no. 1 (December 2019)
</details>

- [ENHANCE](https://github.com/yanailab/enhance-R), an algorithm that denoises single-cell RNA-Seq data by first performing nearest-neighbor aggregation and then inferring expression levels from principal components. Variance-stabilizing normalization of the data before PCA. Implements its own simulation procedure for simulating sampling noise. Outperforms MAGIC, SAVER, ALRA. [Python](https://github.com/yanailab/enhance), and [R](https://github.com/yanailab/enhance-R) implementations. <details>
    <summary>Paper</summary>
    Wagner, Florian, Dalia Barkley, and Itai Yanai. "ENHANCE: Accurate Denoising of Single-Cell RNA-Seq Data" https://doi.org/10.1101/655365  Preprint. Bioinformatics, June 3, 2019. 
</details>

- [scHinter](https://github.com/BMILAB/scHinter) - imputation for small-size scRNA-seq datasets. Three modules: voting-based ensemble distance for learning cell-cell similarity, a SMOTE-based random interpolation module for imputing dropout events, and a hierarchical model for multi-layer random interpolation. [RNA-seq blog](https://www.rna-seqblog.com/schinter-imputing-dropout-events-for-single-cell-rna-seq-data-with-limited-sample-size/). <details>
    <summary>Paper</summary>
    Ye, Pengchao, Wenbin Ye, Congting Ye, Shuchao Li, Lishan Ye, Guoli Ji, and Xiaohui Wu. "ScHinter: Imputing Dropout Events for Single-Cell RNA-Seq Data with Limited Sample Size" https://doi.org/10.1093/bioinformatics/btz627  Bioinformatics, August 8, 2019. 
</details>

- [netNMF-sc](https://github.com/raphael-group/netNMF-sc) - scRNA-seq nonnegative matrix factorization for imputation and dimensionality reduction for improved clustering. Uses gene-gene interaction network to constrain W gene matrix on prior knowledge (graph regularized NMF). Added penalization for dropouts. Tested on simulated and experimental data, compared with several imputation and clustering methods. <details>
    <summary>Paper</summary>
    Elyanow, Rebecca, Bianca Dumitrascu, Barbara E Engelhardt, and Benjamin J Raphael. "NetNMF-Sc: Leveraging Gene-Gene Interactions for Imputation and Dimensionality Reduction in Single-Cell Expression Analysis" https://doi.org/10.1101/544346  BioRxiv, February 8, 2019.
</details>

- [scRMD](https://github.com/ChongC1990/scRMD) - dropout imputation in scRNA-seq via robust matrix decomposition into true expression matrix (further decomposed into a matrix of means and gene's random deviation from its mean) minus dropout matrix plus error matrix. A function to estimate the matrix of means and dropouts. Comparison with MAGIC, scImpute. <details>
    <summary>Paper</summary>
    Chen, Chong, Changjing Wu, Linjie Wu, Yishu Wang, Minghua Deng, and Ruibin Xi. "ScRMD: Imputation for Single Cell RNA-Seq Data via Robust Matrix Decomposition" https://doi.org/10.1101/459404 November 4, 2018
</details>

- [SAVER](https://github.com/mohuangx/SAVER) (single-cell analysis via expression recovery) - scRNA-seq imputation (UMI matrix) utilizing gene-to-gene relationship. Recover missing gene expression, removes technical variation. Assumes gene counts follow a negative binomial distribution, estimates the prior parameters in an empirical Bayes-like approach with as Poisson LASSO regression, using the expression of other genes as predictors. Tested using RNA FISH data as a reference, better recover gene expression using Drop-seq data. Outperforms [MAGIC](#magic) and [scImpute](#scimpute). <details>
    <summary>Paper</summary>
    Huang, Mo, Jingshu Wang, Eduardo Torre, Hannah Dueck, Sydney Shaffer, Roberto Bonasio, John I. Murray, Arjun Raj, Mingyao Li, and Nancy R. Zhang. “SAVER: Gene Expression Recovery for Single-Cell RNA Sequencing.” Nature Methods 15, no. 7 (July 2018): 539–42. https://doi.org/10.1038/s41592-018-0033-z.
</details>

- [netSmooth](https://github.com/BIMSBbioinfo/netSmooth) - network diffusion-based method that uses priors for the covariance structure of gene expression profiles to smooth scRNA-seq experiments. Incorporates prior knowledge (i.e. protein-protein interaction networks) for imputation. Note that dropout applies to whole transcriptome. Compared with MAGIC, scImpute. Improves clustering, biological interpretation. <details>
    <summary>Paper</summary>
    Ronen, Jonathan, and Altuna Akalin. "NetSmooth: Network-Smoothing Based Imputation for Single Cell RNA-Seq" https://doi.org/10.12688/f1000research.13511.3  F1000Research 7 (July 10, 2018)
</details>

- [DCA](https://github.com/theislab/dca) - A deep count autoencoder network to denoise scRNA-seq data. Zero-inflated negative binomial model. Current approaches - scimpute, MAGIC, SAVER. Benchmarking by increased correlation between bulk and scRNA-seq data, between protein and RNA levels, between key regulatory genes, better DE concordance in bulk and scRNA-seq, improved clustering. <details>
    <summary>Paper</summary>
    Eraslan, Gökcen, Lukas M. Simon, Maria Mircea, Nikola S. Mueller, and Fabian J. Theis. "Single Cell RNA-Seq Denoising Using a Deep Count Autoencoder" https://doi.org/10.1101/300681  April 13, 2018. 
</details>

- [kNN-smoothing](https://github.com/yanailab/knn-smoothing) of scRNA-seq data, aggregates information from similar cells, improves signal-to-noise ratio. Based on observation that gene expression in technical replicates are Poisson distributed. Freeman-Tukey transform to minimize variability of low expressed genes. Tested using real and simulated data. Improves clustering, PCA, Selection of k is critical, discussed. <details>
    <summary>Paper</summary>
    Wagner, Florian, Yun Yan, and Itai Yanai. "K-Nearest Neighbor Smoothing for High-Throughput Single-Cell RNA-Seq Data" https://doi.org/10.1101/217737  BioRxiv, April 9, 2018. 
</details>

- <a name="scimpute">[scImpute](https://github.com/Vivianstats/scImpute)</a> - imputation of scRNA-seq data. Methodology: 1) Determine K subpopulations using PCA, remove outliers; 2) Mixture model of gene i in subpopulation k as gamma and normal distributions, estimate dropout probability d; 3) Impute dropout values by splitting the subpopulation into A (dropout larger than threshold t) and B (smaller). Information from B is used to impute A. Better than MAGIC, SAVER. <details>
    <summary>Paper</summary>
    Li, Wei Vivian, and Jingyi Jessica Li. "An Accurate and Robust Imputation Method ScImpute for Single-Cell RNA-Seq Data" https://doi.org/10.1038/s41467-018-03405-7  Nature Communications 9, no. 1 (08 2018)

- [LATE](https://github.com/audreyqyfu/LATE) - Learning with AuToEncoder to imputescRNA-seq data. `TRANSLATE` (TRANSfer learning with LATE) uses reference (sc)RNA-seq dataset to learn initial parameter estimates. TensorFlow implementation for GPU and CPU. ReLu as an activation function. Various optimization techniques. Comparison with MAGIC, scVI, DCA, SAVER. Links to data. <details>
    <summary>Paper</summary>
    Badsha, Md. Bahadur, Rui Li, Boxiang Liu, Yang I. Li, Min Xian, Nicholas E. Banovich, and Audrey Qiuyan Fu. "Imputation of Single-Cell Gene Expression with an Autoencoder Neural Network" https://doi.org/10.1101/504977  BioRxiv, January 1, 2018
</details>

- <a name="magic">[MAGIC](https://github.com/KrishnaswamyLab/MAGIC)</a> - Markov Affinity-based Graph Imputation of Cells. Only \~5-15% of scRNA-seq data is non-zero, the rest are drop-outs. Use the diffusion operator to discover the manifold structure and impute gene expression. Detailed methods description. In real (bone marrow and retinal bipolar cells) and synthetic datasets, Imputed scRNA-seq data clustered better, enhances gene interactions, restores expression of known surface markers, trajectories. scRNA-seq data is preprocessed by library size normalization and PCA (to retain 70% of variability). Comparison with SVD-based low-rank data approximation (LDA) and Nuclear-Norm-based Matrix Completion (NNMC). [GitHub](https://github.com/KrishnaswamyLab/MAGIC). <details>
    <summary>Paper</summary>
    Van Dijk, David, Roshan Sharma, Juozas Nainys, Kristina Yim, Pooja Kathail, Ambrose J. Carr, Cassandra Burdziak et al. "Recovering gene interactions from single-cell data using data diffusion." Cell 174, no. 3 (2018): 716-729. https://doi.org/10.1016/j.cell.2018.05.061
</details>


## Dimensionality reduction

- [MultiMAP](https://github.com/Teichlab/MultiMAP) dimensionality reduction algorithm. Works with dataset-specific features (does not require features to be shared across datasets, e.g., 20K-gene scRNA-seq and 100K-peak scATAC-seq datasets). Generalizes the UMAP algorithm to data with different dimensions, constructs a nonlinear manifold, constructs a joint graph on the manifold (MultiGraph), cross-entropy minimization to optimize the low-dimensional embedding of the manifold and data. Allows to specify the influence of each dataset on the embedding. Tested on synthetic and experimental data, including spatial transcriptomics datasets, outperforms Seurat, LIGER, iNMF, Conos, GLUER, significantly faster and scalable. <details>
    <summary>Paper</summary>
    Jain, M.S., Polanski, K., Conde, C.D. et al. MultiMAP: dimensionality reduction and integration of multimodal data. Genome Biol 22, 346 (2021). https://doi.org/10.1186/s13059-021-02565-y
</details>

- [Poincare maps](https://github.com/facebookresearch/PoincareMaps) for two-dimensional scRNA-seq data representation. Preserves local and global distances, hierarchy, the center of the Poincare disk can be considered as a root node. Three-step procedure: 1) k-nearest-neighbor graph, 2) global geodesic distances from the kNN graph, 3) two-dimensional embeddings in the Poincare disk with hyperbolic distances preserve the inferred geodesic distances. Compared with t-SNE, UMAP, PCA, Monocle 2, SAUCIE and several other visualization and lineage detection methods. Two metrics to compare embeddings, Qlocal and Qglobal. References to several public datasets used for reanalysis. <details>
    <summary>Paper</summary>
    Klimovskaia, Anna, David Lopez-Paz, Léon Bottou, and Maximilian Nickel. "Poincaré Maps for Analyzing Complex Hierarchies in Single-Cell Data" https://doi.org/10.1038/s41467-020-16822-4  Nature Communications 11, no. 1 (December 2020)
</details>

- [scHPF](https://github.com/simslab/schpf) - single-cell hierarchical Poisson Factorization for discovering patterns of gene expressions and cells. A Bayesian factorization method, does not require normalization, explicitly models sparsity across cells and genes. Compared with PCA, NMF, FA, ZIFA, ZINB-WaVE on three datasets, it better captures statistical and biological properties of scRNA-seq data. Python implementation. <details>
    <summary>Paper</summary>
    Levitin, Hanna Mendes, Jinzhou Yuan, Yim Ling Cheng, Francisco JR Ruiz, Erin C Bush, Jeffrey N Bruce, Peter Canoll, et al. "De Novo Gene Signature Identification from Single‐cell RNA‐seq with Hierarchical Poisson Factorization" https://doi.org/10.15252/msb.2018855  Molecular Systems Biology 15, no. 2 (February 2019)
</details>

- [SAUCIE](https://github.com/KrishnaswamyLab/SAUCIE) - deep neural network with regularization on layers to improve interpretability. Denoising, batch removal, imputation, visualization of low-dimensional representation. Extensive comparison on simulated and real data. <details>
    <summary>Paper</summary>
    Amodio, Matthew, David van Dijk, Krishnan Srinivasan, William S Chen, Hussein Mohsen, Kevin R Moon, Allison Campbell, et al. "Exploring Single-Cell Data with Deep Multitasking Neural Networks" https://doi.org/10.1101/237065  August 27, 2018. 
</details>

- [UMAP](http://github.com/lmcinnes/umap) - Uniform Manifold Approximation and Projection, dimensionality reduction using machine learning. Detailed statistical framework. Compared with t-SNE, better preserves global structure. [R implementation](https://github.com/jlmelville/uwot). <details>
    <summary>Paper</summary>
    McInnes, Leland, and John Healy. "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction](http://arxiv.org/abs/1802.03426  ArXiv:1802.03426 [Cs, Stat], February 9, 2018
</details>

- [CIDR](https://github.com/VCCRI/CIDR) - Clustering through Imputation and Dimensionality Reduction. Impute dropouts. Explicitly deconvolve Euclidean distance into distance driven by complete, partially complete, and dropout pairs. Principal Coordinate Analysis. <details>
    <summary>Paper</summary>
    Lin, Peijie, Michael Troup, and Joshua W. K. Ho. "CIDR: Ultrafast and Accurate Clustering through Imputation for Single-Cell RNA-Seq Data" https://doi.org/10.1186/s13059-017-1188-0  Genome Biology 18, no. 1 (December 2017). 
</details>

- [VASC](https://github.com/wang-research/VASC) - deep variational autoencoder for scRNA-seq data for dimensionality reduction and visualization. Tested on twenty datasets vs PCA, tSNE, ZIFA, and SIMLR. Four metrics to assess clustering performance: NMI (normalized mutual information score), ARI (adjusted rand index), HOM (homogeneity) and COM (completeness). No filtering, only log transformation. Keras implementation. [Datasets](https://hemberg-lab.github.io/scRNA.seq.datasets/). <details>
    <summary>Paper</summary>
    Wang, Dongfang, and Jin Gu. "VASC: Dimension Reduction and Visualization of Single Cell RNA Sequencing Data by Deep Variational Autoencoder" https://doi.org/10.1101/199315  October 6, 2017.
</details>

- [ZINB-WAVE](https://bioconductor.org/packages/zinbwave/) - Zero-inflated negative binomial model for normalization, batch removal, and dimensionality reduction. Extends the RUV model with more careful definition of "unwanted" variation as it may be biological. Good statistical derivations in Methods. Refs to real and simulated scRNA-seq datasets. <details>
    <summary>Paper</summary>
    Risso, Davide, Fanny Perraudeau, Svetlana Gribkova, Sandrine Dudoit, and Jean-Philippe Vert. "ZINB-WaVE: A General and Flexible Method for Signal Extraction from Single-Cell RNA-Seq Data" https://doi.org/10.1101/125112  BioRxiv, January 1, 2017. 
</details>

- [RobustAutoencoder](https://github.com/zc8340311/RobustAutoencoder) - Autoencoder and robust PCA for gene expression representation, robust to outliers. Main idea - split the input data X into two parts, L (reconstructed data) and S (outliers and noise). Grouped "l2,1" norm - an l2 regularizer within a group and then an l1 regularizer between groups. Iterative procedure to obtain L and S. TensorFlow implementation. <details>
    <summary>Paper</summary>
    Zhou, Chong, and Randy C. Paffenroth. "Anomaly Detection with Robust Deep Autoencoders" https://doi.org/10.1145/3097983.3098052  In Proceedings of the 23rd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining  - KDD ’17, 665–74. Halifax, NS, Canada: ACM Press, 2017. 
</details>

- [ZIFA](https://github.com/epierson9/ZIFA) - Zero-inflated dimensionality reduction algorithm for single-cell data. Single-cell dimensionality reduction. Model dropout rate as double exponential, give less weights to these counts. EM algorithm that incorporates imputation step for the expected gene expression level of drop-outs. <details>
    <summary>Paper</summary>
    Pierson, Emma, and Christopher Yau. "ZIFA: Dimensionality Reduction for Zero-Inflated Single-Cell Gene Expression Analysis" https://doi.org/10.1186/s13059-015-0805-z  Genome Biology 16 (November 2, 2015)
</details>


## Clustering

- [opt-SNE](https://github.com/omiq-ai/Multicore-opt-SNE) - data-driven automated parameter selection for t-SNE clustering. Utilizes Kullback-Leibler divergence evaluation in real time to tailor the early exaggeration and overall number of gradient descent iterations. Evaluated on flow cytometry data. C++/Python implementation. <details>
    <summary>Paper</summary>
    Belkina, Anna C., Christopher O. Ciccolella, Rina Anno, Richard Halpert, Josef Spidlen, and Jennifer E. Snyder-Cappione. “Automated Optimized Parameters for T-Distributed Stochastic Neighbor Embedding Improve Visualization and Analysis of Large Datasets.” Nature Communications 10, no. 1 (December 2019): 5415. https://doi.org/10.1038/s41467-019-13055-y.
</details>

- [scSSA](https://github.com/houtongshuai123/scSSA/) - scRNA-seq clustering based on autoencoder for dimensionality reduction/denoising (improves performance), FastICA to make the data 2D, Caussian mixture clustering. Outperforms Seurat, CIDR, and other methods on datasets from [Hemberg Lab](https://hemberg-lab.github.io/scRNA.seq.datasets/). <details>
    <summary>Paper</summary>
    Zhao, Jian-Ping, Tong-Shuai Hou, Yansen Su, and Chun-Hou Zheng. “ScSSA: A Clustering Method for Single Cell RNA-Seq Data Based on Semi-Supervised Autoencoder.” Methods 208 (December 2022): 66–74. https://doi.org/10.1016/j.ymeth.2022.10.006.
</details>

- [Challenges in scRNA-seq clustering](https://www.rna-seqblog.com/challenges-in-unsupervised-clustering-of-single-cell-rna-seq-data/). Clustering strategies (dimensionality reduction, k-means, agglomerative/divisive hierarchical clustering, discrete vs. continuous clustering). Table 1 - summary of 15 clustering methods. <details>
    <summary>Paper</summary>
    Kiselev, Vladimir Yu, Tallulah S. Andrews, and Martin Hemberg. "Challenges in Unsupervised Clustering of Single-Cell RNA-Seq Data" https://doi.org/10.1038/s41576-018-0088-9  Nature Reviews Genetics, January 7, 2019. 
</details>

- [Recommendations to properly use t-SNE on large omics datasets](https://github.com/berenslab/rna-seq-tsne) (scRNA-seq in particular) to preserve global geometry. Overview of t-SNE, PCA, MDS, UMAP, their similarities, differences, strengths and weaknesses. PCA initialization (first two components are OK), a high learning rate of n/12, and multi-scale similarity kernels. For very large data, increase exagerration. Strategies to align new points on an existing t-SNE plot, aligning two t-SNE visualizations. Extremely fast implementation is [FIt-SNE](https://github.com/KlugerLab/FIt-SNE). [Code to illustrate the use of t-SNE](https://github.com/berenslab/rna-seq-tsne). <details>
    <summary>Paper</summary>
    Kobak, Dmitry, and Philipp Berens. "The Art of Using T-SNE for Single-Cell Transcriptomics" https://doi.org/10.1038/s41467-019-13056-x  Nature Communications 10, no. 1 (December 2019)
</details>

- [BAMM-SC](https://github.com/CHPGenetics/BAMMSC) - scRNA-seq clustering. A Bayesian hierarchical Dirichlet multinomial mixture model, accounts for batch effect, operates on raw counts. Outperforms K-means, TSCAN, Seurat corrected for batch using MNN or CCA in simulated and experimental settings. <details>
    <summary>Paper</summary>
    Sun, Zhe, Li Chen, Hongyi Xin, Yale Jiang, Qianhui Huang, Anthony R. Cillo, Tracy Tabib, et al. "A Bayesian Mixture Model for Clustering Droplet-Based Single-Cell Transcriptomic Data from Population Studies" https://doi.org/10.1038/s41467-019-09639-3  Nature Communications 10, no. 1 (December 2019)
</details>

- [Spectrum](https://cran.r-project.org/web/packages/Spectrum/) - a spectral clustering method for single- or multi-omics datasets. Self-tuning kernel that adapts to local density of the graph. Tensor product graph data integration method. Implementation of fast spectral clustering method (single dataset only). Finds optimal number of clusters using eigenvector distribution analysis. References to previous methods. Excellent methods description. Compared with M3C, CLEST, PINSplus, SNF, iClusterPlus, CIMLR, MUDAN. [GitHub](https://github.com/crj32/spectrum_manuscript). <details>
    <summary>Paper</summary>
    John, Christopher R., David Watson, Michael R. Barnes, Costantino Pitzalis, and Myles J. Lewis. "Spectrum: Fast Density-Aware Spectral Clustering for Single and Multi-Omic Data" https://doi.org/10.1093/bioinformatics/btz704  Bioinformatics (Oxford, England), September 10, 2019
</details>

- [SAUCIE](#saucie) - a regularized autoencoder for scRNA-seq data denoising, batch correction, low-dimensional representation and clustering.

- [PanoView](https://github.com/mhu10/scPanoView) - scRNA-seq iterative clustering in an evolving 3D PCA space, Ordering Local Maximum by Convex hull (OLMC) to identify clusters of varying density. PCA on most variable genes, finding most optimal largest cluster within first 3 PCs, repeat PCA for the remaining cells etc. Tested on multiple simulated and experimental scRNA-seq datasets, compared with 9 methods, the Adjusted Rand Index as performance metric. <details>
    <summary>Paper</summary>
    Hu, Ming-Wen, Dong Won Kim, Sheng Liu, Donald J. Zack, Seth Blackshaw, and Jiang Qian. "PanoView: An Iterative Clustering Method for Single-Cell RNA Sequencing Data" https://doi.org/10.1371/journal.pcbi.1007040  Edited by Qing Nie. PLOS Computational Biology 15, no. 8 (August 30, 2019)
</details>

- [FIt-SNE](https://github.com/KlugerLab/FIt-SNE) - accelerated version of t-SNE clustering for visualizing thousands/milions of cells. [t-SNE-Heatmaps](https://github.com/KlugerLab/t-SNE-Heatmaps) - discretized t-SNE clustering representation as a heatmap. [Detailed methods](https://gauss.math.yale.edu/~gcl22/blog/numerics/low-rank/t-sne/2018/01/11/low-rank-kernels.html). <details>
    <summary>Paper</summary>
    Linderman, George C., Manas Rachh, Jeremy G. Hoskins, Stefan Steinerberger, and Yuval Kluger. "Fast Interpolation-Based t-SNE for Improved Visualization of Single-Cell RNA-Seq Data" https://doi.org/10.1038/s41592-018-0308-4  Nature Methods, February 11, 2019. 
</details>

- [TooManyCells](https://github.com/GregorySchwartz/tooManyCellsR) - divisive hierarchical spectral clustering of scRNA-seq data. Uses truncated singular vector decomposition to bipartition the cells. Newman-Girvain modularity Q to assess whether bipartition is significant or should be stopped. [BirchBeer](https://github.com/faryabiLab/birch-beer) visualization. Outperforms Phenograph, Seurat, Cellranger, Monocle, the latter is second in performance. Excels for rare populations. Normalization marginally affects performance. <details>
    <summary>Paper</summary>
    Schwartz, Gregory W, Jelena Petrovic, Maria Fasolino, Yeqiao Zhou, Stanley Cai, Lanwei Xu, Warren S Pear, Golnaz Vahedi, and Robert B Faryabi. "TooManyCells Identifies and Visualizes Relationships of Single-Cell Clades" https://doi.org/10.1101/519660  BioRxiv, January 13, 2019. 
</details>

- [scClustViz](https://github.com/BaderLab/scClustViz) - assessment of scRNA-seq clustering using differential expression (Wilcoxon test) as a guide. Testing for two differences: difference in detection rate (dDR) and log2 gene expression ratio (logGER). Two hypothesis testing: one cluster vs. all, each cluster vs. another cluster. accepts SincleCellExperiment and Seurat objects (log2-transformed data), needs a data frame with different cluster assignments. Analysis within R, save as RData, visualize results in R Shiny app. <details>
    <summary>Paper</summary>
    Innes, BT, and GD Bader. "ScClustViz - Single-Cell RNAseq Cluster Assessment and Visualization" https://doi.org/10.12688/f1000research.16198.2  F1000Research 7, no. 1522 (2019). 
</details>

- [SHARP](https://github.com/shibiaowan/SHARP) - an ensemble random projection (RP)-based algorithm. Scalable, allows for clustering of 1.3 million cells (splitting the matrix into blocks, RP on each, then weighted ensemble clustering. Outperforms SC3, SIMLR, hierarchical clustering, tSNE + k-means. Tested on 17 public datasets. Robust to dropouts. Compatible with (UMI-based) counts (per million), FPKM/RPKM, TPM. Methods detailing four algorithm steps (data partition, RP, weighted ensemble clustering, similarity-based meta-clustering). <details>
    <summary>Paper</summary>
    Wan, Shibiao, Junil Kim, and Kyoung Jae Won. "SHARP: Single-Cell RNA-Seq Hyper-Fast and Accurate Processing via Ensemble Random Projection" https://doi.org/10.1101/461640  Preprint. Bioinformatics, November 4, 2018
</details>

- [Performance evaluation of 14 scRNA-seq clustering algorithms](ttps://bioconductor.org/packages/DuoClustering2018/) using nine experimental and three simulated datasets. SC3 and Seurat perform best overall. Normalized Shannon entropy, adjusted Rand index for performance evaluation. Ensemble clustering doesn't help. [R scripts](https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison) and a [data package](https://bioconductor.org/packages/DuoClustering2018/) for clustering benchmarking with preprocessed and experimental scRNA-seq datasets. <details>
    <summary>Paper</summary>
    Duò, Angelo, Mark D. Robinson, and Charlotte Soneson. "A Systematic Performance Evaluation of Clustering Methods for Single-Cell RNA-Seq Data" https://doi.org/10.12688/f1000research.15666.2  F1000Research 7 (September 10, 2018)
</details>

- [clusterExperiment](https://bioconductor.org/packages/clusterExperiment/) R package for scRNA-seq data visualization. Resampling-based Sequential Ensemble Clustering (RSEC) method. clusterMany - makeConsensus - makeDendrogram - mergeClusters pipeline. Biomarker detection by differential expression analysis between clusters. <details>
    <summary>Paper</summary>
    Risso, Davide, Liam Purvis, Russell B. Fletcher, Diya Das, John Ngai, Sandrine Dudoit, and Elizabeth Purdom. "ClusterExperiment and RSEC: A Bioconductor Package and Framework for Clustering of Single-Cell and Other Large Gene Expression Datasets" https://doi.org/10.1371/journal.pcbi.1006378  Edited by Aaron E. Darling. PLOS Computational Biology 14, no. 9 (September 4, 2018)
</details>

- [PHATE](https://github.com/KrishnaswamyLab/PHATE) (Potential of Heat-diffusion for Affinity-based Transition Embedding) - low-dimensional embedding, denoising, and visualization, applicable to scRNA-seq, microbiome, SNP, Hi-C (as affinity matrices) and other data. Preserves biological structures and branching better than PCA, tSNE, diffusion maps. Robust to noise and subsampling. Detailed methods description and graphical representation of the algorithm. [Tweetorial](https://twitter.com/KrishnaswamyLab/status/1201935823056199680?s=20). <details>
    <summary>Paper</summary>
    Moon, Kevin R., David van Dijk, Zheng Wang, Scott Gigante, Daniel Burkhardt, William Chen, Antonia van den Elzen, et al. "Visualizing Transitions and Structure for Biological Data Exploration" https://doi.org/10.1101/120378  June 28, 2018. 
</details>

- [Bisquit](https://github.com/sandhya212/BISCUIT_SingleCell_IMM_ICML_2016) - a Bayesian clustering and normalization method. <details>
    <summary>Paper</summary>
    Azizi, Elham, Ambrose J. Carr, George Plitas, Andrew E. Cornish, Catherine Konopacki, Sandhya Prabhakaran, Juozas Nainys, et al. "Single-Cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment" https://doi.org/10.1016/j.cell.2018.05.060  Cell, June 2018. 
</details>

- [scVAE](https://github.com/chgroenbech/scVAE) - Variational auroencoder frameworks for modelling raw RNA-seq counts, denoising the data to improve biologically plausible grouping in scRNA-seq data. Improvement in Rand index. <details>
    <summary>Paper</summary>
    Grønbech, Christopher Heje, Maximillian Fornitz Vording, Pascal N Timshel, Casper Kaae Sønderby, Tune Hannes Pers, and Ole Winther. "ScVAE: Variational Auto-Encoders for Single-Cell Gene Expression Data" https://doi.org/10.1101/318295  May 16, 2018. 
</details>

- [Conos](https://github.com/hms-dbmi/conos) - clustering of scRNA-seq samples by joint graph construction. Seurat or pagoda2 for data preprocessing, selection of hypervariable genes, initial clustering (KNN, or dimensionality reduction), then joint clustering. R package. <details>
    <summary>Paper</summary>
    Barkas, Nikolas, Viktor Petukhov, Daria Nikolaeva, Yaroslav Lozinsky, Samuel Demharter, Konstantin Khodosevich, and Peter V Kharchenko. “Wiring Together Large Single-Cell RNA-Seq Sample Collections.” BioRxiv, January 1, 2018. 
</details>

- [MetaCell](https://bitbucket.org/tanaylab/metacell/src/default/) - partitioning scRNA-seq data into metacells - disjoint and homogeneous/compact groups of cells exhibiting only sampling variance. Most variable genes to cell-to-cell similarity matrix (PCC on  to Knn similarity graph that is partitioned by bootstrapping to obtain subgraphs. Tested on several [10X datasets](https://support.10xgenomics.com/single-cell-gene-expression/datasets). <details>
    <summary>Paper</summary>
    Baran, Yael, Arnau Sebe-Pedros, Yaniv Lubling, Amir Giladi, Elad Chomsky, Zohar Meir, Michael Hoichman, Aviezer Lifshitz, and Amos Tanay. "MetaCell: Analysis of Single Cell RNA-Seq Data Using k-NN Graph Partitions" https://doi.org/10.1101/437665  BioRxiv, January 1, 2018
</details>

- [SIMLR](https://github.com/BatzoglouLabSU/SIMLR) - scRNA-seq dimensionality reduction, clustering, and visualization based on multiple kernel-learned distance metric. Comparison with PCA, t-SNE, ZIFA. Seven datasets. R and Matlab implementation. <details>
    <summary>Paper</summary>
    Wang, Bo, Junjie Zhu, Emma Pierson, Daniele Ramazzotti, and Serafim Batzoglou. "Visualization and Analysis of Single-Cell RNA-Seq Data by Kernel-Based Similarity Learning" https://doi.org/10.1101/460246  Nature Methods 14, no. 4 (April 2017)
</details>

- [SC3](https://bioconductor.org/packages/SC3/) - single-cell clustering. Multiple clustering iterations, consensus matrix, then hierarhical clustering. Benchmarking against other methods. <details>
    <summary>Paper</summary>
    Kiselev, Vladimir Yu, Kristina Kirschner, Michael T Schaub, Tallulah Andrews, Andrew Yiu, Tamir Chandra, Kedar N Natarajan, et al. "SC3: Consensus Clustering of Single-Cell RNA-Seq Data" https://doi.org/10.1038/nmeth.4236  Nature Methods 14, no. 5 (March 27, 2017)
</details>

- [destiny](https://bioconductor.org/packages/destiny/) - R package for diffusion maps-based visualization of single-cell data. <details>
    <summary>Paper</summary>
    Haghverdi, Laleh, Florian Buettner, and Fabian J. Theis. "Diffusion Maps for High-Dimensional Single-Cell Analysis of Differentiation Data" https://doi.org/10.1093/bioinformatics/btv325  Bioinformatics 31, no. 18 (September 15, 2015) - Introduction of other methods, Table 1 compares them. Methods details. Performance is similar to PCA and tSNE. 
</details>

- [PhenoGraph](https://github.com/jacoblevine/PhenoGraph) - discovers subpopulations in scRNA-seq data. High-dimensional space is modeled as a nearest-neighbor graph, then the Louvain community detection algorithm. No assumptions about the size, number, or form of subpopulations. <details>
    <summary>Paper</summary>
    Levine, Jacob H., Erin F. Simonds, Sean C. Bendall, Kara L. Davis, El-ad D. Amir, Michelle D. Tadmor, Oren Litvin, et al. "Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells That Correlate with Prognosis" https://doi.org/10.1016/j.cell.2015.05.047  Cell 162, no. 1 (July 2015)
</details>

- [SNN-Cliq](http://bioinfo.uncc.edu/SNNCliq/) - shared nearest neighbor clustering of scRNA-seq data, represented as a graph. Similarity between two data points based on the ranking of their shared neighborhood. Automatically determine the number of clusters, accomodates different densities and shapes. Compared with K-means and DBSCAN using Purity, Adjusted Rand Indes, F1-score. Matlab, Python, R implementation. <details>
    <summary>Paper</summary>
    Xu, Chen, and Zhengchang Su. "Identification of Cell Types from Single-Cell Transcriptomes Using a Novel Clustering Method" https://doi.org/10.1093/bioinformatics/btv088  Bioinformatics (Oxford, England) 31, no. 12 (June 15, 2015)
</details>

- [viSNE](https://cran.r-project.org/web/packages/Rtsne/) - the Barnes-Hut implementation of the t-SNE algorithm, improved and tailored for the analysis of single-cell data. [Details of tSNE](https://www.denovosoftware.com/site/manual/visne.htm), and the [Rtsne R package](https://cran.r-project.org/web/packages/Rtsne/). <details>
    <summary>Paper</summary>
    Amir, El-ad David, Kara L Davis, Michelle D Tadmor, Erin F Simonds, Jacob H Levine, Sean C Bendall, Daniel K Shenfeld, Smita Krishnaswamy, Garry P Nolan, and Dana Pe’er. "ViSNE Enables Visualization of High Dimensional Single-Cell Data and Reveals Phenotypic Heterogeneity of Leukemia" https://doi.org/10.1038/nbt.2594  Nature Biotechnology 31, no. 6 (June 2013)
</details>

- [celda](https://github.com/compbiomed/celda) - CEllular Latent Dirichlet Allocation. Simultaneous clustering of cells into subpopulations and genes into transcriptional states. [Tutorials](https://compbiomed.github.io/celda_tutorials/). No preprint yet.




### Time, trajectory inference

- [CellRank](https://cellrank.readthedocs.io/en/stable/) - single-cell fate mapping combining trajectory inference and RNA velocity directionality (scVelo), accounting for the stochastic nature of fate decisions and uncertainty in velocity vectors. Velocity alone is insufficient. Detects the initial, terminal and intermediate cell states and computes a global map of fate potentials. State transitions are modeled using a Markov chain. Stability index to automatically identify terminal states. Outperforms Palantir, STEMNET and FateID in diverse scenarious (development, regeneration, reprogramming, disease), fast, less memory, scalable. Input - (imputed) gene count matrix and velocity matrix (any vector field). Python, installable in Conda environment, Jupyter notebooks. [Tutorial](https://cellrank.readthedocs.io/en/stable/cellrank_basics.html), [Code to reproduce the results](https://github.com/theislab/cellrank_reproducibility). [Tweet](https://twitter.com/dana_peer/status/1481658478296907780?s=20) by Dana Pe'er. <details>
    <summary>Paper</summary>
    Lange, Marius, Volker Bergen, Michal Klein, Manu Setty, Bernhard Reuter, Mostafa Bakhti, Heiko Lickert, et al. “CellRank for Directed Single-Cell Fate Mapping.” Nature Methods, January 13, 2022. https://doi.org/10.1038/s41592-021-01346-6
</details>

- [STREAM](https://stream.pinellolab.partners.org/) - trajectory analysis in both single-cell transcriptomic (scRNA-seq) and epigenomic (scATAC-seq) data. Ability to map new cells on reference trajectories. Exploration of cell type composition, relevant genes, TF binding dynamics. Modified Locally Linear Embedding (MLLE, on top PC components or variable genes) which preserves distances within local neighborhoods. Infers trajectories using Elastic Principal Graph implementation ([ElPiGraph](https://github.com/sysbio-curie/ElPiGraph.R)). Allows for setting root point. Compared (four metrics) with 10 trajectory inference methods out of over 50, supplementary tables have methods details. For scATAC-seq, filter the data and focus on variable chromatin regions. [Bioconda](https://bioconda.github.io/recipes/stream/README.html) recipe, [Github](https://github.com/pinellolab/STREAM) - command-line, Jupyter notebook. <details>
    <summary>Paper</summary>
    Chen, Huidong, Luca Albergante, Jonathan Y. Hsu, Caleb A. Lareau, Giosuè Lo Bosco, Jihong Guan, Shuigeng Zhou, et al. “Single-Cell Trajectories Reconstruction, Exploration and Mapping of Omics Data with STREAM.” Nature Communications 10, no. 1 (December 2019): 1903. https://doi.org/10.1038/s41467-019-09670-4.
</details>

- [DYNO: A collection of 57 trajectory inference methods](https://github.com/dynverse/dynmethods#list-of-included-methods). <details>
    <summary>Paper</summary>
    Saelens, Wouter, Robrecht Cannoodt, Helena Todorov, and Yvan Saeys. "A Comparison of Single-Cell Trajectory Inference Methods: Towards More Accurate and Robust Tools" https://doi.org/10.1101/276907  March 5, 2018. - Review of trajectory 29 inference methods for single-cell RNA-seq (out of 57 methods collected). Slingshot, TSCAN and Monocle DDRTree perform best overall. https://github.com/dynverse/dynverse
</details>

- [Single-cell RNA-seq pseudotime estimation algorithms](https://github.com/agitter/single-cell-pseudotime), by Anthony Gitter. References and descriptions of many algorithms.

- [Supplementary Table 1 of comparison of 11 trajectory inference methods](Pseudotime_tools.png) from Huidong Chen et al., "Single-Cell Trajectories Reconstruction, Exploration and Mapping of Omics Data with STREAM" https://doi.org/10.1038/s41467-019-09670-4  Nature Communications 10, no. 1 (April 23, 2019)

- [LACE](https://github.com/BIMIB-DISCo/LACE) - R package for Longitudinal Analysis of Cancer Evolution using mutations called from scRNA-seq data. Input - binary matrix with 1 indicating somatic variants, called using GATK. Boolean matrix factorization solved via exhaustive search or via MCMC. Output - the maximum likelihood clonal tree describing the longitudinal evolution of a tumor. Compared with CALDER, SCITE, TRaIT. <details>
    <summary>Paper</summary>
    Ramazzotti, Daniele, Fabrizio Angaroni, Davide Maspero, Gianluca Ascolani, Isabella Castiglioni, Rocco Piazza, Marco Antoniotti, and Alex Graudenzi. "Longitudinal Cancer Evolution from Single Cells" https://doi.org/10.1101/2020.01.14.906453  Preprint. Bioinformatics, January 15, 2020.
</details>

- [Tempora](https://github.com/BaderLab/Tempora) - pseudotime reconstruction for scRNA-seq data, considering time series information. Matrix of scRNA-seq gene expression to clusters to pathways (single-cell enrichment), to graph (ARACNE), incorporate time scores, refine cell types. Outperforms Monocle2, TSCAN. <details>
    <summary>Paper</summary>
    Tran, Thinh N., and Gary D. Bader. "Tempora: Cell Trajectory Inference Using Time-Series Single-Cell RNA Sequencing Data" https://doi.org/10.1101/846907  Preprint. Bioinformatics, November 18, 2019.
</details>

- [STREAM](https://github.com/pinellolab/STREAM) (Single-cell Trajetrories Reconstruction, Exploration And Mapping) - reconstruction of pseudotime trajectory with branching from scRNA-seq or scATAC-seq data. Principal graph, modified locally linear embedding, Elastic Principal Graph (previously published). Illustrated on several datasets. Compared with ten methods (Monocle2, scTDA, Wishbone, TSCAN, SLICER, DPT, GPFates, Mpath, SCUBA, PHATE). [Web-version, visualization of published datasets](http://stream.pinellolab.org/), [Jupyter notebooks, Docker container](https://github.com/pinellolab/STREAM). <details>
    <summary>Paper</summary>
    Chen, Huidong, Luca Albergante, Jonathan Y. Hsu, Caleb A. Lareau, Giosuè Lo Bosco, Jihong Guan, Shuigeng Zhou, et al. "Single-Cell Trajectories Reconstruction, Exploration and Mapping of Omics Data with STREAM" https://doi.org/10.1038/s41467-019-09670-4  Nature Communications 10, no. 1 (April 23, 2019)
</details>

- [Monocle 2](https://github.com/cole-trapnell-lab/monocle-release) - Reversed graph embedding (DDRTree), finding low-dimensional mapping of differential genes while learning the graph in this reduced space. Allows for the selection of root. Compared with Monocle 1, Wishbone, Diffusion Pseudotime, SLICER. Code, https://github.com/cole-trapnell-lab/monocle-release, [analysis scripts](https://github.com/cole-trapnell-lab/monocle2-rge-paper). <details>
    <summary>Paper</summary>
    Qiu, Xiaojie, Qi Mao, Ying Tang, Li Wang, Raghav Chawla, Hannah A Pliner, and Cole Trapnell. "Reversed Graph Embedding Resolves Complex Single-Cell Trajectories" https://doi.org/10.1038/nmeth.4402  Nature Methods 14, no. 10 (August 21, 2017)
</details>

- [Slingshot](https://github.com/kstreet13/slingshot) - Inferring multiple developmental lineages from single-cell gene expression. Clustering by gene expression, then inferring cell lineage as an ordered set of clusters -minimum spanning tree through the clusters using Mahalanobis distance. Initial state and terminal state specification. Principal curves to draw a path through the gene expression space of each lineage. <details>
    <summary>Paper</summary>
    Street, Kelly, Davide Risso, Russell B Fletcher, Diya Das, John Ngai, Nir Yosef, Elizabeth Purdom, and Sandrine Dudoit. "Slingshot: Cell Lineage and Pseudotime Inference for Single-Cell Transcriptomics" https://doi.org/10.1186/s12864-018-4772-0  BioRxiv, January 1, 2017. 
</details>

- [SCITE](https://github.com/cbg-ethz/SCITE) - a stochastic search algorithm to identify the evolutionary history of a tumor from mutation patterns in scRNA-seq data. MCMC to compute the maximum-likelihood mutation history. Accounts for noise and dropouts. Input - Boolean mutation matrix, output - maximum-likelihood-inferred mutation tree. Compared with Kim & Simon approach, BitPhylogeny. <details>
    <summary>Paper</summary>
    Jahn, Katharina, Jack Kuipers, and Niko Beerenwinkel. "Tree Inference for Single-Cell Data" https://doi.org/10.1186/s13059-016-0936-x  Genome Biology 17, no. 1 (December 2016)
</details>

- [DPT](https://theislab.github.io/destiny/) - diffusion pseudotime, arrange cells in the pseudotemporal order. Random-walk-based distance that is computed based on Euclidean distance in the diffusion map space. Weighted nearest neighborhood of the data, probabilities of transitioning to each other cell using random walk, DTP is the euclidean distance between the two vectors, stored in a transition matrix. Robust to noise and sparsity. Method compared with Monocle, Wishbone, Wanderlust. <details>
    <summary>Paper</summary>
    Haghverdi, Laleh, Maren Büttner, F. Alexander Wolf, Florian Buettner, and Fabian J. Theis. "Diffusion Pseudotime Robustly Reconstructs Lineage Branching" https://doi.org/10.1038/nmeth.3971  Nature Methods 13, no. 10 (2016)
</details>

- [TSCAN](https://bioconductor.org/packages/TSCAN/) - pseudo-time reconstruction for scRNA-seq. Clustering first, then minimum spanning tree over cluster centers. Cells are projected to the tree (PCA) to determine their pseudo-time and order. R package that includes GUI. <details>
    <summary>Paper</summary>
    Ji, Zhicheng, and Hongkai Ji. "TSCAN: Pseudo-Time Reconstruction and Evaluation in Single-Cell RNA-Seq Analysis" https://doi.org/10.1093/nar/gkw430  Nucleic Acids Research 44, no. 13 (27 2016)
</details>

- [Wishbone](https://github.com/ManuSetty/wishbone) - ordering scRNA-seq along bifurcating developmental trajectories. nearest-heighbor graphs to capture developmental distances using shortest paths. Solves short-circuits by low-dimensional projection using diffusion maps. Waypoints as guides for building the trajectory. Detailed and comprehensive Methods description. Supersedes Wanderlust. Comparison with SCUBA, Monocle. <details>
    <summary>Paper</summary>
    Setty, Manu, Michelle D. Tadmor, Shlomit Reich-Zeliger, Omer Angel, Tomer Meir Salame, Pooja Kathail, Kristy Choi, Sean Bendall, Nir Friedman, and Dana Pe’er. "Wishbone Identifies Bifurcating Developmental Trajectories from Single-Cell Data" https://doi.org/10.1038/nbt.3569  Nature Biotechnology 34, no. 6 (2016)
</details>

- [cellTree](https://bioconductor.org/packages/cellTree/) - hierarchical tree inference and visualization. Latent Dirichlet Allocation (LDA). Cells are analogous to text documents, discretized gene expression levels replace word frequencies. The LDA model represents topic distribution for each cell, analogous to low-dimensional embedding of the data where similarity is measured with chi-square distance. Fast and precise. <details>
    <summary>Paper</summary>
    duVerle, David A., Sohiya Yotsukura, Seitaro Nomura, Hiroyuki Aburatani, and Koji Tsuda. "CellTree: An R/Bioconductor Package to Infer the Hierarchical Structure of Cell Populations from Single-Cell RNA-Seq Data" https://doi.org/10.1186/s12859-016-1175-6  BMC Bioinformatics 17, no. 1 (December 2016). 
</details>

- [Monocle](https://cole-trapnell-lab.github.io/monocle-release/) - Temporal ordering of single cell gene expression profiles. Independent Component Analysis to reduce dimensionality, Minimum Spanning Tree on the reduced representation and the longest path through it. <details>
    <summary>Paper</summary>
    Trapnell, Cole, Davide Cacchiarelli, Jonna Grimsby, Prapti Pokharel, Shuqiang Li, Michael Morse, Niall J. Lennon, Kenneth J. Livak, Tarjei S. Mikkelsen, and John L. Rinn. "The Dynamics and Regulators of Cell Fate Decisions Are Revealed by Pseudotemporal Ordering of Single Cells" https://doi.org/10.1038/nbt.2859  Nature Biotechnology 32, no. 4 (April 2014)
</details>

- [SCUBA](https://github.com/gcyuan/SCUBA) - single-cell clustering using bifurcation analysis. Cells may differentiate in a monolineage manner or may differentiate into multiple cell lineages, which is the bifurcation event - two new lineages. Methods. Matlab code. <details>
    <summary>Paper</summary>
    Marco, Eugenio, Robert L. Karp, Guoji Guo, Paul Robson, Adam H. Hart, Lorenzo Trippa, and Guo-Cheng Yuan. "Bifurcation Analysis of Single-Cell Gene Expression Data Reveals Epigenetic Landscape" https://doi.org/10.1073/pnas.1408993111  Proceedings of the National Academy of Sciences of the United States of America 111, no. 52 (December 30, 2014)
</details>


### Networks

- [hdWGCNA](https://github.com/smorabit/hdWGCNA) - high-dimensional WGCNA for co-expression network analysis in scRNA-seq data. Integrates with Seurat framework. Modules for gene module identification, network visualization, enrichment analysis, more. The concept of "metacell" (combining highly similar scRNA-seq profiles) to reduce sparsity. Module eigengenes (first PCs of the module's gene expression matrix), or alternatives (UCell, Seurat's algoriths). "Metaspots" for spatial transcriptomics. The EGAD neighbor-voting algorithm to predict pathways associated with gene modules. Applied to spatial transcriptomics of the mouse brain, isoform coexpression, autism spectrum disorders, Alzheimer disease. [Detailed documentation and vignettes](https://smorabit.github.io/hdWGCNA/). <details>
    <summary>Paper</summary>
    Morabito, Samuel, Fairlie Reese, Negin Rahimzadeh, Emily Miyoshi, and Vivek Swarup. “HdWGCNA Identifies Co-Expression Networks in High-Dimensional Transcriptomics Data.” Cell Reports Methods 3, no. 6 (June 2023): 100498. https://doi.org/10.1016/j.crmeth.2023.100498.
</details>

- [scHumanNet](https://github.com/netbiolab/scHumanNet) - scRNA-seq reference-guided (HumanNet-XC) network analysis (SCINET algorithm). Data preprocessing using the ACTIONet package (imputation, transformation, normalization). Compared with five reference-free methods (rawPCC, MetaCell, SAVER, GRNboost2, bigSCale2) and one reference-guided method (SCINET). Better defines cell type-specific genes, improves network properties. Differential network analysis. Visualization. Conda, R. <details>
    <summary>Paper</summary>
    Cha, Junha, Jiwon Yu, Jae-Won Cho, Martin Hemberg, and Insuk Lee. “ScHumanNet: A Single-Cell Network Analysis Platform for the Study of Cell-Type Speciﬁcity of Disease Genes,” n.d., 14.
</details>

- Benchmarking of six scRNA-seq network inference methods (GENIE3, GRNBoost2, PPCOR, PIDC, GeneNet, CLR) based on their reproducibility (percentage of intersection, weighted Jaccard similarity) when applied to two independent datasets for the same biilogical conditions (human retina, T-cells in colorectal cancer, human hematopoiesis, Table 1). Other benchmarking studies: [Chen and Mar 2018" https://doi.org/10.1126/science.aam8999) and [Pratapa 2020" https://doi.org/10.1038/s41592-019-0690-6). Brief description of each method. For large (up to 100,000 links) networks, GENIE2 and GRNBoost2 perform best. For smaller (100-1000 links), GRNBoost2 and CLR perform best. [scNET Jupyter notebook](https://github.com/ComputationalSystemsBiology/scNET) implementing all analyses. <details>
    <summary>Paper</summary>
    Kang, Yoonjee, Denis Thieffry, and Laura Cantini. “Evaluating the Reproducibility of Single-Cell Gene Regulatory Network Inference Algorithms.” Frontiers in Genetics 12 (March 22, 2021): 617282. https://doi.org/10.3389/fgene.2021.617282. 
    Cantini, Laura, Pooya Zakeri, Celine Hernandez, Aurelien Naldi, Denis Thieffry, Elisabeth Remy, and Anaïs Baudot. "Benchmarking Joint Multi-Omics Dimensionality Reduction Approaches for the Study of Cancer" https://doi.org/10.1038/s41467-020-20430-7  Nature Communications, (December 2021)
</details>

- Benchmarking of 11 scRNA-seq network inference methods. Top performers (PEARSON, PIDC, MERLIN, SCENIC), middle (Inferelator, SCODE, LEAP, Scribe) and bottom (knnDREMI, SILGGM). Simple correlation works well. Imputation did not benefit network inference, Human, mouse, yeast data, using scRNA-seq and bulk data (minimal performance differences). Brief description of methods, gold standard, evaluation metrics. <details>
    <summary>Paper</summary>
    Stone, Matthew, Sunnie Grace McCalla, Alireza Fotuhi Siahpirani, Viswesh Periyasamy, Junha Shin, and Sushmita Roy. "Identifying Strengths and Weaknesses of Methods for Computational Network Inference from Single Cell RNA-Seq Data" https://doi.org/10.1101/2021.06.01.446671  Preprint. Bioinformatics, June 2, 2021.
</details>

- [PAGA](https://github.com/theislab/paga) - graph-like representation of scRNA-seq data. The kNN graph is partitioned using Louvain community detection algorithm, discarding spurious edged (denoising). Much faster than UMAP. Part of Scanpy pipeline. <details>
    <summary>Paper</summary>
    Wolf, F. Alexander, Fiona K. Hamey, Mireya Plass, Jordi Solana, Joakim S. Dahlin, Berthold Göttgens, Nikolaus Rajewsky, Lukas Simon, and Fabian J. Theis. "PAGA: Graph Abstraction Reconciles Clustering with Trajectory Inference through a Topology Preserving Map of Single Cells" https://doi.org/10.1186/s13059-019-1663-x  Genome Biology 20, no. 1 (March 19, 2019)
</details>

- [SCIRA](https://bioconductor.org/packages/SEPIRA/) - infer tissue-specific regulatory networks using large-scale bulk RNA-seq, estimate regulatory activity. SEPIRA uses a greedy partial correlation framework to infer a regulatory network from GTeX data, TF-specific regulons used as target profiles in a linear regression model framework. Compared against SCENIC. Works even for small cell populations. Tested on three scRNA-seq datasets. A part of SEPIRA R package. <details>
    <summary>Paper</summary>
    Wang, Ning, and Andrew E Teschendorff. "Leveraging High-Powered RNA-Seq Datasets to Improve Inference of Regulatory Activity in Single-Cell RNA-Seq Data" https://doi.org/10.1101/553040  BioRxiv, February 22, 2019. 
</details>

- [GraphDDP](https://github.com/fabriziocosta/GraphEmbed) - combines user-guided clustering and transition of differentiation processes between clusters. Shortcomings of PCA, MDS, t-SNE. Tested on several datasets to improve interpretability of clustering, compared with other methods (Monocle2, SPRING, TSCAN). Detailed methods. <details>
    <summary>Paper</summary>
    Costa, Fabrizio, Dominic Grün, and Rolf Backofen. "GraphDDP: A Graph-Embedding Approach to Detect Differentiation Pathways in Single-Cell-Data Using Prior Class Knowledge" https://doi.org/10.1038/s41467-018-05988-7  Nature Communications 9, no. 1 (December 2018)
</details>

- [SCENIC](https://github.com/aertslab/SCENIC) - single-cell network reconstruction and cell-state identification. Three modules: 1) [GENIE3](https://github.com/aertslab/GENIE3) - connect co-expressed genes and TFs using random forest regression; 2) RcisTarget - Refine them using cis-motif enrichment; 3) [AUCell](https://github.com/aertslab/AUCell) - assign activity scores for each network in each cell type. The R implementation of GENIE3 does not scale well with larger datasets; use arboreto instead, which is a much faster Python implementation. <details>
    <summary>Paper</summary>
    Aibar, Sara, Carmen Bravo González-Blas, Thomas Moerman, Vân Anh Huynh-Thu, Hana Imrichova, Gert Hulselmans, Florian Rambow, et al. "SCENIC: Single-Cell Regulatory Network Inference and Clustering" https://doi.org/10.1038/nmeth.4463  Nature Methods 14, no. 11 (November 2017)
    Davie et al. "[A Single-Cell Transcriptome Atlas of the Aging Drosophila Brain" https://doi.org/10.1016/j.cell.2018.05.057  Cell, 2018 
</details>

- [pySCENIC, Python fast implementation of single-cell network reconstruction, detailed computational protocol](https://github.com/aertslab/SCENICprotocol). `loom` format for scRNA-seq data, minimal filtering, network inference, module generation, motif enrichment, interactive visualization via SCope. Docker/Singularity images available. Nextflow implementation. <details>
    <summary>Paper</summary>
    Van de Sande, Bram, Christopher Flerin, Kristofer Davie, Maxime De Waegeneer, Gert Hulselmans, Sara Aibar, Ruth Seurinck, et al. "A Scalable SCENIC Workflow for Single-Cell Gene Regulatory Network Analysis" https://doi.org/10.1038/s41596-020-0336-2  Nature Protocols 15, no. 7 (July 2020)
</details>

- [SINCERA](https://research.cchmc.org/pbge/sincera.html) - identification of major cell types, the corresponding gene signatures and transcription factor networks. Pre-filtering (expression filter, cell specificity filter) improves inter-sample correlation and decrease inter-sample distance. Normalization: per-sample z-score, then trimmed mean across cells. Clustering (centered Pearson for distance, average linkage), and other metrics, permutation to assess clustering significance. Functional enrichment, cell type enrichment analysis, identification of cell signatures. TF networks and their parameters (disruptive fragmentation centrality, disruptive connection centrality, disruptive distance centrality). Example analysis of mouse lung cells at E16.5, Fluidigm, 9 clusters, comparison with SNN-Cliq, scLVM, SINGuLAR Analysis Toolset. Web-site: https://research.cchmc.org/pbge/sincera.html, [GitHub](https://github.com/xu-lab/SINCERA), [Data](https://lungmap.net/breath-entity-page/?entityType=none&entityId=&entityLabel=&experimentTypes[]=LMXT0000000016). <details>
    <summary>Paper</summary>
    Guo, Minzhe, Hui Wang, S. Steven Potter, Jeffrey A. Whitsett, and Yan Xu. "SINCERA: A Pipeline for Single-Cell RNA-Seq Profiling Analysis" https://doi.org/10.1371/journal.pcbi.1004575  PLoS Computational Biology 11, no. 11 (November 2015)
</details>

### RNA velocity

- [RNA velocity software](https://docs.google.com/spreadsheets/d/1IPe2ozb1Mny8sLvJaSE57RJr3oruiBoSudAVhSH-O8M/edit#gid=586567397) - Google doc collecting references, by @AlbertVilella. [Tweet](https://twitter.com/AlbertVilella/status/1400396034505908226?s=20)

- [VeloViz](https://github.com/JEFworks-Lab/veloviz) - scRNA-seq 2D embedding using RNA velocity. R package, several tutorials. <details>
    <summary>Paper</summary>
    Atta L, Sahoo A, Fan J. [VeloViz: RNA-velocity informed embeddings for visualizing cellular trajectories" https://doi.org/10.1093/bioinformatics/btab653), Bioinformatics, 10 September 2021
</details>

- [RNA velocity estimation is affected by spliced/unspliced quantification approaches, exon/intron definitions](https://github.com/csoneson/rna_velocity_quant). Technical insights into RNA velosity calculations. Tested alevin, kallisto/bus, starsolo, cellranger for spliced/unspliced abundance estimation, on four public 10X genomics datasets. Large differences in inferred velocities, affecting biological interpretation. alevin_sep_gtr (index built on spliced transcripts and introns, introns extracted separately for each transcript isoform. Joint exonic/intronic quantification) yields results most aligned with expectations. [Code](https://github.com/csoneson/rna_velocity_quant). <details>
    <summary>Paper</summary>
    Soneson, Charlotte, Avi Srivastava, Rob Patro, and Michael B Stadler. "Preprocessing Choices Affect RNA Velocity Results for Droplet ScRNA-Seq Data" https://doi.org/10.1101/2020.03.13.990069  Preprint. Bioinformatics, March 13, 2020.
</details>

- [RNA velocity in bulk RNA-seq data](https://github.com/praneet1988/Inferring-and-Visualizing-RNA-Velocity-in-Bulk-RNA-SEQ)

- [scVelo](https://scvelo.readthedocs.io/) - more precise estimation of RNA velocity by solving the full transcriptional dynamics of splicing kinetics using a likelihood-based dynamical model. Description of steady-state (original), dynamical, and stochastic models. Ten-fold faster. 

- [velociraptor](https://bioconductor.org/packages/velociraptor/) - an R wrapper of scVelo. [Tweet](https://twitter.com/BABSBioinformat/status/1339989828742594564?s=20). <details>
    <summary>Paper</summary>
    Bergen, Volker, Marius Lange, Stefan Peidli, F. Alexander Wolf, and Fabian J. Theis. "Generalizing RNA Velocity to Transient Cell States through Dynamical Modeling" https://doi.org/10.1101/820936  Preprint. Bioinformatics, October 29, 2019. 
</details>

- [velocyto](http://velocyto.org/) - RNA velocity, the time derivative of the gene expression state, estimated by the balance of spliced and unspliced mRNAs, and the mRNA degradation, in scRNA-seq (10X, inDrop, SMART-seq2, STRT/C1 protocols). Demonstrated on several datasets. [Brief video tutorial](https://youtu.be/EPTgF4EA2zY). Python and R implementation. <details>
    <summary>Paper</summary>
    La Manno, Gioele, Ruslan Soldatov, Amit Zeisel, Emelie Braun, Hannah Hochgerner, Viktor Petukhov, Katja Lidschreiber, et al. "RNA Velocity of Single Cells" https://doi.org/10.1038/s41586-018-0414-6  Nature 560, no. 7719 (August 2018)
</details>


## Differential expression

- Differential expression methods developed for bulk RNA-seq work better in scRNA-seq in accounting for between-sample variability. 18 gold standard scRNA-seq datasets with matched bulk RNA-seq. Detailed data description in Methods. 14 DE methods tested, area under the concordance curve. edgeR DESeq2, limma perform best. Wilcoxon and MAST - medium to low performance. [GitHub](https://github.com/neurorestore/DE-analysis). <details>
    <summary>Paper</summary>
    Squair, Jordan W., Matthieu Gautier, Claudia Kathe, Mark A. Anderson, Nicholas D. James, Thomas H. Hutson, Rémi Hudelle, et al. “Confronting False Discoveries in Single-Cell Differential Expression.” Nature Communications 12, no. 1 (September 28, 2021): 5692. https://doi.org/10.1038/s41467-021-25960-2.
</details>

- Comparison of 11 differential gene expression detection methods for scRNA-seq data. Variable performance, poor overlap. Brief description of the statistics of each method. Bulk RNA-sec methods perform well, edgeR is good and fast. [Table 1. Software tools for identifying DE genes using scRNA-seq data](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2599-6/tables/1). <details>
    <summary>Paper</summary>
    Wang, Tianyu, Boyang Li, Craig E. Nelson, and Sheida Nabavi. "Comparative Analysis of Differential Gene Expression Analysis Tools for Single-Cell RNA Sequencing Data" https://doi.org/10.1186/s12859-019-2599-6  BMC Bioinformatics 20, no. 1 (December 2019)
</details>

- [Comparison of 36 differential detection methods in scRNA-seq data](http://imlspenticton.uzh.ch:3838/conquer/). Prefiltering of low expressed genes is important. edgeRQLFDetRate performs best. The conquer database, scRNA-seq datasets as RDS objects. <details>
    <summary>Paper</summary>
    Soneson, Charlotte, and Mark D Robinson. "Bias, Robustness and Scalability in Single-Cell Differential Expression Analysis" https://doi.org/10.1038/nmeth.4612  Nature Methods, February 26, 2018. 
</details>

- [waddR](https://bioconductor.org/packages/waddR/) - R package for testind distributional differences, using 2-Wasserstein distance, decomposed into change in location, shape, size, proportion of zeros (CMH test). Methods for one/multiple replicates per condition. Compared with scDD, SigEMD. Fast. In real data, detects similar genes to edgeR plus additional revealing new biology. [GitHub](https://github.com/goncalves-lab/waddR). <details>
    <summary>Paper</summary>
    Schefzik, Roman, Julian Flesch, and Angela Goncalves. "Fast Identification of Differential Distributions in Single-Cell RNA-Sequencing Data with WaddR" https://doi.org/10.1093/bioinformatics/btab226  Bioinformatics, 01 April 2021
</details>

- [singleCellHaystack](https://cran.r-project.org/web/packages/singleCellHaystack/index.html) - prediction of differentially expressed genes in scRNA-seq data in a multi-dimensional space, without explicit clustering. Grid in a multi-dimensional space, estimation of a reference distribution, comparing gene distribution at each grid point vs. reference using Kullback-Leibler divergence. Permutation to estimate significance. Compared with DEsingle, EMDo,ocs. scDD, edgeR, Monocle2, Seurat on synthetic (Splatter) and experimental scRNA-seq data (Tabula Muris), including spatial transcriptomics. Input - binarized expression matrix. Very fast. An R package, [GitHub](https://github.com/alexisvdb/singleCellHaystack). <details>
    <summary>Paper</summary>
    Vandenbon, Alexis, and Diego Diez. "A Clustering-Independent Method for Finding Differentially Expressed Genes in Single-Cell Transcriptome Data" https://doi.org/10.1038/s41467-020-17900-3  Nature Communications 11, no. 1 (2020): 1–10.
</details>

- [DEsingle](https://bioconductor.org/packages/DEsingle/) - an R package for detecting three types of differentially expressed genes in scRNA-seq data. Using Zero Inflated Negative Binomial distribution to distinguish true zeros from dropouts. Differential expression status (difference in the proportion of real zeros), DE abundance (typical DE without proportion of zeros change), DE general (both DE and proportion of zeros change). The majority of information is in [supplementary material](https://academic.oup.com/bioinformatics/article/34/18/3223/4983067#supplementary-data). <details>
    <summary>Paper</summary>
    Miao, Zhun, Ke Deng, Xiaowo Wang, and Xuegong Zhang. "DEsingle for Detecting Three Types of Differential Expression in Single-Cell RNA-Seq Data" https://doi.org/10.1093/bioinformatics/bty332  Edited by Bonnie Berger. Bioinformatics 34, no. 18 (September 15, 2018)
</details>

- [scDD](https://github.com/kdkorthauer/scDD) R package to identify differentially expressed genes in single cell RNA-seq data. Accounts for unobserved data. Four types of differential expression (DE, DP, DM, DB, see paper). <details>
    <summary>Paper</summary>
    Korthauer, Keegan D., Li-Fang Chu, Michael A. Newton, Yuan Li, James Thomson, Ron Stewart, and Christina Kendziorski. "A Statistical Approach for Identifying Differential Distributions in Single-Cell RNA-Seq Experiments" https://doi.org/10.1186/s13059-016-1077-y  Genome Biology 17, no. 1 (December 2016). 
</details>

- [MAST](https://github.com/RGLab/MAST) - scRNA-seq DEG analysis. CDR - the fraction of genes that are detectably expressed in each cell - added to the hurdle model that explicitly parameterizes distributions of expressed and non-expressed genes. Generalized linear model, log2(TPM+1), Gaussian. Regression coeffs are estimated using Bayesian approach. Variance shrinkage, gamma distribution. <details>
    <summary>Paper</summary>
    Finak, Greg, Andrew McDavid, Masanao Yajima, Jingyuan Deng, Vivian Gersuk, Alex K. Shalek, Chloe K. Slichter, et al. "MAST: A Flexible Statistical Framework for Assessing Transcriptional Changes and Characterizing Heterogeneity in Single-Cell RNA Sequencing Data" https://doi.org/10.1186/s13059-015-0844-5  Genome Biology 16 (December 10, 2015)
</details>

- [SCDE](https://hms-dbmi.github.io/scde/index.html) - Stochasticity of gene expression, high drop-out rate. A mixture model of two processes - detected expression and drop-out failure modeled as low-magnitude Poisson. Drop-out rate depends on the expected expression and can be approximated by logistic regression. <details>
    <summary>Paper</summary>
    Kharchenko, Peter V., Lev Silberstein, and David T. Scadden. "Bayesian Approach to Single-Cell Differential Expression Analysis" https://doi.org/10.1038/nmeth.2967  Nature Methods 11, no. 7 (July 2014)
</details>

### Differential abundance

- [Milo](https://github.com/MarioniLab/miloR) - an R package for differential abundance testing on scRNA-seq data between two groups or multiple conditions. Building a graph on the first 40 components of PCA, defining neighborhoods using a graph sampling algorithm. Each neighborhood (partially overlapping, in contrast to discrete clustering) contains cells from different conditions - differential abundance is tested using a negative binomial GLM. Tested on simulated datasets (dyntoy), a time course of mouse thymic epithelial cells development, liver cirrhosis analysis. Replicated datasets needed, batch corrected. Competitors: [DA-seq](https://github.com/KlugerLab/DAseq), [Cydar](https://bioconductor.org/packages/cydar/). [Code to reproduce results for the paper](https://github.com/MarioniLab/milo_analysis_2020). <details>
    <summary>Paper</summary>
    Dann, Emma, Neil C. Henderson, Sarah A. Teichmann, Michael D. Morgan, and John C. Marioni. "Milo: Differential Abundance Testing on Single-Cell Data Using k-NN Graphs" https://doi.org/10.1101/2020.11.23.393769  BioRxiv, January 1, 2020
</details>

## Downstream analysis

- [irGSEA](https://github.com/chuiqin/irGSEA) - R package for GSEA analysis of scRNA-seq data, differentially expressed genes between clusters. Input - Seurat object.

- [Hotspot](https://github.com/YosefLab/Hotspot) - a Python package for gene module identification in scRNA-seq data. Graph-based procedure, informative modules from co-variable genes in groups of cells (biological proximity, local autocorrelation within a KNN similarity graph, can be defined using multiple metrics, including spatial). Three steps: compute a similarity map, select informative genes exhibiting non-random variation within the similarity map, group these genes into modules. Benchmarked against WGCRNA and GRNBoost. <details>
    <summary>Paper</summary>
    DeTomaso, David, and Nir Yosef. “Hotspot Identifies Informative Gene Modules across Modalities of Single-Cell Genomics.” Cell Systems 12, no. 5 (May 2021): 446-456.e9. https://doi.org/10.1016/j.cels.2021.04.005.
</details>



## CNV

- [SCEVAN](https://github.com/AntonioDeFalco/SCEVAN) - R package that automatically classifies the cells in the scRNA data by segregating non-malignant cells of tumor microenviroment from the malignant cells. It also infers the copy number profile of malignant cells, identifies subclonal structures and analyses the specific and shared alterations of each subpopulation. <details>
    <summary>Paper</summary>
    De Falco, Antonio, Francesca Caruso, Xiao-Dong Su, Antonio Iavarone, and Michele Ceccarelli. “A Variational Algorithm to Detect the Clonal Copy Number Substructure of Tumors from scRNA-Seq Data.” Nature Communications 14, no. 1 (February 25, 2023): 1074. https://doi.org/10.1038/s41467-023-36790-9.
</details>

- [Numbat](https://github.com/kharchenkolab/numbat) - Numbat is a haplotype-aware CNV caller from single-cell and spatial transcriptomics data. It integrates signals from gene expression, allelic ratio, and population-derived haplotype information to accurately infer allele-specific CNVs in single cells and reconstruct their lineage relationship. <details>
    <summary>Paper</summary>
    Teng Gao, Ruslan Soldatov, Hirak Sarkar, Adam Kurkiewicz, Evan Biederstedt, Po-Ru Loh, Peter Kharchenko. Haplotype-aware analysis of somatic copy number variations from single-cell transcriptomes. Nature Biotechnology (2022). https://doi.org/10.1038/s41587-022-01468-y
</details>

- [CopyKAT](https://github.com/navinlabcode/copykat) - CopyKAT (Copynumber Karyotyping of Tumors) is a computational tool using integrative Bayesian approaches to identify genome-wide aneuploidy at 5MB resolution in single cells to separate tumor cells from normal cells, and tumor subclones using high-throughput sc-RNAseq data.  <details>
    <summary>Paper</summary>
    Ruli Gao, Shanshan Bai, Ying C. Henderson, Yiyun Lin, Aislyn Schalck, Yun Yan, Tapsi Kumar, Min Hu, Emi Sei, Alexander Davis, Fang Wang, Simona F. Shaitelman, Jennifer Rui Wang, Ken Chen, Stacy Moulder, Stephen Y. Lai and Nicholas E. Navin (2021). Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes. Nature Biotechnology (2021). https://doi.org/10.1038/s41587-020-00795-2
</details>

- [CaSpER](https://github.com/akdess/CaSpER) - identification of CNVs from  RNA-seq data, bulk and single-cell (full-transcript only, like SMART-seq). Utilized multi-scale smoothed global gene expression profile and B-allele frequency (BAF) signal profile, detects concordant shifts in signal using a 5-state HMM (homozygous deletion, heterozygous deletion, neutral, one-copy-amplification, high-copy-amplification). Reconstructs subclonal CNV architecture for scRNA-seq data. Tested on GBM scRNA-seq, TCGA, other. Compared with HoneyBADGER. R code and tutorials. <details>
    <summary>Paper</summary>
    Serin Harmanci, Akdes, Arif O. Harmanci, and Xiaobo Zhou. "CaSpER Identifies and Visualizes CNV Events by Integrative Analysis of Single-Cell or Bulk RNA-Sequencing Data" https://doi.org/10.1038/s41467-019-13779-x  Nature Communications 11, no. 1 (December 2020)
</details>

- [chromoscope](https://github.com/hms-dbmi/chromoscope) - multi-omics circular visualization of cancer genomes. A browser-based application (Python, Jupyter notebook support). Multi-sample view, genome view, variant view, breakpoint view. Implemented using the Gosling library, Python package available. Input: BEDPE (structural variants), TSV (CNVs, driver mutations), VCF, BAM. [Web demo](https://chromoscope.bio/), [Documentation](https://chromoscope.bio/docs/#/). <details>
    <summary>Paper</summary>
    L’Yi, Sehi, Dominika Maziec, Victoria Stevens, Trevor Manz, Alexander Veit, Michele Berselli, Peter J Park, Dominik Glodzik, and Nils Gehlenborg. “Chromoscope: Interactive Multiscale Visualization for Structural Variation in Human Genomes.” Preprint. Open Science Framework, May 10, 2023. https://doi.org/10.31219/osf.io/pyqrx.
</details>

- CNV estimation algorithm in scRNA-seq data - moving 100-gene window, deviation of expression from the chromosome average. Details in Methods. <details>
    <summary>Paper</summary>
    Tirosh, Itay, Andrew S. Venteicher, Christine Hebert, Leah E. Escalante, Anoop P. Patel, Keren Yizhak, Jonathan M. Fisher, et al. "Single-Cell RNA-Seq Supports a Developmental Hierarchy in Human Oligodendroglioma" https://doi.org/10.1038/nature20123  Nature 539, no. 7628 (November 2016)
</details>

- [infercnv](https://github.com/broadinstitute/infercnv) - Inferring copy number alterations from tumor single cell RNA-Seq data, as compared with a set of reference normal cells. Positional expression intensity comparison. [Documentation](https://github.com/broadinstitute/inferCNV/wiki).

## Splicing

- [SICILIAN](https://github.com/salzmanlab/SICILIAN) (SIngle Cell precIse spLice estImAtioN) - splice junction from bulk and scRNA-seq data. Deconvolves biochemical and computational noise, uses generalized linear modeling with various read mapping features as predictors. Input: a BAM file from a splice-aware aligner (e.g., STAR). Three steps: (1) assign a statistical score to each junctional read’s alignment to quantify the likelihood that the read alignment is truly from RNA expression; (2) aggregate read scores to summarize the likelihood that a given junction is a true positive; and (3) report single-cell resolved junction expression quantification, corrected for multiple hypotheses testingHigher accuracy on simulated data, improves concordance between matched scRNA-seq and bulk datasets, biological replicates. Python and R implementation. <details>
    <summary>Paper</summary>
    Dehghannasiri, Roozbeh, Julia Eve Olivieri, and Julia Salzman. “Specific Splice Junction Detection in Single Cells with SICILIAN.” Preprint. Bioinformatics, April 15, 2020. https://doi.org/10.1101/2020.04.14.041905.
</details>

## Annotation, subpopulation identification

- [HiCAT](https://github.com/combio-dku/HiCAT) - hierarchical, marker-based cell-type identifier using gene set analysis for statistical marker scoring. Major-minor-subtype classification as a three-level taxonomy tree. Best for immune-cell annotations. Six steps: 1. PCA and clustering, 2. marker counting and gene set analysis scoring, 3. unknown cluster detection, 4. Gaussian mixture model-based correction, 5. rejecting unclear cells, 6. kNN-based correction. Tested on eight normal and cancer single-cell datasets, including BRCA. Input: cell x gene matrix and taxonomy tree & subset markers. [Cell marker database](https://www.rndsystems.com/resources/cell-markers). Python, [R wrapper](https://github.com/combio-dku/HiCAT). <details>
    <summary>Paper</summary>
    Lee, Joongho, Minsoo Kim, Keunsoo Kang, Chul-Su Yang, and Seokhyun Yoon. “Hierarchical Cell-Type Identifier Accurately Distinguishes Immune-Cell Subtypes Enabling Precise Profiling of Tissue Microenvironment with Single-Cell RNA-Sequencing,” n.d.
</details>

- [SCimilarity](https://github.com/Genentech/scimilarity) - cell type annotation based on SCimilarity score. A deep learning model using one fully connected encoder (reduced dimensionality to 128) and one decoder, reusing the network three times per training triplet, heavy dropout regularization. Optimizing two objectives: a supervised triplet loss function used for embedding and an unsupervised mean squared error loss function for presenving variation across gene expression profiles. Trained on 22.7M cells from 399 10x Genomics scRNA-seq studies. Uses Cell Ontology terms, only normal cells (cancer cells and cell lines excluded). Can be used for data integration, outperforms Harmony, scVI, scArches. Query cells can be matched to the full spectrum of cells, or a tissue-specific subset. [Detailed documentatio and tutorial](https://genentech.github.io/scimilarity/index.html), using pretrained model does not require GPU. <details>
    <summary>Paper</summary>
    Heimberg, Graham, Tony Kuo, Daryle DePianto, Tobias Heigl, Nathaniel Diamant, Omar Salem, Gabriele Scalia, et al. “Scalable Querying of Human Cell Atlases via a Foundational Model Reveals Commonalities across Fibrosis-Associated Macrophages.” Preprint. Bioinformatics, July 19, 2023. https://doi.org/10.1101/2023.07.18.549537.
</details>

- [CellTypist](https://www.celltypist.org/) - machine learning tool for precise cell type annotation, immune cell types. Trained on 20 tissues with harmonized cell type labels, hierarchy of 45 cell types. L2-regularized logistic regression, machine learning framework, gradient descent, 30 epoch. Scanpy pipeline, batch correction using [bbknn](https://github.com/Teichlab/bbknn), markers detection using [rbcde](https://github.com/Teichlab/rbcde). [GitHub](https://github.com/Teichlab/celltypist). <details>
    <summary>Paper</summary>
    Domínguez Conde, C, C Xu, Lb Jarvis, T Gomes, Sk Howlett, Db Rainbow, O Suchanek, et al. "Cross-Tissue Immune Cell Analysis Reveals Tissue-Specific Adaptations and Clonal Architecture in Humans" https://doi.org/10.1101/2021.04.28.441762  Preprint. Immunology, April 28, 2021. 
</details>

- [tricycle](https://github.com/hansenlab/tricycle) - transfer learning approach to learn cell cycle PCA projections from a reference dataset and project new data on it. Combining the biology of the cell cycle, the mathematical properties of PCA of unimodal periodicity of genes associated with cell cycle. [Tweet](https://twitter.com/loyalgoff/status/1379806742167527427?s=20). <details>
    <summary>Paper</summary>
    Zheng, Shijie C., Genevieve Stein-O’Brien, Jonathan J. Augustin, Jared Slosberg, Giovanni A. Carosso, Briana Winer, Gloria Shin, Hans T. Bjornsson, Loyal A. Goff, and Kasper D. Hansen. "[Universal prediction of cell cycle position using transfer learning" https://doi.org/10.1101/2021.04.06.438463   bioRxiv April 11, 2021
</details>

- [igrabski/scRNAseq-cell-type](https://github.com/igrabski/scRNAseq-cell-type) - A statistical approach for cell type annotation from scRNA-seq data. Considers all genes, uses a latent variable model to define cell-type-specific barcodes and probabilistically annotates cell type identity, while accounting for batch effects. Methods, modeling gene-specific distribution using off-low, off-high, on states from scRNA expression bimodal distribution. Train the model using reference data from the PanglaoDB https://panglaodb.se/, tested the method on the PBMC, colon, and brain scRNA-seq datasets. Clustering methods tend to overcluster, marker genes are unreliable due to sparsity. Outperforms scmap, CaSTLe, SingleR, Garnett, CellAssign. [R code](https://github.com/igrabski/scRNAseq-cell-type), [Tweet](https://twitter.com/rafalab/status/1341771349254135808?s=20). <details>
    <summary>Paper</summary>
    Isabella N. Grabski, Rafael A. Irizarry "A Probabilistic Gene Expression Barcode for Annotation of Cell-Types from Single Cell RNA-Seq Data" https://doi.org/10.1101/2020.01.05.895441  biorXiv, December 22, 2020
</details>

- [scTyper](https://github.com/omicsCore/scTyper) - cell- or cluster typing analysis of scRNA-seq data. Includes scTyper.db and CellMarker databases with literature-curated 213 cell marker sets (malignant cells, cancer-associated fibroblasts, 149 immune cells). Three enrichment methods (Nearest template prediction NTP, gene set enrichment analysis GSEA, average expression values). DNA copy number inference (inferCNV, removing gene clusters) for malignant cell typing. Supports 10X Genomics CellRanges and the Seurat package using the corresponding modules: “QC”, “Cell Ranger”, “Seurat processing”, “cell typing”, and “malignant cell typing”. <details>
    <summary>Paper</summary>
    Choi, Ji-Hye, Hye In Kim, and Hyun Goo Woo. “ScTyper: A Comprehensive Pipeline for the Cell Typing Analysis of Single-Cell RNA-Seq Data.” BMC Bioinformatics 21, no. 1 (December 2020): 342. https://doi.org/10.1186/s12859-020-03700-5.
</details>

- [Azimuth](https://azimuth.hubmapconsortium.org/) - Mapping query scRNA-seq dataset to multimodal references and assigning cell types. Supervised principal component analysis to identify a projection of the query dataset that maximally captures the structure defined by the WNN graph. Combined with the anchor-based framework, allows projection on the previously defined reference UMAP visualization. Human PBMC, motor cortex, pancreas, mouse motor cortex references. [Online apps](https://azimuth.hubmapconsortium.org/). <details>
    <summary>Paper</summary>
    Hao, Yuhan, Stephanie Hao, Erica Andersen-Nissen, William M. Mauck, Shiwei Zheng, Andrew Butler, Maddie Jane Lee, et al. "Integrated Analysis of Multimodal Single-Cell Data" https://doi.org/10.1101/2020.10.12.335331  Preprint. Genomics, October 12, 2020. 
</details>

- [MARS](https://github.com/snap-stanford/mars) - a meta-learning approach for identifying known and new cell types in scRNA-seq data. Constructs a meta-dataset from experiments with annotated cell types (used to learn the cell type landmarks in the embedding space) and an unannotated experiment (mathed to the embedded landmarks). The embedding space and objective function are defined such that cells (annotated and unannotated) embed close to their cell-type landmarks, cell type landmarks are most distinct. Autoencoder with 1000 and 100 neurons, input - all 22.9K genes. Applied to the Tabula Muris Senis dataset, several others. Significantly outperform ScVi, SIMLR, Scanpy and Seurat on adjusted Rand index, adjusted MI and other metrics. <details>
    <summary>Paper</summary>
    Brbić, Maria, Marinka Zitnik, Sheng Wang, Angela O. Pisco, Russ B. Altman, Spyros Darmanis, and Jure Leskovec. "MARS: Discovering Novel Cell Types across Heterogeneous Single-Cell Experiments" https://doi.org/10.1038/s41592-020-00979-3  Nature Methods, October 19, 2020.
</details>

- [Garnett](https://cole-trapnell-lab.github.io/garnett/) - annotating cells in scRNA-seq data. Hierarchy of cell types and their markers should be pre-defined using a markup language. A classifier is trained to classify additional datasets. Trained on cells from one organisms, can be applied to different organisms. Pre-trained classifiers available. R-based. <details>
    <summary>Paper</summary>
    Pliner, Hannah A., Jay Shendure, and Cole Trapnell. "Supervised Classification Enables Rapid Annotation of Cell Atlases" https://doi.org/10.1038/s41592-019-0535-3  Nature Methods 16, no. 10 (October 2019)
</details>

- [COMET](https://github.com/MSingerLab/COMETSC) - a computational framework (Python 3.6) for the identification of candidate marker panels consisting of one or more genes for cell populations of interest identified with scRNA-seq data. Uses the XL-minimal HyperGeometric test (XL-mHG test) to binarize each gene as being expressed in a cluster-specific manner. Input: a gene-by-cell expression matrix, cluster assignments, and visualization coordinates. Output: gene markers most positively and negatively predictive of a cluster membership. [Documentation](https://hgmd.readthedocs.io/en/latest/). <details>
    <summary>Paper</summary>
    Delaney, Conor, Alexandra Schnell, Louis V Cammarata, Aaron Yao‐Smith, Aviv Regev, Vijay K Kuchroo, and Meromit Singer. “Combinatorial Prediction of Marker Panels from Single‐cell Transcriptomic Data.” Molecular Systems Biology 15, no. 10 (October 2019). https://doi.org/10.15252/msb.20199005.
</details>

- [Single-Cell Signature Explorer](https://sites.google.com/site/fredsoftwares/products/single-cell-signature-explorer) - scoring (sum of UMIs in a signature over total UMIs in a cell) of gene-set signatures at the single cell level and their visualization using t-SNE or UMAP. Contains four successive tools: Scorer (computes signature scores for each cell), Merger (collates score table with t-SNE and UMAP coordinates), Viewer (visualizes signature scores, one at a time), Combiner (arithmetic combination of two signature scores, and visualization). Utilizes sctransform normalization. Implemented in Go, with viewers implemented in R/Shiny. Outperforms ssGSEA, GSVA, AutoCompare_SES. [Signature databases](https://sites.google.com/site/fredsoftwares/products/databases) for multiple databases. <details>
    <summary>Paper</summary>
    Pont, Frédéric, Marie Tosolini, and Jean J Fournié. “Single-Cell Signature Explorer for Comprehensive Visualization of Single Cell Signatures across ScRNA-Seq Datasets.” Nucleic Acids Research 47, no. 21 (December 2, 2019): e133–e133. https://doi.org/10.1093/nar/gkz601.
</details>

- [CellAssign](https://github.com/irrationone/cellassign) - R package for scRNA-seq cell type inference. Probabilistic graphical model to assign cell type probabilities to single cells using known marker genes (binarized matrix), including "unassigned" categorization. Insensitive to batch- or sample-specific effects. Outperforms Seurat, SC3, PhenoGraph, densityCut, dynamicTreeCut, scmap-cluster, correlation-based methods, SCINA. Applied to delineate the composition of the tumor microenvironment. Built using TensorFlow. <details>
    <summary>Paper</summary>
    Zhang, Allen W., Ciara O’Flanagan, Elizabeth A. Chavez, Jamie L. P. Lim, Nicholas Ceglia, Andrew McPherson, Matt Wiens, et al. "Probabilistic Cell-Type Assignment of Single-Cell RNA-Seq for Tumor Microenvironment Profiling" https://doi.org/10.1038/s41592-019-0529-1  Nature Methods, September 9, 2019. 
</details>

- [SCINA](https://lce.biohpc.swmed.edu/scina/) - semi-supervised cell annotation algorithm for scRNA-seq data. Uses previously established gene signatures (one or more genes) and an expectation-maximization (EM) algorithm. Assumes bimodal distribution for each signature gene, with higher mode corresponding to the cell type of interest. Cells with overall low signature expression (all signatures) are labeled unknown. Simulated and experimental signatures/data, outperforms K-means clustering, Seurat, SINCERA, PhenoGraph, robust to noise. [GitHub](https://github.com/jcao89757/SCINA), [R package](https://CRAN.R-project.org/package=SCINA) and a [web-server](https://lce.biohpc.swmed.edu/scina/analysis.php). <details>
    <summary>Paper</summary>
    Zhang, Luo, Zhong, Choi, Ma, Wang, Mahrt, et al. “SCINA: Semi-Supervised Analysis of Single Cells in Silico.” Genes 10, no. 7 (July 12, 2019): 531. https://doi.org/10.3390/genes10070531.
</details>

- [SingleCellNet](http://github.com/pcahan1/singleCellNet/) - quantitative cell type annotation. Top-scoring pair transformation to match query and reference datasets. Compared with SCMAP, binary cell type classifier based on correlation. Benchmarked on 12 scRNA-seq datasets, provided in the GitHub repo, http://github.com/pcahan1/singleCellNet/. [Blog post](https://www.rna-seqblog.com/singlecellnet-a-computational-tool-to-classify-single-cell-rna-seq-data-across-platforms-and-across-species/). <details>
    <summary>Paper</summary>
    Tan, Yuqi, and Patrick Cahan. "SingleCellNet: A Computational Tool to Classify Single Cell RNA-Seq Data Across Platforms and Across Species" https://doi.org/10.1016/j.cels.2019.06.004  Cell Systems, July 2019. 
</details>

- [scPopCorn](https://github.com/ncbi/scPopCorn/) - subpopulation identification across scRNA-seq experiments. Identifies shared and unique subpopulations. Joint network of two graphs. First, graphs are built for each experiment using co-expression to identify subpopulations. Second, the corresponsence of the identified subpopulations is refined using Google's PageRank algorithm to identify subpopulations. Compared with Seurat alignment + Louvain, mutual nearest neighbor (MNN) method, and MNN + Louvain. Several assessment metrics. Tested on pancreatic, kidney cells, healthy brain and glioblastoma scRNA-seq data. Sankey diagrams showing how subpopulation assignment change. <details>
    <summary>Paper</summary>
    Wang, Yijie, Jan Hoinka, and Teresa M. Przytycka. "Subpopulation Detection and Their Comparative Analysis across Single-Cell Experiments with ScPopCorn" https://doi.org/10.1016/j.cels.2019.05.007  Cell Systems 8, no. 6 (June 2019)
</details>

- [matchSCore2](https://github.com/elimereu/matchSCore2) - classifying cell types based on reference data. <details>
    <summary>Paper</summary>
    Mereu, Elisabetta, Atefeh Lafzi, Catia Moutinho, Christoph Ziegenhain, Davis J. MacCarthy, Adrian Alvarez, Eduard Batlle, et al. "Benchmarking Single-Cell RNA Sequencing Protocols for Cell Atlas Projects" https://doi.org/10.1101/630087  Preprint. Genomics, May 13, 2019. 
</details>

- [Single-Cell Signature Explorer](https://sites.google.com/site/fredsoftwares/products/single-cell-signature-explorer) - gene signature (\~17,000 from MSigDb, KEGG, Reactome) scoring (sum of UMIs in in a gene signature over the total UMIs in a cell) for single cells, and visualization on top of a t-SNE plot. Optional Noise Reduction (Freeman-Tuckey transform to stabilize technical noise). Four consecutive tools (Go language, R/Shiny). Comparison with Seurat's Cell CycleScore module and AUCell from SCENIC. Very fast. <details>
    <summary>Paper</summary>
    Pont, Frédéric, Marie Tosolini, and Jean Jacques Fournié. "Single-Cell Signature Explorer for Comprehensive Visualization of Single Cell Signatures across ScRNA-Seq Data Sets" https://doi.org/10.1101/621805  Preprint. Bioinformatics, April 29, 2019. 
</details>

- [SingleR](https://github.com/dviraran/SingleR/) - scRNA-seq cell type assignment by (Spearman) correlating to reference bulk RNA-seq data of pure cell types. Validated on ImmGen data. The package provides Human Primary Cell Atlas data, Blueprint and ENCODE consortium data, ImmGen, three others as a reference. Post-Seurat analysis. [Web tool](http://comphealth.ucsf.edu/SingleR/), that takes SingleR objects, instructions are on GitHub, https://github.com/dviraran/SingleR/. [Example analysis](http://comphealth.ucsf.edu/sample-apps/SingleR/SingleR.MCA.html), [Bioconductor package](https://bioconductor.org/packages/devel/bioc/html/SingleR.html), [Twitter](https://twitter.com/dvir_a/status/1170117086711930880?s=03). <details>
    <summary>Paper</summary>
    Aran, Dvir, Agnieszka P. Looney, Leqian Liu, Esther Wu, Valerie Fong, Austin Hsu, Suzanna Chak, et al. "Reference-Based Analysis of Lung Single-Cell Sequencing Reveals a Transitional Profibrotic Macrophage" https://doi.org/10.1038/s41590-018-0276-y  Nature Immunology 20, no. 2 (February 2019)
</details>

- [TooManyCells](https://github.com/GregorySchwartz/tooManyCellsR) - divisive hierarchical spectral clustering of scRNA-seq data. Uses truncated singular vector decomposition to bipartition the cells. Newman-Girvain modularity Q to assess whether bipartition is significant or should be stopped. BirchBeer visualization. Outperforms Phenograph, Seurat, Cellranger, Monocle, the latter is second in performance. Excels for rare populations. Normalization marginally affects performance. <details>
    <summary>Paper</summary>
    Schwartz, Gregory W, Jelena Petrovic, Maria Fasolino, Yeqiao Zhou, Stanley Cai, Lanwei Xu, Warren S Pear, Golnaz Vahedi, and Robert B Faryabi. "TooManyCells Identifies and Visualizes Relationships of Single-Cell Clades" https://doi.org/10.1101/519660  BioRxiv, January 13, 2019. 
</details>

- [VISION](https://github.com/YosefLab/VISION) - functional annotation of scRNA-seq data using gene signatures (Geary's C statistics), unsupervised and supervised. Operates downstream of dimensionality reduction, clustering. A continuation of FastProject. <details>
    <summary>Paper</summary>
    DeTomaso, David, Matthew Jones, Meena Subramaniam, Tal Ashuach, Chun J Ye, and Nir Yosef. "Functional Interpretation of Single-Cell Similarity Maps" https://doi.org/10.1101/403055  August 29, 2018. 
</details>

### Cell markers

- Benchmarking study of 41 computational methods for marker gene selection in scRNA-seq data. 10 experimental and over 170 simulated datasets (splatter R package). CellMarker and PangaloDB marker gene set databases, [Supplementary tables S1-S4](https://www.biorxiv.org/content/10.1101/2022.05.09.490241v1.supplementary-material) - expert-annotated marker gene sets pbmc3k immune cells, Lawlor pancreatic cells, Zeisel brain cells, Smart-seq3 immune cells, [CSV](data/Pullin_McCarthy_2022.csv). Poor overlap among methods. Wilcoxon test performs best, Student t-test (not Welch) and logistic resression also OK. [Scripts](https://gitlab.svi.edu.au/biocellgen-public/mage_2020_marker-gene-benchmarking/) wrapped in a Snakemake pipeline. <details>
    <summary>Paper</summary>
    Pullin, Jeffrey M., and Davis J. McCarthy. "A comparison of marker gene selection methods for single-cell RNA sequencing data." bioRxiv (May 10, 2022). https://doi.org/10.1101/2022.05.09.490241
</details>

- [UCell](https://github.com/carmonalab/UCell) - an R package for gene signature enrichment in scRNA-seq data based on the Mann-Whitney U statistics. Integrates with the Seurat pipeline. Annotates each cell with signature enrichments. Requires gene lists as signature definitions. <details>
    <summary>Paper</summary>
    Andreatta, Massimo, and Santiago J. Carmona. “UCell: Robust and Scalable Single-Cell Gene Signature Scoring.” Computational and Structural Biotechnology Journal 19 (2021): 3796–98. https://doi.org/10.1016/j.csbj.2021.06.043.
</details>

- [clustermole](https://github.com/igordot/clustermole) - blindly digging for cell types in scRNA-seq clusters. Cell type prediction based on marker genes, cell type prediction based on a full expression matrix, a database of cell type markers. https://github.com/igordot/clustermole, [CRAN](https://cran.r-project.org/web/packages/clustermole/index.html)

- [scMatch](https://github.com/asrhou/scMatch) - Python tool for annotating scRNA-seq cells by their closest match (Spearman, Pearson correlation) in large reference datasets (FANTOM5, SingleR, Xena Cancer browser). <details>
    <summary>Paper</summary>
    Hou, Rui, Elena Denisenko, and Alistair R. R. Forrest. "ScMatch: A Single-Cell Gene Expression Profile Annotation Tool Using Reference Datasets" https://doi.org/10.1093/bioinformatics/btz292  Bioinformatics (Oxford, England), April 26, 2019. 
</details>

- [scGeneFit](https://github.com/solevillar/scGeneFit) - selection of hierarchical gene markers, contrasted with one-vs-all gene selection. MATLAB implementation. <details>
    <summary>Paper</summary>
    Dumitrascu, Bianca, Soledad Villar, Dustin G. Mixon, and Barbara E. Engelhardt. "Optimal Gene Selection for Cell Type Discrimination in Single Cell Analyses" https://doi.org/10.1101/599654  BioRxiv, April 4, 2019. 
</details>

- [CancerSEA](http://biocc.hrbmu.edu.cn/CancerSEA/home.jsp) - cancer scRNA-seq studies. [Download](http://biocc.hrbmu.edu.cn/CancerSEA/goDownload) individual studies, as well as gene signatures of 14 states (fincluding stemness, invasion, metastasis, proliferation, EMT, angiogenesis, apoptosis, cell cycle, differentiation, DNA damage, DNA repair, hypoxia, inflammation and quiescence). ssGSEA on each cell. [Supplementary Table S3](https://academic.oup.com/nar/article/47/D1/D900/5133662#supplementary-data) - metastasis-associated genes across all cancers. <details>
    <summary>Paper</summary>
    Yuan, Huating, Min Yan, Guanxiong Zhang, Wei Liu, Chunyu Deng, Gaoming Liao, Liwen Xu et al. "CancerSEA: a cancer single-cell state atlas." Nucleic acids research 47, no. D1 (08 January 2019) https://doi.org/10.1093/nar/gky939
</details>

- [CellMarker](http://xteam.xbio.top/CellMarker/) - cell type-specific markers, human, mouse, single-cell data. 13,605 cell markers of 467 cell types in 158 human tissues/sub-tissues and 9,148 cell makers of 389 cell types in 81 mouse tissues/sub-tissues. [Browse](http://xteam.xbio.top/CellMarker/browse.jsp) signatures for various cell/tissue types. [Search](http://xteam.xbio.top/CellMarker/search.jsp) individual markers. Tab-separated [download](http://xteam.xbio.top/CellMarker/download.jsp) of human, mouse, and single-cell markers. [Ehrichment analysis using CellMarker's data](https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html). <details>
    <summary>Paper</summary>
    Zhang, Xinxin, Yujia Lan, Jinyuan Xu, Fei Quan, Erjie Zhao, Chunyu Deng, Tao Luo, et al. “CellMarker: A Manually Curated Resource of Cell Markers in Human and Mouse.” Nucleic Acids Research 47, no. D1 (January 8, 2019): D721–28. https://doi.org/10.1093/nar/gky900.
</details>

#### Immune markers

See [Immuno_notes, immune-markers](https://github.com/mdozmorov/Immuno_notes#immune-markers)

#### Brain markers

- [Brain immune atlas](http://www.brainimmuneatlas.org/index.php) scRNA-seq resource. Border-associated macrophages from discrete mouse brain compartments, tissue-specific transcriptional signatures. [GitHub](https://github.com/saeyslab/brainimmuneatlas/). <details>
    <summary>Paper</summary>
    Van Hove, Hannah, Liesbet Martens, Isabelle Scheyltjens, Karen De Vlaminck, Ana Rita Pombo Antunes, Sofie De Prijck, Niels Vandamme, et al. “A Single-Cell Atlas of Mouse Brain Macrophages Reveals Unique Transcriptional Identities Shaped by Ontogeny and Tissue Environment.” Nature Neuroscience, May 6, 2019. https://doi.org/10.1038/s41593-019-0393-4.
</details>


## Phylogenetic inference

- [TAR-scRNA-seq](https://github.com/fw262/TAR-scRNA-seq) - snakemake pipeline for transcript annotation/prediction from scRNA-seq data of poorly annotated genomes. Transcriptionally Active Regions (TARs) detected using groHMM, a hidden Markov model to divide the genome into two states (transcribed and untranscribed, 50bp window). Unannotated TADs (uTARs) preserve cell-type-dependent clusters. Annotated using homology (BLASTn). Spatial transcriptomics confirms coexpression. Analysis of human, mouse, chicken, mole rat, lemur and sea urchin genomes (10X Genomics data). [Summary of data availability](https://www.nature.com/articles/s41467-021-22496-3#data-availability). <details>
    <summary>Paper</summary>
    Wang, Michael F. Z., Madhav Mantri, Shao-Pei Chou, Gaetano J. Scuderi, David W. McKellar, Jonathan T. Butcher, Charles G. Danko, and Iwijn De Vlaminck. “Uncovering Transcriptional Dark Matter via Gene Annotation Independent Single-Cell RNA Sequencing Analysis.” Nature Communications 12, no. 1 (December 2021): 2158. https://doi.org/10.1038/s41467-021-22496-3.
</details>

- [OncoNEM](https://bitbucket.org/edith_ross/onconem/src/master/) (oncogenetic nested effects model) - tumor evolution inference from single cell data from somatic SNPs of single cells. Identifies homogeneous subpopulations and infers their genotypes and phylogenetic tree. Probabilistically accounts for noise in the observed genotypes, allele dropouts, unobserved subpopulations. Input - binary genotype matrix, false positive and negative rates. Output - inferred tumor subpopulations, evolutionary tree, posterior probabilities of mutations. Assessed in simulation studies, outperforms similar methods. Robust to the selection of parameters. <details>
    <summary>Paper</summary>
    Ross, Edith M., and Florian Markowetz. "OncoNEM: Inferring Tumor Evolution from Single-Cell Sequencing Data" https://doi.org/10.1186/s13059-016-0929-9  Genome Biology, (December 2016)
</details>

## Immuno-analysis

- [TCR-BCR-seq-analysis](https://github.com/crazyhottommy/TCR-BCR-seq-analysis) - 
T/B cell receptor sequencing analysis notes by Ming Tang. Also, [23 tools to work with (single-cell) TCR/BCR-seq immune repertoire data](https://crazyhottommy.blogspot.com/2022/12/23-tools-to-work-with-single-cell.html) by Ming, [Tweet](https://twitter.com/tangming2005/status/1600868620128718848?s=20&t=EzFDzKyBvH1mykAomlahGg).

- Overview of B-cell receptor development/affinity maturation, sequencing B-cell repertoire, clonal lineage assignment and clustering, somatic hypermutation analysis, challenges. <details>
    <summary>Paper</summary>
    Hoehn, Kenneth B., Anna Fowler, Gerton Lunter, and Oliver G. Pybus. “The Diversity and Molecular Evolution of B-Cell Receptors during Infection.” Molecular Biology and Evolution 33, no. 5 (May 2016): 1147–57. https://doi.org/10.1093/molbev/msw015.
</details>

- [DALI](https://github.com/vibscc/DALI) - R-package for the analysis of single-cell TCR/BCR data and scRNA-seq (10X Genomics) in the Seurat ecosystem. Read10X_vdj reads Cellranger multi data, Interactive_VDJ launches Shiny app. Input - scRNA-seq Seurat object (.rds), and vdj data. [Demo data](https://cloud.irc.ugent.be/public/index.php/s/9ys5czsaNtNQtSd). <details>
    <summary>Paper</summary>
    Verstaen, Kevin, Inés Lammens, Jana Roels, Yvan Saeys, Bart N Lambrecht, Niels Vandamme, and Stijn Vanhee. “DALI (Diversity AnaLysis Interface): A Novel Tool for the Integrated Analysis of Multimodal Single Cell RNAseq Data and Immune Receptor Profiling.” Preprint. Bioinformatics, December 7, 2021. https://doi.org/10.1101/2021.12.07.471549.
</details>

- [TRUST4](https://github.com/liulab-dfci/TRUST4) - reconstructing immune repertoires of T- and B cells from bulk and scRNA-seq (10X Genomics 5'). Input: FASTQ or BAM. Faster, more sensitive (MiXCR, CATT, TRUST3). Uses the international ImmunoGeneTics (IMGT) database. Methods describe all algorithms. [Scripts to reproduce paper](https://github.com/liulab-dfci/TRUST4_manuscript_evaluation). <details>
    <summary>Paper</summary>
    Song, Li, David Cohen, Zhangyi Ouyang, Yang Cao, Xihao Hu, and X. Shirley Liu. “TRUST4: Immune Repertoire Reconstruction from Bulk and Single-Cell RNA-Seq Data.” Nature Methods 18, no. 6 (June 2021): 627–30. https://doi.org/10.1038/s41592-021-01142-2.
</details>

- [scRepertoire](https://github.com/ncborcherding/scRepertoire) - R package for 10x Genomics Chromium Immune Profiling for both T-cell receptor (TCR) and immunoglobulin (Ig) enrichment data. Clonotype analysis: cdr3 distribution, clonotype calling, proportion, repertoire overlap, diversity. Comparison between samples. Rich visualization capabilities. [Tweet](https://twitter.com/theHumanBorch/status/1218998979247120386?s=20). <details>
    <summary>Paper</summary>
    Borcherding, N, NL Bormann, and G Kraus. “ScRepertoire: An R-Based Toolkit for Single-Cell Immune Receptor Analysis [Version 2; Peer Review: 2 Approved].” F1000Research 9, no. 47 (2020). https://doi.org/10.12688/f1000research.22139.2.
</details>

- [enclone](https://10xgenomics.github.io/enclone/) - Accurate and user-friendly computational tool for clonal grouping to study the adaptive immune system. Analyzes 10x Genomics Chromium Single Cell V(D)J data containing B cell receptor (BCR) and T cell receptor (TCR) RNA sequences. [GitHub](https://github.com/10XGenomics/enclone)

- [immunarch](https://immunarch.com) - Exploration of Single-cell and Bulk T-cell/Antibody Immune Repertoires in R

- [Immcantation](https://immcantation.readthedocs.io) - a start-to-finish analytical ecosystem for high-throughput AIRR-seq (adaptive immune receptor repertoire) datasets. Beginning from raw reads, Python and R packages are provided for pre-processing, population structure determination, and repertoire analysis (pRESTO, Change-O, Alakazam, SHazaM, TIgGER, SCOPer, dowser, RDI, RAbHIT, IgPhyML, sumrep). Pipelines for various technologies, inclusing 10X Genomics. [Tutorials](https://immcantation.readthedocs.io/en/stable/tutorials/tutorials.html). [Introduction to B cell repertoire analysis using the Immcantation framework](https://kleinstein.bitbucket.io/tutorials/intro-lab/index.html) and the [Jupyter notebook](https://bitbucket.org/kleinstein/immcantation/src/master/training/). [10x Genomics VDJ Sequence Analysis Tutorial](https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html) with the Docker image

## Cell-cell interactions

- [LIANA](https://github.com/saezlab/liana) - LIgand-receptor ANalysis frAmework for cell-cell communication analysis of scRNA-seq clusters. Comparing 16 resources and 7 methods (Table 1: CellChat, CellPhoneDB, Connectome, Crosstalk scores, logFC meam, NATMI, SingleCellSignalR, Consensus). Decoupling methods and resources, allowing for combination of any of them. [Website](https://saezlab.github.io/liana/index.html) with several tutorials. R package. <details>
    <summary>Paper</summary>
    Dimitrov, Daniel, Dénes Türei, Martin Garrido-Rodriguez, Paul L. Burmedi, James S. Nagai, Charlotte Boys, Ricardo O. Ramirez Flores, et al. “Comparison of Methods and Resources for Cell-Cell Communication Inference from Single-Cell RNA-Seq Data.” Nature Communications 13, no. 1 (June 9, 2022): 3224. https://doi.org/10.1038/s41467-022-30755-0.
</details>

- [Connectome](https://github.com/msraredon/Connectome) - calculation and exploration of cell-cell signaling network topologies (centrality, hub scores, differential connectome) in scRNA-seq data (Seurat objects). Can be used with any ligand-receptor data. Visualization capabilities. [Website](https://msraredon.github.io/Connectome/) with several tutorials. R package. <details>
    <summary>Paper</summary>
    Raredon, Micha Sam Brickman, Junchen Yang, James Garritano, Meng Wang, Dan Kushnir, Jonas Christian Schupp, Taylor S. Adams, et al. “Computation and Visualization of Cell–Cell Signaling Topologies in Single-Cell Systems Data Using Connectome.” Scientific Reports 12, no. 1 (March 9, 2022): 4187. https://doi.org/10.1038/s41598-022-07959-x.
</details>

- [OmniPath](https://omnipathdb.org/) - database of intra- and extracellular signaling, as well as transcriptional and post-transcriptional regulation. Combining [over 100 resources](https://omnipathdb.org/info). Human data, converted to mouse and rat via homology. [OmnipathR](https://github.com/saezlab/OmnipathR) [Bioconductor](https://bioconductor.org/packages/OmnipathR/) R package, [omnipath](https://github.com/saezlab/omnipath) Python client, [Cytoscape](https://apps.cytoscape.org/apps/omnipath) plugin. [Code used in paper](https://github.com/saezlab/omnipath_analysis). <details>
    <summary>Paper</summary>
    Türei, Dénes, Alberto Valdeolivas, Lejla Gul, Nicolàs Palacio‐Escat, Michal Klein, Olga Ivanova, Márton Ölbei, et al. “Integrated Intra‐ and Intercellular Signaling Knowledge for Multicellular Omics Analysis.” Molecular Systems Biology 17, no. 3 (March 2021). https://doi.org/10.15252/msb.20209923.
</details>

- [CellChat](http://www.cellchat.org/) - web tool and R package for cell-cell communication network inference from single cell DNA-seq. Predicts signaling inputs and outputs using network analysis and pattern recognition (graph theory). Comparative analysis of signaling pathways between datasets. Visualization (hierarchical plot, circle and bubble plots). At the heart - manually curated CellChatDB. Other methods: SingleCellSignalR, iTALK, NicheNet, CellPhoneDB. Demonstrated on mouse and human skin datasets. [GitHub](https://github.com/sqjin/CellChat). <details>
    <summary>Paper</summary>
    Jin, S., C. F. Guerrero-Juarez, L. Zhang, I. Chang, R. Ramos, C. H. Kuan, P. Myung, M. V. Plikus, and Q. Nie. "Inference and analysis of cell-cell communication using CellChat."" Nat. Commun. 12, https://doi.org/10.1038/s41467-021-21246-9 (17 February 2021).
</details>

- [CellChat](https://github.com/jinworks/CellChat) v2 R protocol. CellChat quantifies the signaling communication probability between two cell groups using a simplified mass action-based model. Three unique features: (i) incorporation of soluble and membrane-bound stimulatory and inhibitory cofactors, (ii) classification of ligand-receptor pairs into functionally related signaling pathways, and (iii) rich annotations of each ligand-receptor pair. Incorporates spatial location of cells, if available. Visualization using CellChat Explorer. Input: gene expression matrix and cell type annotations, from text files, AnnData (scanpy) or Seurat objects. [Data for the protocol](https://figshare.com/projects/Example_data_for_cell-cell_communication_analysis_using_CellChat/157272). <details>
    <summary>Paper</summary>
    Jin, Suoqin, Maksim V. Plikus, and Qing Nie. “CellChat for Systematic Analysis of Cell-Cell Communication from Single-Cell and Spatially Resolved Transcriptomics.” Preprint. Bioinformatics, November 5, 2023. https://doi.org/10.1101/2023.11.05.565674.
</details>

- [CellPhoneDB](https://www.cellphonedb.org/) - database and a tool for cell-cell communication analysis. Contains ligands, receptors, their interactions (978 proteins, 1396 interactions). Input - annotated scRNA-seq data, also protein expression, . Permutation-based comparison of mean expression of receptor-ligand coding genes. Examples how to reformat a Seurat object, scanpy adata. The paper contains full tutorial on using the [cellphonedb Python package](https://github.com/ventolab/CellphoneDB. v5 update includes non-protein ligands (endocrine hormones and GPRCs), a differential expression-based methodology, a scoring methodology to prioritize specific cell-cell interactions utilizing other single-cell modalities such as spatial information of TF activities (CellSign module), visualization using [CellPhoneDBViz](https://github.com/datasome/cellphonedbviz/) (D3). [Raw data](https://github.com/ventolab/CellphoneDB-data), [Documentation](https://cellphonedb.readthedocs.io/en/latest/). <details>
    <summary>Paper</summary>
    Efremova, Mirjana. "CellPhoneDB: Inferring Cell–Cell Communication from Combined Expression of Multi-Subunit Ligand–Receptor Complexes" https://doi.org/10.1038/s41596-020-0292-x  NATURE PROTOCOLS, (2020)
    Troulé, Kevin, Robert Petryszak, Martin Prete, James Cranley, Alicia Harasty, Zewen Kelvin Tuong, Sarah A. Teichmann, Luz Garcia-Alonso, and Roser Vento-Tormo. "CellPhoneDB v5: inferring cell-cell communication from single-cell multiomics data." arXiv preprint arXiv:2311.04567 (2023). https://arxiv.org/pdf/2311.04567
</details>




## Simulation

- [scDesign](https://github.com/Vivianstats/scDesign) - scRNA-seq data simulator and statistical framework to access experimental design for differential gene expression analysis. Gamma-Normal mixture model better fits scRNA-seq data, accounts for dropout events (Methods describe step-wise statistical derivations). Single- or double-batch sequencing scenarios. Comparable or superior performance to simulation methods `splat`, `powsimR`, `scDD`, Lun et al. method. DE tested using t-test. Applications include DE methods evaluation, dimensionality reduction testing. <details>
    <summary>Paper</summary>
    Li, Wei Vivian, and Jingyi Jessica Li. "A Statistical Simulator ScDesign for Rational ScRNA-Seq Experimental Design" https://doi.org/10.1093/bioinformatics/btz321  Bioinformatics 35, no. 14 (July 15, 2019)
</details>

- [Splatter](https://github.com/Oshlack/splatter) - scRNA-seq simulator and pre-defined differential expression. 6 methods, description of each. Issues with scRNA-seq data - dropouts, zero inflation, proportion of zeros, batch effect. Negative binomial for simulation. No simulation is perfect. <details>
    <summary>Paper</summary>
    Zappia, Luke, Belinda Phipson, and Alicia Oshlack. "Splatter: Simulation Of Single-Cell RNA Sequencing Data" https://doi.org/10.1186/s13059-017-1305-0  July 24, 2017. 
</details>

### Power

- [SCOPIT](https://alexdavisscs.shinyapps.io/scs_power_multinomial/) - Shiny app for estimating the number of cells that must be sequenced to observe cell types in a single-cell sequencing experiment. By Alexander Davis

- [How many cells](https://satijalab.org/howmanycells) do we need to sample so that we see at least n cells of each type. By Satija's lab.

- [scPower](https://github.com/heiniglab/scPower) - an R package for power calculation for single-cell RNA-seq studies. Estimates power of differential expression and eQTLs using zero-inflated negative binomial distribution. Also, power to detect rare cell types. Figure 1 shows the dependence among experimental design parameters. Tested on several datasets, generalizes well across technologies. [GitHub](https://github.com/heiniglab/scPower) and [Shiny app](http://scpower.helmholtz-muenchen.de/). <details>
    <summary>Paper</summary>
    Schmid, Katharina T., Barbara Höllbacher, Cristiana Cruceanu, Anika Böttcher, Heiko Lickert, Elisabeth B. Binder, Fabian J. Theis, and Matthias Heinig. "ScPower Accelerates and Optimizes the Design of Multi-Sample Single Cell Transcriptomic Studies" https://doi.org/10.1038/s41467-021-26779-7  Nature Communications, (December 2021)
</details>

- [powsimR](https://github.com/bvieth/powsimR) - an R package for simulating scRNA-seq datasets and assess performance of differential analysis methods. Supports Poisson, Negative Binomial, and zero inflated NB, or estimates parameters from user-provided data. Simulates differential expression with pre-defined fold changes, estimates power, TPR, FDR, sample size, and for the user-provided dataset. <details>
    <summary>Paper</summary>
    Vieth, Beate, Christoph Ziegenhain, Swati Parekh, Wolfgang Enard, and Ines Hellmann. "PowsimR: Power Analysis for Bulk and Single Cell RNA-Seq Experiments" https://doi.org/10.1093/bioinformatics/btx435  Edited by Ivo Hofacker. Bioinformatics 33, no. 21 (November 1, 2017)
</details>

### Benchmarking

- [CellBench](https://bioconductor.org/packages/CellBench/) - an R package for benchmarking of scRNA-seq analysis pipelines. Simulated datasets using mixtures of either cells of RNA from five cancer cell lines, dilution series, ERCC spike-in controls. Four technologies. Methods: normalization, imputation, clustering, trajectory analysis, data integration. Evaluation metrics: silhouette width, correlations, others. Best performers: Normalization - Linnorm, scran, scone; Imputation - kNN, DrImpute; Clustering - all methods are OK, Seurat performs well; Trajectory - Slingshot and Monocle2. [Processed datasets used for the analysis](https://github.com/LuyiTian/sc_mixology), [GitHub](https://github.com/Shians/CellBench). <details>
    <summary>Paper</summary>
    Tian, Luyi, Xueyi Dong, Saskia Freytag, Kim-Anh Lê Cao, Shian Su, Abolfazl JalalAbadi, Daniela Amann-Zalcenstein, et al. "Benchmarking Single Cell RNA-Sequencing Analysis Pipelines Using Mixture Control Experiments" https://doi.org/10.1038/s41592-019-0425-8  Nature Methods, May 27, 2019. 
</details>

## Deep learning

- [scGPT](https://github.com/bowang-lab/scGPT/tree/main) - a generative pre-trained transformer for scRNA-seq gene and cell embedding analyses (cell type annotation, multi-batch and multi-omics integration, genetic perturbation prediction, gene network inference). Gene tokens are analogous to word tokens, same for cell tokens, condition tokens. Flash-Attention accelerated self-attention implementation. Pre-trained on over 33 million cells (whole human), as well as on brain, blood, pancreas, lung, heart, kidney, intestine datasets. The pre-trained model parameters can be fine-tuned on new data, improving performance. Benchmarked against transformer-based methods TOSICA and scBert (for cell annotations), scVI, Seurat, Harmony (for batch integration), scGLUE and Seurat (for multi-omics integration). Many datasets used for benchmarking. Pytorch, GPU. [Documentation](https://scgpt.readthedocs.io/en/latest/index.html). <details>
    <summary>Paper</summary>
    Cui, Haotian, Chloe Wang, Hassaan Maan, Kuan Pang, Fengning Luo, and Bo Wang. “ScGPT: Towards Building a Foundation Model for Single-Cell Multi-Omics Using Generative AI.” Preprint. Bioinformatics, May 1, 2023. https://doi.org/10.1101/2023.04.30.538439.
</details>

- [CeLEry](https://github.com/QihuangZhang/CeLEry) (Cell Location recovEry) - spatial location prediction from scRNA-seq data. A feed-forward deep neural network trained on scRNA-seq and spatial transcriptomics. Data augmentation (variational autoencoder) improves performance. Logistic loss or rank-consistent logistic loss function (Methods). Tested on 10x Visium data ([LIBD](http://research.libd.org/spatialLIBD/) data), other human and mouse brain data obtained with different technologies, human liver cancer (MERSCOPE), also a HER2-positive, ESR1-positive, and PGR-negative breast tumor dataset obtained from the 10x Xenium platform. Benchmarked against Tangram, spaOTsc, novoSpaRc. <details>
    <summary>Paper</summary>
    Zhang, Qihuang, Shunzhou Jiang, Amelia Schroeder, Jian Hu, Kejie Li, Baohong Zhang, David Dai, Edward B. Lee, Rui Xiao, and Mingyao Li. “Leveraging Spatial Transcriptomics Data to Recover Cell Locations in Single-Cell RNA-Seq with CeLEry.” Nature Communications 14, no. 1 (July 8, 2023): 4050. https://doi.org/10.1038/s41467-023-39895-3.
</details>

- [scTranslator](https://github.com/TencentAILabHealthcare/scTranslator) - a pre-trained large language model for translating single-cell transcriptome to proteome. Transformer-based, Fast Attention Via positive Orthogonal Random features approach (FAVOR+), details on implementation. Trained on bulk and scRNA-seq data paired with proteomics (TCGA, CPTAC, others). Benchmarked against [cTP-Net](https://github.com/zhouzilu/ctpnetpy/tree/master/ctpnet), [sciPENN](https://github.com/jlakkis/sciPENN). Python (Pytorch), Jupyter, in Docker, GPU required. <details>
    <summary>Paper</summary>
    Liu, Linjing, Wei Li, Ka-chun Wong, Fan Yang, and Jianhua Yao. "A pre-trained large language model for translating single-cell transcriptome to proteome." bioRxiv (2023): 2023-07. https://doi.org/10.1101/2023.07.04.547619
</details>

- [scBERT](https://github.com/TencentAILabHealthcare/scBERT) - a pretrained deep neural network, single-cell bidirectional encoder representations from transformers model, for cell type annotation of scRNA-seq data. Gene embeddings are obtained from gene2vec (binned gene expression); such embeddings capture the co-expression similarity (analogous to semantic similarity). Performer to allow for over 16,000 gene inputs. Benchmarked against marker-based (SCINA, Garnett, scSorter), correlation-based (Seurat, SingleR, scmap_cell, scmap_cluster, Cell_ID), and machine learning methods (SciBet, scNym). Robust to class imbalance, performs well across cohorts and organs, discovers novel cell types. <details>
    <summary>Paper</summary>
    Yang, Fan, Wenchuan Wang, Fang Wang, Yuan Fang, Duyu Tang, Junzhou Huang, Hui Lu, and Jianhua Yao. “ScBERT as a Large-Scale Pretrained Deep Language Model for Cell Type Annotation of Single-Cell RNA-Seq Data.” Nature Machine Intelligence 4, no. 10 (September 26, 2022): 852–66. https://doi.org/10.1038/s42256-022-00534-z.
</details>

- <a name="saucie">[SAUCIE](https://github.com/KrishnaswamyLab/SAUCIE/)</a> - a regularized autoencoder for scRNA-seq data denoising, batch correction, low-dimensional representation and clustering. Three encoding layers (512, 256, 128 neurons), an embedding layer (two dimensional), three decoding layers. Information dimension (ID) regularization encourages binarizable activation of the neurons (helps clustering). A maximal mean discrepancy (MMD) to penalize differences between probability distributions of internal activations of samples (for batch correction across samples). Outperforms many clustering methods, scalable. <details>
    <summary>Paper</summary>
    Amodio, Matthew, David van Dijk, Krishnan Srinivasan, William S. Chen, Hussein Mohsen, Kevin R. Moon, Allison Campbell, et al. “Exploring Single-Cell Data with Deep Multitasking Neural Networks.” Nature Methods, October 7, 2019. https://doi.org/10.1038/s41592-019-0576-7.
</details>

- [scScope](https://github.com/AltschulerWu-Lab/scScope) - for extracting informative representations, clustering and analysing cell type composition. Recurrent neural network architecture, variable number of recurrent steps (at step = 1, the architecture is standard autoencoder). Tested on simulated data (Splatter, SIMLR) and four experimental scRNA-seq dtasets, at different sparsity levels, rare subpopulation fractions. Compared with PCA, ZINB-WaVE, scVI, others. Multi-GPU Python, TensorFlow, numpy, scikit-learn implementation. <details>
    <summary>Paper</summary>
    Deng, Yue, Feng Bao, Qionghai Dai, Lani F. Wu, and Steven J. Altschuler. “Scalable Analysis of Cell-Type Composition from Single-Cell Transcriptomics Using Deep Recurrent Learning.” Nature Methods, March 18, 2019. https://doi.org/10.1038/s41592-019-0353-7
</details>

- [Solo](https://github.com/calico/Solo) - semi-supervised deep learning for doublet identification. Variational autoencoder (scVI) followed by a classifier to detect doublets. Compared with Scrubled and DoubletFinder, improves area under the precision-recall curve. <details>
    <summary>Paper</summary>
    Bernstein, Nicholas J., Nicole L. Fong, Irene Lam, Margaret A. Roy, David G. Hendrickson, and David R. Kelley. "Solo: Doublet Identification in Single-Cell RNA-Seq via Semi-Supervised Deep Learning" https://doi.org/10.1016/j.cels.2020.05.010  Cell Systems 11, no. 1 (July 2020)
</details>

- [scover](https://github.com/jacobhepkema/scover) - de novo identification of regulatory motifs and their cell type-specific importance from scRNA-seq or scATAC-seq data. Shallow convolutional neural network on one-hot encoded sequence data, k-fold training and selecting most optimal network, extracting motifs from convolutional filters, cluster them, matching with motifs, associating with peak strength/gene expression. application for human kidney scRNA-seq data, Tabula Muris, mouse cerebral cortex SNARE-seq data. [Docs](https://scover.readthedocs.io/en/latest/), [Tweet](https://twitter.com/jhepkema/status/1332444937493897216?s=20). <details>
    <summary>Paper</summary>
    Hepkema, Jacob, Nicholas Keone Lee, Benjamin J Stewart, Siwat Ruangroengkulrith, Varodom Charoensawan, Menna R Clatworthy, and Martin Hemberg. "Predicting the Impact of Sequence Motifs on Gene Regulation Using Single-Cell Data" https://doi.org/10.1101/2020.11.26.400218 
</details>

- [SAVER-X](https://github.com/jingshuw/SAVERX) - denoising scRNA-seq data using deep autoencoder with a Bayesian model. Decomposes the variation into three components: 1) predictable, 2) unpredictable, 3) technical noise. Pretrained on the Human Cell Atlas project, 10X Genomics immune cells, allows for human-mouse cross-species learning. Improves clustering and the detection of differential genes. Outperforms downsampling, MAGIC, DCA, scImpute. <details>
    <summary>Paper</summary>
    Littmann, Maria, Katharina Selig, Liel Cohen-Lavi, Yotam Frank, Peter Hönigschmid, Evans Kataka, Anja Mösch, et al. "Validity of Machine Learning in Biology and Medicine Increased through Collaborations across Fields of Expertise" https://doi.org/10.1038/s42256-019-0139-8  Nature Machine Intelligence, January 13, 2020. 
</details>

- [scVI](https://github.com/YosefLab/scVI) - low-dimensional representation of scRNA-seq data used for batch correction, imputation, clustering, differential expression. Deep neural networks to approximate the distribution that underlie observed expression values. Zero-inflated negative binomial distribution conditioned on the batch annotation and unobserved random variables. Compared with DCA, ZINB-WAVE on simulated and real large and small datasets. [Perspective by Way & Greene](https://www.nature.com/articles/s41592-018-0230-9). <details>
    <summary>Paper</summary>
    Lopez, Romain, Jeffrey Regier, Michael B Cole, Michael Jordan, and Nir Yosef. "Bayesian Inference for a Generative Model of Transcriptome Profiles from Single-Cell RNA Sequencing" https://doi.org/10.1101/292037  September 23, 2018. 
</details>

## Spatial transcriptomics

- [awesome_spatial_omics](https://github.com/crazyhottommy/awesome_spatial_omics) - tools and notes for spatial omics, by Ming Tang

- [PrinciplesSTA](https://github.com/lmweber/PrinciplesSTA) - Principles of Spatial Transcriptomics Analysis with Bioconductor" book, [web version](https://lmweber.org/PrinciplesSTA/devel/)

- Overview of spatial technologies. Table 1 - comparison of in situ capture-based technologies, imaging-based, and region of interest-based. Table 2 - tools for all inclusive analysis, deconvolution, CNV inference, cellular interaction, deep learning. Cellular architecture of tumor microenvironment, clinical use. <details>
    <summary>Paper</summary>
    Chen, Julia, Ludvig Larsson, Alexander Swarbrick, and Joakim Lundeberg. “Spatial Landscapes of Cancers: Insights and Opportunities.” Nature Reviews Clinical Oncology 21, no. 9 (September 2024): 660–74. https://doi.org/10.1038/s41571-024-00926-7.
</details>

- [Sopa](https://github.com/gustaveroussy/sopa) - a technology-invariant, memory-efficient pipeline with a unified visualizer for all image-based spatial omics. Includes segmentation, multi-level annotation, spatial statistics, niche geometry analysis. SpatialData object structure. Output compatible with the Xenium Explorer. Benchmarked on two spatial transcriptomics (MERISCOPE, Xenium) and two multiplex imaging technologies (PhenoCycler, MACSima). Implemented as a Snakemake pipeline, CLI, Python API. <details>
    <summary>Paper</summary>
    Blampey, Quentin, Kevin Mulder, Charles-Antoine Dutertre, Margaux Gardet, Fabrice Andre, Florent Ginhoux, and Paul-Henry Cournede. “Sopa: A Technology-Invariant Pipeline for Analyses of Image-Based Spatial-Omics,” 2024, https://doi.org/10.1038/s41467-024-48981-z
</details>

- [10x Xenium Human Breast Tumors Explorer](https://www.10xgenomics.com/products/xenium-in-situ/human-breast-dataset-explorer) - interactive zoomable explorer of Xenium In Situ spatial plots. Single-cell spatial transcriptomics of the breast cancer tumor microenvironment, FFPE tissue. FFPE-compatible single cell gene expression workflow (Chromium Fixed RNA Profiling; scFFPE-seq), spatial transcriptomics (Visium CytAssist), and automated microscopy-based in situ technology using a 313-plex gene panel (Xenium In Situ). Xenium description, can be scaled to 1000 gene panel, area up to 2.8cm^2, non-destructive, allows for post-Xenium immunofluorescence and H&E analysis, more sensitive than scFFPE. 17 different cell types, Xenium resolved them with single cell resolution. Integrative analysis across technologies identified three tumor subtypes (low- and high-grade ductal carcinoma in situ, invasive carcinoma). Standalone [Xenium explorer](https://www.10xgenomics.com/support/software/xenium-explorer/downloads), [tutorials](https://www.10xgenomics.com/support/software/xenium-explorer/tutorials). [Downloadable datasets](https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast). <details>
    <summary>Paper</summary>
    Janesick, Amanda, Robert Shelansky, Andrew D. Gottscho, Florian Wagner, Morgane Rouault, Ghezal Beliakoff, Michelli Faria De Oliveira, et al. “High Resolution Mapping of the Breast Cancer Tumor Microenvironment Using Integrated Single Cell, Spatial and in Situ Analysis of FFPE Tissue.” Preprint. Cancer Biology, October 7, 2022. https://doi.org/10.1101/2022.10.06.510405.
</details>

- [SpatialDecon](https://bioconductor.org/packages/SpatialDecon/) - cell type deconvolution of spatial gene expression (Nanostring GeoMx platform). Log-normal regression instead of the least squares regression, performs better that four other methods. [Supplementary Data S1](https://www.nature.com/articles/s41467-022-28020-5#Sec39) - gene markers for immune and stromal cell types. [GitHub](https://github.com/Nanostring-Biostats/SpatialDecon). <details>
    <summary>Paper</summary>
    Danaher, Patrick, Youngmi Kim, Brenn Nelson, Maddy Griswold, Zhi Yang, Erin Piazza, and Joseph M. Beechem. “Advances in Mixed Cell Deconvolution Enable Quantification of Cell Types in Spatial Transcriptomic Data.” Nature Communications 13, no. 1 (January 19, 2022): 385. https://doi.org/10.1038/s41467-022-28020-5.
</details>

- [CosMx spatial molecular imaging technology](https://nanostring.com/products/cosmx-spatial-molecular-imager/single-cell-imaging-overview/) for spatial profiling of more than 1,000 RNA and 64 proteins at single-cell level. XYZ coordinates. Highly correlated with scRNA-seq. Fresh-frozen and paraffin-embedded. [Human liver and non-small lung cancer data](https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/). <details>
    <summary>Paper</summary>
    He, Shanshan, Ruchir Bhatt, Carl Brown, Emily A. Brown, Derek L. Buhr, Kan Chantranuvatana, Patrick Danaher, et al. “High-Plex Imaging of RNA and Proteins at Subcellular Resolution in Fixed Tissue by Spatial Molecular Imaging.” Nature Biotechnology 40, no. 12 (December 2022): 1794–1806. https://doi.org/10.1038/s41587-022-01483-z.
</details>

- [SpatialData](https://github.com/scverse/spatialdata-notebooks) - a community standards-based framework to store, process and annotate multi-platform spatial data. Handles larger-than-memory data and transformations, Zarr file format. Demo on multi-modal Xenium and Visium breast cancer study. Python implementation, integrates with other multi-omics ecosystems (Squidpy, Scanpy, MONAI, scvi-tools). [Documentation](https://spatialdata.scverse.org/en/latest/index.html). <details>
    <summary>Paper</summary>
    Marconato, Luca, Giovanni Palla, Kevin A. Yamauchi, Isaac Virshup, Elyas Heidari, Tim Treis, Marcella Toth, et al. “SpatialData: An Open and Universal Data Framework for Spatial Omics.” Preprint. Bioinformatics, May 8, 2023. https://doi.org/10.1101/2023.05.05.539647.
</details>

- [STELLAR](http://snap.stanford.edu/stellar/) (SpaTial cElL LeARning) - a geometric deep learning method for cell-type discovery and identification in spatially resolved single-cell datasets. Input - annotated spatially-resolved single-cell dataset. Graph convolutional network learns latent low-dimensional cell representations spatially and molecularly close. Transfers learned representations on an unannotated dataset. Compared with baseline methods XGBoost, SVM, RF, AdaBoost, as well as single-cell annotation baselines Seurat, scNym, scANVI, outperforms all (Seurat and scNym perform poorly). Robust to normalization strategies, noise, missing values. PyTorch implementation, [GitHub](https://github.com/snap-stanford/stellar). <details>
    <summary>Paper</summary>
    Brbić, Maria, Kaidi Cao, John W. Hickey, Yuqi Tan, Michael P. Snyder, Garry P. Nolan, and Jure Leskovec. “Annotation of Spatially Resolved Single-Cell Data with STELLAR.” Nature Methods 19, no. 11 (November 2022): 1411–18. https://doi.org/10.1038/s41592-022-01651-8.
</details>

- [SpatialPCA](https://github.com/shangll123/SpatialPCA) - spatially aware dimention reduction method, a low dimensional representation of the spatial transcriptomics data with biological signal and preserved spatial correlation structure. Used for spatial domain detection. Outperforms BayesSpace, SpaGCN, HMRF, PCA, NMF in simulation study, better spatial continuity and smoothness, ARI, CHAOS, LISI coefficients. [Code to reproduce the paper](https://github.com/shangll123/SpatialPCA_analysis_codes). [Tutorials](https://lulushang.org/SpatialPCA_Tutorial/index.html). <details>
    <summary>Paper</summary>
    Shang, Lulu, and Xiang Zhou. “Spatially Aware Dimension Reduction for Spatial Transcriptomics.” Nature Communications 13, no. 1 (November 23, 2022): 7203. https://doi.org/10.1038/s41467-022-34879-1.
</details>

- [Stereo-seq](https://db.cngb.org/stomics/mosta/stereo.seq.html) spatial transcriptomics technology. Combines DNA nanoball-patterned (DNB) arrays (random barcodes on a patterned attay so each DNB contains the coordinate identity) and tissue RNA capture, achieves cellular resolution. Generates the [mouse organogenesis spatiotemporal transcriptomic atlas (MOSTA)](https://db.cngb.org/stomics/mosta/) to map the spatiotemporal transcriptomic dynamics of the developing mouse embryo. [Raw data](https://db.cngb.org/search/project/CNP0001543/). [SAW](https://github.com/BGIResearch/SAW) - Workflow for analyzing Stereo-seq transcriptomic data. <details>
    <summary>Paper</summary>
    Chen, Ao, Sha Liao, Mengnan Cheng, Kailong Ma, Liang Wu, Yiwei Lai, Xiaojie Qiu, et al. “Spatiotemporal Transcriptomic Atlas of Mouse Organogenesis Using DNA Nanoball-Patterned Arrays.” Cell 185, no. 10 (May 2022): 1777-1792.e21. https://doi.org/10.1016/j.cell.2022.04.003.
</details>

- [NovoSpaRc](https://pypi.org/project/novosparc/) Python package and protocol for reconstruction of spatial positioning of scRNA-seq data. Structural correspondence hypothesis (optimal transport), cells in physical proximity share similar expression profiles. Reference atlas is optional and improves the spatial reconstruction. Input: gene expression matrix and a target space (1D, 2D, or 3D coordinates of the physical space, defaults available). Reference atlas expression is optional. Calculates three cost matrices (cell-cell, location-location, reference atlas), outputs a transport matrix (probabilistic mapping of cells onto the target locations) and the inferred gene expression over the target space. Table 1 - comparison with Seurat, DistMap, Perler, Tangram, CSOmap. Application demo, Python code. [GitHub](https://github.com/rajewsky-lab/novosparc). <details>
    <summary>Paper</summary>
    Moriel, Noa, Enes Senel, Nir Friedman, Nikolaus Rajewsky, Nikos Karaiskos, and Mor Nitzan. “NovoSpaRc: Flexible Spatial Reconstruction of Single-Cell Gene Expression with Optimal Transport.” Nature Protocols 16, no. 9 (September 2021): 4177–4200. https://doi.org/10.1038/s41596-021-00573-7.

    Nitzan, Mor, Nikos Karaiskos, Nir Friedman, and Nikolaus Rajewsky. “Gene Expression Cartography.” Nature 576, no. 7785 (December 5, 2019): 132–37. https://doi.org/10.1038/s41586-019-1773-3.
</details>

- [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) - R/Bioconductor package providing data infrastructure (S4 class) for spatial transcriptomics data. Intro into spot-based and molecule-based spatial transcriptomics technologies. `read10xVisium()` reads in SpaceRanger-processed data. Example datasets and visualization tools in the [STexampleData](https://github.com/lmweber/STexampleData), [TENxVisiumData](https://github.com/HelenaLC/TENxVisiumData), and [ggspavis](https://github.com/lmweber/ggspavis) packages. <details>
    <summary>Paper</summary>
    Righelli, Dario, Lukas M Weber, Helena L Crowell, Brenda Pardo, Leonardo Collado-Torres, Shila Ghazanfar, Aaron TL Lun, Stephanie C Hicks, and Davide Risso. “SpatialExperiment: Infrastructure for Spatially Resolved Transcriptomics Data in R Using Bioconductor.” Preprint. Bioinformatics, January 27, 2021. https://doi.org/10.1101/2021.01.27.428431.
</details>

- [Spatial transcriptomics analysis and tools](https://htmlpreview.github.io/?https://github.com/drieslab/awesome-spatial-data-analysis/blob/main/review_spat_trns_methods.html) review. Commercial technologies (10X Visium, NanoString GeoMx). Cell type annotation, deconvolution, mapping annotations to spatial coordinates, spatial segmentation of transcriptomic profiles. cell-cell communication, integrative approaches, spatial visualization. Software tools/packages for each analysis, brief description of each. [GitHub: Awesome Spatial transcriptomics methods](https://github.com/drieslab/awesome-spatial-data-analysis). <details>
    <summary>Paper</summary>
    Dries, Ruben, Jiaji Chen, Natalie del Rossi, Mohammed Muzamil Khan, Adriana Sistig, and Guo-Cheng Yuan. "Advances in Spatial Transcriptomic Data Analysis" https://doi.org/10.1101/gr.275224.121  Genome Research, (October 2021)
</details>

- Review of spatial single-cell transcriptomics technologies. Categorized as 1)microdissection-based, 2) in situ hybridization, 3) in situ sequencing, 4) in situ capturing, 5) in silico reconstruction. Timeline (Figure 1), summary of technologies (Table 1), details of each technology, studies, illustrations. <details>
    <summary>Paper</summary>
    Asp, Michaela, Joseph Bergenstråhle, and Joakim Lundeberg. "Spatially Resolved Transcriptomes—Next Generation Tools for Tissue Exploration" https://doi.org/10.1002/bies.201900221  BioEssays, May 4, 2020
</details>

- [spatialLIBD](https://bioconductor.org/packages/spatialLIBD/) - an R package for 10X Visium spatial transcriptomics data manipulation and visualization. Can handle multiple samples, in contrast to Loupe and Giotto pipelines. Includes demo data. Integrates with Bioconductor via SpatialExperiment class. [Shiny app](http://spatial.libd.org/spatialLIBD/). <details>
    <summary>Paper</summary>
    Pardo, Brenda, Abby Spangler, Lukas M Weber, Stephanie C Hicks, Andrew E Jaffe, Keri Martinowich, and Kristen R Maynard. "SpatialLIBD: An R/Bioconductor Package to Visualize Spatially-Resolved Transcriptomics Data" https://doi.org/10.1101/2021.04.29.440149  biorXiv, April 30, 2021
</details>

- [PASTE](https://github.com/raphael-group/paste) - Probabilistic Alignment of Spatial Transcriptomics Experiments, to align and integrate ST data across adjacent tissue slices. Uses spatial coordinate and transcriptomic similarity. Fused Gromov-Wasserstein Optimal Transport, NMF, methodological. <details>
    <summary>Paper</summary>
    Zeira, Ron, Max Land, and Benjamin J Raphael. "Alignment and Integration of Spatial Transcriptomics Data" https://doi.org/10.1101/2021.03.16.435604  Preprint, March 16, 2021
</details>

- [Giotto](http://spatial.rc.fas.harvard.edu/) - an R framework for spatial transcriptomics analysis pipeline. Two modules, analysis and visualization. Custom S4 data structure. QC, preprocessing, feature selection, dimensionality reduction, clustering, marker gene identification, spatial grid and neighborhood networks HMRF for spatial gene expression patterns detection. Visualization using tSNE, physical, physicalsimple panels (figures in the paper). Applied to the seqFISH+ dataset with 10K genes profiled in 913 cells. <details>
    <summary>Paper</summary>
    Dries R, Zhu Q, Dong R, Eng CH, Li H, Liu K, Fu Y, Zhao T, Sarkar A, Bao F, George RE. "[Giotto: a toolbox for integrative analysis and visualization of spatial expression data. Genome biology" https://doi.org/10.1186/s13059-021-02286-2   08 March 2021
</details>

- [stereoscope](https://github.com/almaan/stereoscope) - Spatial transcriptomics deconvolution using scRNA-seq data. New framework for proportion estimation of each cell type at every capture location. Uses negative binomial distribution for both spatial and scRNA-seq data, parameters estimation from data. Utilizes complete expression profiles instead of marker genes. Input: scRNA-seq matrix, annotations, spatial data matrix (paired data not required). Applied to [mouse brain](http://mousebrain.org/) and human developmental heart spatial data. Tested on synthetic data, modelled after experimental, outperforms [DWLS](https://github.com/dtsoucas/DWLS) and [deconSeq](https://github.com/rosedu1/deconvSeq). [Scripts and data for analyses](https://zenodo.org/record/3951884). Python implementation. <details>
    <summary>Paper</summary>
    Andersson, Alma, Joseph Bergenstråhle, Michaela Asp, Ludvig Bergenstråhle, Aleksandra Jurek, José Fernández Navarro, and Joakim Lundeberg. “Single-Cell and Spatial Transcriptomics Enables Probabilistic Inference of Cell Type Topography.” Communications Biology 3, no. 1 (December 2020): 565. https://doi.org/10.1038/s42003-020-01247-y.
</details>

- [Spatial transcriptomics of the human dorsolateral prefrontal cortex](http://research.libd.org/spatialLIBD/) (12 samples from 3 DLPFC regions). 10X Genomics Visium platform. Layers with distinct expression signatures. Integrated with other scRNA-seq data, neuropsychiatric disorder gene sets. All data are public and accessible via a [Shiny app](http://research.libd.org/spatialLIBD/) and through an [spatialLIBD R package](https://github.com/LieberInstitute/spatialLIBD). [Other code](https://github.com/LieberInstitute/HumanPilot). [Tweet](https://twitter.com/fellgernon/status/1233661576433061888?s=20). <details>
    <summary>Paper</summary>
    Maynard, Kristen E, Leonardo Collado-Torres, Lukas M. Weber, Cedric Uytingco, Brianna K. Barry, Stephen R. Williams, Joseph L. Catallini, et al. "Transcriptome-Scale Spatial Gene Expression in the Human Dorsolateral Prefrontal Cortex" https://doi.org/10.1101/2020.02.28.969931  Preprint. Neuroscience, February 28, 2020. 
</details>

- [Seurat](http://satijalab.org/seurat/) - single-cell RNA-seq for spatial cellular localization, using in situ RNA patterns as a reference. Impute landmark genes, relate them to the reference map. [Tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html), and [Dave Tang notes](https://davetang.org/muse/2017/08/01/getting-started-seurat/). <details>
    <summary>Paper</summary>
    Satija, Rahul, Jeffrey A. Farrell, David Gennert, Alexander F. Schier, and Aviv Regev. "Spatial Reconstruction of Single-Cell Gene Expression Data" https://doi.org/10.1038/nbt.3192  Nature Biotechnology 33, no. 5 (May 2015)
</details>

- [spatialomics.net](http://spatialomics.net/) - Zoom/recorded seminar series on spatial transcriptomics

- Analysis and visualization of spatial transcriptomics data using scanpy, 10X Genomics Visium and MERFISH data, [Jupyter notebook](https://nbviewer.jupyter.org/github/theislab/scanpy-tutorials/blob/master/analysis-visualization-spatial.ipynb), [Tweet by Fabian Theis](https://twitter.com/fabian_theis/status/1224741146242572289?s=20)

- [brainmapr](https://github.com/hms-dbmi/brainmapr) - R package to infer spatial location of neuronal subpopulations within the developing mouse brain by integrating single-cell RNA-seq data with in situ RNA patterns from the Allen Developing Mouse Brain Atlas.

- [Spatial Gene Expression, Space Ranger by 10X Genomics](https://support.10xgenomics.com/spatial-gene-expression/software/overview/welcome)

- [STUtility](https://github.com/jbergenstrahle/STUtility) - an R package for visualization and analysis of spatial transcriptomics data (10X Visium). Seurat-compatible. [Website](https://ludvigla.github.io/STUtility_web_site/index.html) with tutorials and walkthroughs.

## Technology

- [Collections of library structure and sequence of popular single cell genomic methods](https://teichlab.github.io/scg_lib_structs/) from Sarah Teichmann's group. [GitHub](https://github.com/Teichlab/scg_lib_structs)

- Comparison of six scRNA-seq technologies (CEL-seq2, Drop-seq, MARS-seq, SCRB- seq, Smart-seq, and Smart-seq2) using mouse ESCs. Description of each technology, UMI-based and full length. Power, sensitivity analysis. Table 1 compares TPR, FDR, number of cells, library cost factor. Drop-seq is the cheapest at the expense of higher FDR. <details>
    <summary>Paper</summary>
    Ziegenhain, Christoph, Beate Vieth, Swati Parekh, Björn Reinius, Amy Guillaumet-Adkins, Martha Smets, Heinrich Leonhardt, Holger Heyn, Ines Hellmann, and Wolfgang Enard. “Comparative Analysis of Single-Cell RNA Sequencing Methods.” Molecular Cell 65, no. 4 (February 16, 2017): 631-643.e4. https://doi.org/10.1016/j.molcel.2017.01.023.
</details>

- **SPARC** (Single-Cell Protein And RNA Co-profiling) - combines scRNA-seq (improved Smart-seq2 protocol) with proximity extension assays (PEA) to simultaneously measure global mRNA and 89 intracellular proteins in individual cells. Description of PEA - a proximity-based assay that require two binding events to generate a DNA reporter molecule (Figure 1). Applied to human embryonic stem cells following directed neural induction (from 0h to 48h). mRNA data highly correlate with regular Smart-seq2 data, not so with protein expression. Temporal reconstruction using mRNA and protein-based profiles produces similar results, suggesting that protein changes follow mRNA changes. [Raw](https://doi.org/10.17044/scilifelab.14207462) and [processed](https://doi.org/10.17044/scilifelab.14207909) data at SciLifeLab Data Repository. [GitHub](https://github.com/b97jre/SPARC). <details>
    <summary>Paper</summary>
    Reimegård, Johan, Marcel Tarbier, Marcus Danielsson, Jens Schuster, Sathishkumar Baskaran, Styliani Panagiotou, Niklas Dahl, Marc R. Friedländer, and Caroline J. Gallant. “A Combined Approach for Single-Cell MRNA and Intracellular Protein Expression Analysis.” Communications Biology 4, no. 1 (December 2021): 624. https://doi.org/10.1038/s42003-021-02142-w
</details>

- **SHARE-seq** - simultaneous profiling of scRNA-seq and sc-ATAC-seq from the same cells. Built upon SPLiT-seq, a combinatorial indexing method. Confirmed by separate scRNA-seq and scATAC-seq datasets. Lineage-determining genes are marked by domains of regulatory chromatin (DORCs). Chromatin accessibility at DORCs precedes gene expression. [Visualization of SHARE-seq skin data](https://buenrostrolab.shinyapps.io/skinnetwork/) <details>
    <summary>Paper</summary>
    Ma, Sai, Bing Zhang, Lindsay M. LaFave, Andrew S. Earl, Zachary Chiang, Yan Hu, Jiarui Ding, et al. “Chromatin Potential Identified by Shared Single-Cell Profiling of RNA and Chromatin.” Cell 183, no. 4 (November 2020): 1103-1116.e20. https://doi.org/10.1016/j.cell.2020.09.056.
</details>

- **snmC-seq2** - single-cell DNA methylation, improved read mapping, reduced artifacts, enhanced throughput, library complexity, coverage over previously developed snmC-seq. Uses random primers with A/T/C nucleotides. Protocol in the [supplementary methods](https://www.nature.com/articles/s41467-018-06355-2#Sec15). Raw data on [GSE112471](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112471), [Data visualization](http://neomorph.salk.edu/snmC-eq2.php). Code on [GitHub](https://github.com/zhoujt1994/Luo2018_ncomm). <details>
    <summary>Paper</summary>
    Luo, Chongyuan, Angeline Rivkin, Jingtian Zhou, Justin P. Sandoval, Laurie Kurihara, Jacinta Lucero, Rosa Castanon, et al. “Robust Single-Cell DNA Methylome Profiling with SnmC-Seq2.” Nature Communications 9, no. 1 (December 2018): 3824. https://doi.org/10.1038/s41467-018-06355-2.
</details>

- **Drop-seq technology** - single cells encapsulated in lipid droplets with nanoparticles with cell- and UMI barcodes. Barcoding strategy, "split-and-pool" synthesis cycles to synthesize 12bp cell barcodes, then 8bp UMI synthesis (Figure 1). Majority of droplets are empty, doublets depend on initial cell concentration. Example on a mixture of 589 human HEK and 412 mouse 3T3 cells. Expression profiles from 49,300 retinal cells profiled using Drop-seq. 13,155 largest libraries, reduce dimensionality by PCA to 32 components (decided by permutation), tSNE for visualization. 39 clusters matched to known cell types. [GEO GSE63473](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63473). <details>
    <summary>Paper</summary>
    Macosko, Evan Z., Anindita Basu, Rahul Satija, James Nemesh, Karthik Shekhar, Melissa Goldman, Itay Tirosh, et al. "Highly Parallel Genome-Wide Expression Profiling of Individual Cells Using Nanoliter Droplets" https://doi.org/10.1016/j.cell.2015.05.002  Cell, (May 21, 2015)
</details>

### 10X Genomics

- [10X Genomics list of publications](https://www.10xgenomics.com/resources/publications/), categorized by species, tissue types, journals, technology

- [List of tools and resources related to the 10x Genomics GEMCode/Chromium system](https://github.com/johandahlberg/awesome-10x-genomics)

- Zheng, Grace X. Y., Jessica M. Terry, Phillip Belgrader, Paul Ryvkin, Zachary W. Bent, Ryan Wilson, Solongo B. Ziraldo, et al. "[Massively Parallel Digital Transcriptional Profiling of Single Cells](https://doi.org/10.1038/ncomms14049)"  Nature Communications 8 (January 16, 2017) - 10X technology. Details of each wet-lab step, sequencing, and basic computational analysis. Calling SNPs from scRNA-seq data. Reduce dimensionality with PCA (50 PCs), K-means to assign cluster labels, visualizing with tSNE. [Code for the paper](https://github.com/10XGenomics/single-cell-3prime-paper). [scRNA-Seq dataset of 3000 peripheral blood mononuclear cells (PBMCs)](https://support.10xgenomics.com/single-cell-gene-expression/datasets)

- [Cell Ranger, Loupe Cell Browser software download](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest), 10X Genomics

- [Cell Ranger R kit](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/rkit) by 10X Genomics

- [Getting started with Cell Ranger](https://davetang.org/muse/2018/08/09/getting-started-with-cell-ranger/) by Dave Tang

- [An R script for importing the fragments file from the CellRanger output and making a SummarizedExperiment](https://github.com/caleblareau/scATAC_10X_raw_to_kmers/blob/master/example_kmers.R), by Caleb Lareau. 

#### 10X QC

- [10xqc.com](https://10xqc.com/) - Submit your 10X Cell Ranger® report and compare to data from 155 other reports, contributed by other users from across the globe. This tool was developed to allow users of 10X Genomics single-cell 3'mRNA-seq technology to share their experiences. Requires `web_summary.html` upload. [GitHub](https://github.com/s-andrews/10XQC).

- [bxcheck](https://github.com/pd3/bxcheck) - Toolset for QC and processing 10x genomics data. https://github.com/pd3/bxcheck

- [What is sequencing saturation?](https://kb.10xgenomics.com/hc/en-us/articles/115005062366-What-is-sequencing-saturation-) by 10X Genomics

## Data

- [SingleCellMultiModal](https://bioconductor.org/packages/SingleCellMultiModal/) - an R data package with landmark datasets from important single-cell multimodal protocols, including CITE-Seq, ECCITE-Seq, SCoPE2, scNMT, 10X Multiome, seqFISH, and G&T. RNA, protein, open chromatin, methylation, spatial data modalities. Data in MultiAssayExperiment format. Mouse, human. Blood, brain, fibroblast cells. <details>
    <summary>Paper</summary>
    Eckenrode, Kelly B., Dario Righelli, Marcel Ramos, Ricard Argelaguet, Christophe Vanderaa, Ludwig Geistlinger, Aedin C. Culhane, et al. “Curated Single Cell Multimodal Landmark Datasets for R/Bioconductor.” Edited by Mingyao Li. PLOS Computational Biology 19, no. 8 (August 25, 2023): e1011324. https://doi.org/10.1371/journal.pcbi.1011324.
</details>

- [Broad Single Cell Portal](https://singlecell.broadinstitute.org/single_cell) - scRNA-seq data from human, mouse, and more. Description, download

- [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/) - single cell data collections, datasets, gene expression. Downloadable. [CuratedAtlasQueryR](https://github.com/stemangiola/CuratedAtlasQueryR) - an R package for the programmatic exploration and retrieval of the harmonised, curated and reannotated CELLxGENE single-cell human cell atlas. Data can be retrieved at cell, sample, or dataset levels based on filtering criteria.

- [10X Genomics Visium spatial transcriptomics data](https://support.10xgenomics.com/spatial-gene-expression/datasets) - adult human/mouse brain, human breast cancer

- [Detailed description of 9 datasets and download links](https://github.com/JinmiaoChenLab/Batch-effect-removal-benchmarking), from Tran et al. "[A Benchmark of Batch-Effect Correction Methods for Single-Cell RNA Sequencing Data](https://doi.org/10.1186/s13059-019-1850-9)" study. Dataset 1: human dendritic cells, Dataset 2: mouse cell atlas, Dataset 3: simulation, Dataset 4: human pancreas, Dataset 5: human peripheral blood mononuclear cell (PBMC), Dataset 6: cell line, Dataset 7: mouse retina, Dataset 8: mouse brain, Dataset 9: human cell atlas

- [Single-cell portal, Broad Institute](https://portals.broadinstitute.org/single_cell). https://portals.broadinstitute.org/single_cell

- [Conquer DB of scRNA-seq datasets as R MultiAssayExperiment objects](http://imlspenticton.uzh.ch:3838/conquer/).

- [10X Genomics datasets](https://support.10xgenomics.com/single-cell-gene-expression/datasets)

- [Immunogenomics.io](https://immunogenomics.io/) - immunology and genomics data, visualization and download.

- [scRNAseq](https://bioconductor.org/packages/scRNAseq/) - Collection of Public Single-Cell RNA-Seq Datasets, to play with in R

- [muscData](https://bioconductor.org/packages/muscData/) - Multi-sample multi-group scRNA-seq data R package

- [UCSC cell browser](http://cells.ucsc.edu/) - web interface, Python pipeline and Javascript scatter plot library for single-cell datasets visualization. Static HTML/CSS files for visualizing via a local server. Scripts for data import (cbImportCellranger, scImportSeurat, cbImportScanpy, etc.), or for building from text files (expression matrix, cell metadata with cluster definition, cell coordinates). Many pre-build datasets are available. [Documentation](https://cellbrowser.readthedocs.io/en/master/), [GitHub](https://github.com/maximilianh/cellBrowser). <details>
    <summary>Paper</summary>
    Speir, Matthew L, Aparna Bhaduri, Nikolay S Markov, Pablo Moreno, Tomasz J Nowakowski, Irene Papatheodorou, Alex A Pollen, Lucas Seninge, W James Kent, and Maximilian Haeussler. "UCSC Cell Browser: Visualize Your Single-Cell Data" https://doi.org/10.1101/2020.10.30.361162  Preprint. Bioinformatics, October 31, 2020. 
</details>

- [PanglaoDB](https://panglaodb.se/) - a web server for exploration of mouse and human single-cell RNA sequencing data. 10X Genomics, Chromium, Smart-seq2 data processed using a unified computational pipeline. [Search](https://panglaodb.se/search.html) for gene, browse and visualize [studies](https://panglaodb.se/samples.html), explore cell type gene expression [markers](https://panglaodb.se/markers.html?cell_type=%27choose%27), raw and RData [data download](https://panglaodb.se/bulk.html). <details>
    <summary>Paper</summary>
    Franzén, Oscar, Li-Ming Gan, and Johan L M Björkegren. “PanglaoDB: A Web Server for Exploration of Mouse and Human Single-Cell RNA Sequencing Data” Database (Oxford), 2019 Jan 1,https://doi.org/10.1093/database/baz046
</details>

- [CellBench data](https://github.com/LuyiTian/CellBench_data), single cell RNA-seq benchmarking, R SingleCellExperiment object.

- [SingleCellExperiment data](https://hemberg-lab.github.io/scRNA.seq.datasets/), human and mouse, brain, embryo development, embryo stem cells, hematopoietic stem cells, pancreas, retina, other tissues. [GitHub](https://github.com/hemberg-lab/scRNA.seq.datasets)

- [12 scRNA-seq datasets](http://github.com/pcahan1/singleCellNet/) (Tabula Muris, Microwell-Seq, Baron pancreas, Xin pancreas, Segerstolpe pancreas, Murano pancreas, Zheng PBMC, Darminis brain, Zeisel brain, Tasic cortex) processed in the SingleCellNet study.

- [Pan-cancer portal](http://blueprint.lambrechtslab.org/) of genomic blueprint of stromal cell heterogeneity using scRNA-seq ( from lung, colorectal, ovary, and breast cancer and normal tissues.68 stromal cell populations, 46 are shared and 22 are unique. Clustering and characterization of each cluster - marker genes, functional enrichment, transcription factors (SCENIC), trajectory reconstruction (Monocle). [Table S5](https://www.nature.com/articles/s41422-020-0355-0#Sec26) -  markers of cells present in stroma. A web portal for exploring gene expression in each cell subtype, comparison of conditions, raw data is not available. Visualization using SCope. <details>
    <summary>Paper</summary>
    Qian, Junbin, Siel Olbrecht, Bram Boeckx, Hanne Vos, Damya Laoui, Emre Etlioglu, Els Wauters, et al. "A Pan-Cancer Blueprint of the Heterogeneous Tumor Microenvironment Revealed by Single-Cell Profiling" https://doi.org/10.1038/s41422-020-0355-0  Cell Research, June 19, 2020. 
</details>

- [FUMA](https://github.com/Kyoko-wtnb/FUMA_scRNA_data) - Functional Mapping and Annotation of GWAS using 43 scRNA-seq datasets (human and mouse). MAGMA Cell type-specific enrichment. [Applied to 26 GWAS disorders](https://fuma.ctglab.nl/). Processed data and instructions for self-download. <details>
    <summary>Paper</summary>
    Watanabe, Kyoko, Maša Umićević Mirkov, Christiaan A. de Leeuw, Martijn P. van den Heuvel, and Danielle Posthuma. "Genetic Mapping of Cell Type Specificity for Complex Traits" https://doi.org/10.1038/s41467-019-11181-1  Nature Communications 10, no. 1 (December 2019)
</details>

- [A list of scRNA-seq studies](http://www.nxn.se/single-cell-studies/). List of scRNA-seq databases (The Human Cell Atlas, JingleBells, conquer, PangaloDB, the EMBL-EBI Single Cell Expression Atlas, Single Cell Portal, scRNASeqDB). Number of cell types is directly proportional to the number of cells analyzed. <details>
    <summary>Paper</summary>
    Svensson, Valentine, and Eduardo da Veiga Beltrame. "A Curated Database Reveals Trends in Single Cell Transcriptomics" https://doi.org/10.1101/742304  Preprint. Genomics, August 21, 2019. 
</details>

- [SCPortalen](http://single-cell.clst.riken.jp/) - database of scRNA-seq datasets from the International Nucleotide Sequence Database Collaboration (INSDC). Manually curated datasets with metadata,QC's, processed, PCA and tSNE coordinates, FPKM gene expression. <details>
    <summary>Paper</summary>
    Abugessaisa, Imad, Shuhei Noguchi, Michael Böttcher, Akira Hasegawa, Tsukasa Kouno, Sachi Kato, Yuhki Tada, et al. "SCPortalen: Human and Mouse Single-Cell Centric Database" https://doi.org/10.1093/nar/gkx949  Nucleic Acids Research 46, no. D1 (04 2018)
</details>

- [scRNASeqDB](https://bioinfo.uth.edu/scrnaseqdb/) - a human-oriented scRNA-seq database. 38 studies, 200 cell types. <details>
    <summary>Paper</summary>
    Cao, Yuan, Junjie Zhu, Peilin Jia, and Zhongming Zhao. "ScRNASeqDB: A Database for RNA-Seq Based Gene Expression Profiles in Human Single Cells" https://doi.org/10.3390/genes8120368  Genes 8, no. 12 (December 5, 2017). 
</details>

- [scEiaD](https://plae.nei.nih.gov) - scRNA-seq data of the eye. 1.2 million single-cell back of the eye transcriptomes across 28 studies, 18 publications, and 3 species.



### Human

- [CuratedAtlasQueryR](https://stemangiola.github.io/CuratedAtlasQueryR/) - human immune system single-cell atlas (29M cells) across different organs (45 anatomical sites), age, sex, and ethnicity. Detailed analysis of age-, sex-, and ethnicity changes (Tregs are particularly affected), figures. [sccomp](https://github.com/MangiolaLaboratory/sccomp/tree/multilevel) - R package for compositional analysis of single-cell data. <details>
    <summary>Paper</summary>
    “A Multi-Organ Map of the Human Immune System across Age, Sex and Ethnicity,” bioRxiv, April 29, 2024. https://doi.org/10.1101/2023.06.08.542671
</details>

- [Single Cell Atlas (SCA)](https://www.singlecellatlas.org/) - data from five single-cell omics (scRNA-seq, scATAC-seq, scImmune, CyTOF, flow cytometry, [GitHub code for processing](https://github.com/eudoraleer/sca), spatial transcriptomics, and two bulk omics across 125 healthy adult and fetal tissues. Selected multi-omics datasets Strict QC. [Download](https://www.singlecellatlas.org/download) processed data. <details>
    <summary>Paper</summary>
    Pan, Lu, Paolo Parini, Roman Tremmel, Joseph Loscalzo, Volker M. Lauschke, Bradley A. Maron, Paola Paci, et al. “Single Cell Atlas: A Single-Cell Multi-Omics Human Cell Encyclopedia.” Genome Biology 25, no. 1 (April 19, 2024): 104. https://doi.org/10.1186/s13059-024-03246-2.
</details>

- **Multi-tissue human single-cell atlas**, 209,126 snRNA-seq profiles, 8 tissue types, annotated with 43 broad and 74 fine categories. Lipid-assiciated (LAM) macrophages. Raw sequence data are available at [AnVIL](https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_GTEx_V9_hg38); dbGaP accession phs000424. Gene expression matrices are available from the [GTEx Portal](www.gtexportal.org) and the [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP1479). The Jupyter notebooks are available at [Zenodo](https://zenodo.org/records/5775379). <details>
    <summary>Paper</summary>
    Eraslan, Gökcen, Eugene Drokhlyansky, Shankara Anand, Evgenij Fiskin, Ayshwarya Subramanian, Michal Slyper, Jiali Wang, et al. “Single-Nucleus Cross-Tissue Molecular Reference Maps toward Understanding Disease Gene Function.” Science 376, no. 6594 (May 13, 2022): eabl4290. https://doi.org/10.1126/science.abl4290.
</details>

- [HLCA](https://beta.fastgenomics.org/p/hlca) - Human Lunc Cell Atlas. 46 datasets, over 2.2 million cells from 44 individuals. integrated with best performing scANVI. [scRNAseq_pipelines](https://github.com/LungCellAtlas/scRNAseq_pipelines) - the HLCA pipeline for processing of sequencing data to count matrices. [HLCA_reproducibility](https://github.com/LungCellAtlas/HLCA_reproducibility) - code used for the HLCA project. [mapping_data_to_the_HLCA](https://github.com/LungCellAtlas/mapping_data_to_the_HLCA) - code for users to map new data to the HLCA core. [HLCA reference model, embedding and extended models](https://zenodo.org/record/6337966#.Yid5Vi9Q28U). [Data](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293). <details>
    <summary>Paper</summary>
    Sikkema, L, D Strobl, L Zappia, E Madissoon, NS Markov, L Zaragosi, M Ansari, et al. “An Integrated Cell Atlas of the Human Lung in Health and Disease,” March 11, 2022. https://doi.org/10.1101/2022.03.10.483747
</details>

- [Tablua Sapiens](https://tabula-sapiens-portal.ds.czbiohub.org/) - human scRNA-seq reference database. 500,000 cells from 24 organs of 15 normal human subjects. Visualization and exploratory analysis of all and organ-specific datasets. scVI low-dimensional representation, UMAP visualization, cellxgene for cell annotation, SICILIAN to identify splice junctions. [GitHub](https://github.com/czbiohub/tabula-sapiens), Raw data on [Registry of Open Data on AWS](https://registry.opendata.aws/tabula-sapiens/). Processed data on [FigShare](https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219) and [GEO GSE149590](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149590). <details>
    <summary>Paper</summary>
    Quake, Stephen R, Tabula Sapiens Consortium, and others. "The Tabula Sapiens: A Single Cell Transcriptomic Atlas of Multiple Organs from Individual Human Donors" https://doi.org/10.1101/2021.07.19.452956  Biorxiv, 2021.
</details>

- [Human Cell Atlas Preview Datasets](https://preview.data.humancellatlas.org/). Massive amount of data. [HCAData](https://github.com/federicomarini/HCAData) - an R package to access this data is being developed. [Galaxy interface to HCA data analysis](https://humancellatlas.usegalaxy.eu/)

- [Single Cell Immune Profiling Datasets](https://support.10xgenomics.com/single-cell-vdj/datasets). \~150,000 CD8+ T cells from 4 human donors across a highly multiplexed panel of 44 distinct, specific peptide–MHC (pMHC) multimers. 

- Normal lung scRNA-seq data, human and mouse. 58 cell populations, 41 known and 14 unknown. Unknown cell types are replicated in humans, absent in mice. 10X Genomics V2, CellRanger, Seurat. [Counts/UMI tables on Synapse](https://www.synapse.org/#!Synapse:syn21041850/wiki/600865), [online browser](https://hlca.ds.czbiohub.org/), [code and data GitHub repo](https://github.com/krasnowlab/HLCA). <details>
    <summary>Paper</summary>
    Travaglini, Kyle J, Ahmad N Nabhan, Lolita Penland, Rahul Sinha, Astrid Gillich, Rene V Sit, Stephen Chang, et al. "A Molecular Cell Atlas of the Human Lung from Single-Cell RNA Sequencing" https://doi.org/10.1038/s41586-020-2922-4  Nature 587, no. 7835 (November 18, 2020): 619–25.
    [Supplementary Table 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2922-4/MediaObjects/41586_2020_2922_MOESM3_ESM.xlsx) - canonical gene markers of 45 lung cell types
    [Supplementary Table 4](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2922-4/MediaObjects/41586_2020_2922_MOESM6_ESM.xlsx) - Markers of all cell populations, for each cluster
</details>

- [14,039 human PBMCs](https://www.nature.com/articles/nbt.4096#supplementary-information) from eight patients into two groups: one stimulated with interferon-beta (IFN-β) and a culture-matched control. Eight clusters. [Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583). <details>
    <summary>Paper</summary>
    Kang, Hyun Min, Meena Subramaniam, Sasha Targ, Michelle Nguyen, Lenka Maliskova, Elizabeth McCarthy, Eunice Wan, et al. "Multiplexed Droplet Single-Cell RNA-Sequencing Using Natural Genetic Variation" https://doi.org/10.1038/nbt.4042  Nature Biotechnology 36, no. 1 (January 2018)
</details>

- [Nine PBMC datasets](https://genome.med.nyu.edu/results/external/iCellR/data/) provided by the Broad institute to test batch effect, data in text and .rda formats. <details>
    <summary>Paper</summary>
    Khodadadi-Jamayran, Alireza, Joseph Pucella, Hua Zhou, Nicole Doudican, John Carucci, Adriana Heguy, Boris Reizis, and Aristotelis Tsirigos. "ICellR: Combined Coverage Correction and Principal Component Alignment for Batch Alignment in Single-Cell Sequencing Analysis" https://doi.org/10.1101/2020.03.31.019109  Preprint. Bioinformatics, April 1, 2020
</details>

- Paul, Franziska, Ya’ara Arkin, Amir Giladi, Diego Adhemar Jaitin, Ephraim Kenigsberg, Hadas Keren-Shaul, Deborah Winter, et al. "Transcriptional Heterogeneity and Lineage Commitment in Myeloid Progenitors" https://doi.org/10.1016/j.cell.2015.11.013  Cell 163, no. 7 (December 17, 2015) - 3′ MARS (massively parallel RNA single-cell)-Seq protocol with shallow sequencing (1,453 genes/cell) to examine heterogeneity across bone marrow resident myeloid progenitors, 2,686 cells. [Batch-corrected UMI count matrix: 19 clusters](http://compgenomics.weizmann.ac.il/tanay/?page_id=649). Seurat (PMID: 29608179) clusters, [Supplementary Data 2](https://www.nature.com/articles/nbt.4096#supplementary-information). [GEO GSE72857](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72857). <details>
    <summary>Paper</summary>
    [Amit_2015/Supplementary_dataset_1.txt] - Cell metadata for IFNB response analysis. [Source](https://media.nature.com/original/nature-assets/nbt/journal/v36/n5/extref/nbt.4096-S3.txt)
    [Amit_2015/Supplementary_dataset_2.txt] - Cell metadata for murine hematopoiesis analysis. [Source](https://media.nature.com/original/nature-assets/nbt/journal/v36/n5/extref/nbt.4096-S4.txt)
    [Amit_2015] - Cell metadata for cross-species pancreatic islet analysis. [Source](https://media.nature.com/original/nature-assets/nbt/journal/v36/n5/extref/nbt.4096-S5.txt)
    [Amit_2015/generate_clustering.tar.gz] - all the data and scripts required for reproducing the MARS-seq results. 913Gb. [Source](http://www.wisdom.weizmann.ac.il/~arkiny/Paul2015/generate_clustering.tar.gz)
</details>

- Yan, Liying, Mingyu Yang, Hongshan Guo, Lu Yang, Jun Wu, Rong Li, Ping Liu, et al. "[Single-Cell RNA-Seq Profiling of Human Preimplantation Embryos and Embryonic Stem Cells](https://doi.org/10.1038/nsmb.2660)"  Nature Structural & Molecular Biology 20, no. 9 (September 2013) - scRNA-seq of human embryos and embryonic stem cells. 124 cells at different stages of development. Clustering, PCA. Expression of novel genes and de novo assembly of other transcripts, including lncRNAs. [Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36552).

#### Cancer

- [3CA Curated Cancer Cell Atlas](https://www.weizmann.ac.il/sites/3CA/) - Collected, annotated and analyzed cancer scRNA-seq datasets. Expression matrix, cell/gene names, copy number alterations, UMAP coordinates. [Tweet](https://twitter.com/TiroshLab/status/1473366146665623557?s=20&t=9axW_HQgFoZLXejUwjo7lg)

- [CancerSEA](http://biocc.hrbmu.edu.cn/CancerSEA/home.jsp) - cancer scRNA-seq studies. Download individual studies, as well as gene signatures (from Angiogenesis, DNA damage to EMT, metastasis, etc.)

- [scTIME Portal](http://sctime.sklehabc.com/unicellular/home) - a database and an exploration/analysis portal for single cell transcriptomes of tumor immune microenvironment. Cell clusters, expression of selected genes, data/image download. Links to other portals/databases.

- Short-term TGFBR inhibition (mimicking decrease in TGFB pathway activity after pregnancy, p27+ progenitors) of prepubertal ACI inbred and Sprague Dawley outbred rats prevents estrogen- and carcinogen-induced breast cancer by likely increasing and depleting the pool of epithelial subpopulation of secretory basal cells (SBCs) with progenitor features. Bulk, scRNA-seq (10X Genomics, CellRanger, Seurat, monocle3, liana, more methods), SBC signature (336 genes, [Supplementary Data 7](https://www.nature.com/articles/s41467-022-35043-5#Sec27)), reanalysis of [GSE106273](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106273) and [GSE161529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529). [GitHub](https://github.com/csimona/tumor-prevention-rat-scRNAseq), scRNA-seq data at [GSE184095](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/geo/query/acc.cgi?acc=GSE184095), R data object at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7293642.svg)](https://doi.org/10.5281/zenodo.7293642). Great [Tweetorial](https://twitter.com/simocristea/status/1600578512578035733?s=20&t=CgDFM9JAFhVjliIps8QucQ) by [Simona Cristea](https://twitter.com/simocristea). <details>
    <summary>Paper</summary>
    Alečković, Maša, Simona Cristea, Carlos R. Gil Del Alcazar, Pengze Yan, Lina Ding, Ethan D. Krop, Nicholas W. Harper, et al. “Breast Cancer Prevention by Short-Term Inhibition of TGFβ Signaling.” Nature Communications 13, no. 1 (December 7, 2022): 7558. https://doi.org/10.1038/s41467-022-35043-5.
    </details>

- [ScBrAtlas](#scbratlas) scRNA-seq data (R data objects and code) of normal and tumorigenic (ER+, HER2+, and TNBC) human breast tissue. Approx. 430K cells from 69 tissue specimens from 55 patients. 10X genomics, CellRanger, hg38, Seurat processing (custom settings), R objects. Differential expression/abundance/copy number variation analysis, KEGG enrichment analysis (limma::kegga), gene signature of luminal progenitor, mature luminal and stromal cell populations. Raw data at [GSE161529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529). [FigShare](https://figshare.com/s/c584dda937d346cc9a80), [GitHub](https://github.com/yunshun/HumanBreast10X). <details>
    <summary>Paper</summary>
    Chen, Yunshun, Bhupinder Pal, Geoffrey J. Lindeman, Jane E. Visvader, and Gordon K. Smyth. “R Code and Downstream Analysis Objects for the ScRNA-Seq Atlas of Normal and Tumorigenic Human Breast Tissue.” Scientific Data 9, no. 1 (December 2022): 96. https://doi.org/10.1038/s41597-022-01236-2.
</details>

- Multi-omics single-cell analysis of breast cancer. >130K scRNA-seq across 11 ER+, 5 HER2+ and 10 TNBC primary breast tumors. Immunophenotyping by CIRE-seq. 10X Visium Spatial transcriptomics. SCSubtype signatures - subtype classification (Basal, Her2E, LumA, LumB), Supplementary Table 4. Recurrent gene modules (GMs) driving neoplastic cell heterogeneity. Supplementary Table 5 - gene lists for 7 GMs. DScore (BIRC5, CCNB1, CDC20, NUF2, CEP55, NDC80, MKI67, PTTG1, RRM2, TYMS and UBE2C) and proliferation score. The cytotoxic gene list containing effector cytotoxic proteins (GZMA, GZMB, GZMH, GZMK, GZMM, GNLY, PRF1 and FASLG) and cytotoxic T cell activation markers (IFNG, TNF, IL2R and IL2). Bioinformatics methods: inferCNV, Stereoscope, CIBERSORTx, Monocle 2, CITE-seq-Count, DWLS. Code for all tools on [GitHub](https://github.com/Swarbricklab-code/BrCa_cell_atlas). [Supplementary Tables](https://www.nature.com/articles/s41588-021-00911-1#Sec39), [Processed data](https://singlecell.broadinstitute.org/single_cell/study/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers), [GEO GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078), [Spatially resolved transcriptomics data](https://zenodo.org/record/4739739). <details>
    <summary>Paper</summary>
    Wu, S.Z., Al-Eryani, G., Roden, D.L. et al. "A Single-Cell and Spatially Resolved Atlas of Human Breast Cancers" https://doi.org/10.1038/s41588-021-00911-1  Nature Genetics 53 (September 2021)
</details>

- [scRNA-seq of healthy breast, GEO GSE164898](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164898). Five freshly collected samples, 18K cells, 20K genes. 13 epithelial cell clusters. Breast cancers may originate from 3 luminal mature and 1 progenitor subclusters. TBX3 and PDK4 subclassify ER+ breast cancers in at least four subtypes. [Table S2](https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(21)00035-5#supplementaryMaterial) - cluster-specific gene expression. [Matrices in HDF5 format](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164898). <details>
    <summary>Paper</summary>
    Bhat-Nakshatri, Poornima, Hongyu Gao, Liu Sheng, Patrick C. McGuire, Xiaoling Xuei, Jun Wan, Yunlong Liu et al. "[A single-cell atlas of the healthy breast tissues reveals clinically relevant clusters of breast epithelial cells" https://doi.org/10.1016/j.xcrm.2021.100219   Cell Reports Medicine 2, no. 3 (2021): 100219.
</details>

- [Interactive mammary cell gene expression atlas](https://cells.ucsc.edu/?ds=mouse-mammary-epithelium-integrated) - Integrated 50K mouse and 24K human mammary epithelial cell atlases, scRNA-seq. Consensus lineage trajectory - embryonic stem cells differentiate into three epithelial lineages (Basal, luminal hormone-sensing L-Hor, luminal alveolar L-Alv). Integration of four public and one new datasets. Harmony, LIGER, scALIGN for integration. STREAM for lineage tracing. ssGSVA for gene set enrichment. [Supplementary Data](https://www.nature.com/articles/s42003-021-02201-2#Sec35): Supplementary Data 4 - mouse gene signatures of MaSC (mammary stem cells), Basal, LA-Pro (Luminal Alveloar progenitors), L-Alv, LH-Pro (Luminal Hormone-sensing), L-Hor. Supplementary Data 10 - mouse/human-specific and common stem/basal/Alv/Hor lineage genes. <details>
    <summary>Paper</summary>
    Saeki, Kohei, Gregory Chang, Noriko Kanaya, Xiwei Wu, Jinhui Wang, Lauren Bernal, Desiree Ha, Susan L. Neuhausen, and Shiuan Chen. "Mammary Cell Gene Expression Atlas Links Epithelial Cell Remodeling Events to Breast Carcinogenesis" https://doi.org/10.1038/s42003-021-02201-2  Communications Biology, (December 2021)
</details>

- <a name="scbratlas">ScBrAtlas</a> scRNA-seq catalog of breast cancer. >340,000 cells, normal breast, preneoplastic tissue, the major breast cancer subtypes, and pairs of tumors and involved lymph nodes. 34 treatment-naive primary tumors. Transition from preneoplastic to tumor involves immune microenvironment shift. ER+ tumors are different. 10X Genomics, hg39 alignment with CellRanger, integration with Seurat, classification by cell cycle markers, PAM50, immune cell-specific and other signatures. edgeR::read10X to read files, TMM normalization, limma-voom and TREAT for differential analysis. [InferCNV](https://github.com/broadinstitute/inferCNV) to infer CNVs from scRNA-seq data. [GEO GSE161529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529) - scRNA-seq data in MatrixMarket format (edgeR::read10X), 69 samples. [GEO GSE161892](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161892) - Bulk RNA-seq (luminal progenitors (LP), mature luminal (M), basal, stromal populations). [FigShare](https://figshare.com/s/c584dda937d346cc9a80) - R code and downstream analysis objects. <details>
    <summary>Paper</summary>
    Pal, Bhupinder, Yunshun Chen, François Vaillant, Bianca D Capaldo, Rachel Joyce, Xiaoyu Song, Vanessa L Bryant, et al. "A Single‐cell RNA Expression Atlas of Normal, Preneoplastic and Tumorigenic States in the Human Breast" https://doi.org/10.15252/embj.2020107333  The EMBO Journal, May 5, 2021

    Chen, Y. & Smyth, G. K. Data, [R code and output Seurat objects for single cell RNA-seq analysis of human breast tissues" https://doi.org/10.1101/2021.11.30.470523). https://figshare.com/s/c584dda937d346cc9a80
</details>

- scRNA-seq of breast cancer, four women, human mammary epithelial cells. Marker-free algorithm ([LandSCENT](https://github.com/ChenWeiyan/LandSCENT)) that identifies stem-like bipotent state, characterized by YBX1 and ENO1, two modulators of breast cancer risk. [Source data 6B](https://www.nature.com/articles/s42003-019-0554-8#Sec22) - 12- and 72 gene signature of bipotent state, basal-like. [GEO GSE113197](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113197) - scRNA-seq data, annotated. <details>
    <summary>Paper</summary>
    Chen, Weiyan, Samuel J. Morabito, Kai Kessenbrock, Tariq Enver, Kerstin B. Meyer, and Andrew E. Teschendorff. "Single-Cell Landscape in Mammary Epithelium Reveals Bipotent-like Cells Associated with Breast Cancer Risk and Outcome" https://doi.org/10.1038/s42003-019-0554-8  Communications Biology, (December 2019)
</details>

- scRNA-seq of >25K normal human breast epithelial cells from seven individuals. Three cell populations, one basal and two luminal (secretory L1 and hormone-responsive L2). Within luminal L1, three cell states (milk production, secretory, epithelial keratin expression), but the combined analysis reports L1_1 and L1_2 signatures. Fluidigm, 10X Genomics, Seurat, Monocle analyses. [GitHub](https://github.com/kessenbrocklab/Nguyen_Pervolarakis_Nat_Comm_2018), [GEO GSE113197](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113197). [Supplementary Data 2](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-04334-1/MediaObjects/41467_2018_4334_MOESM5_ESM.xlsx) - myeloepithelial gene signature to stratify basal cells into either "Basal" or "Myeloepithelial" grouping; scRNA-seq derived "Basal", "Basal_myoepithelial", "L1_1", "L1_2", "L_2", "Unclassified" signatures; Metabric-derived LumA and LumB signatures. <details>
    <summary>Paper</summary>
    Nguyen, Quy H., Nicholas Pervolarakis, Kerrigan Blake, Dennis Ma, Ryan Tevia Davis, Nathan James, Anh T. Phung, et al. "Profiling Human Breast Epithelial Cells Using Single Cell RNA Sequencing Identifies Cell Diversity" https://doi.org/10.1038/s41467-018-04334-1  Nature Communications 9, no. 1 (December 2018)
</details>

- scRNA-seq of immune cells in BRCA - continuous activation of T cells, no macrophage polarization. inDrop and 10X platforms. 47,016 CD45+ cells from 8 primary breast carcinomas. 83 clusters, tested by cross-validation. [GEO GSE114727](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114727), [GEO GSE114724](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114725), [GEO GSE114724](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114724). <details>
    <summary>Paper</summary>
    Azizi, Elham, Ambrose J. Carr, George Plitas, Andrew E. Cornish, Catherine Konopacki, Sandhya Prabhakaran, Juozas Nainys, et al. "Single-Cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment" https://doi.org/10.1016/j.cell.2018.05.060  Cell, June 2018. 
</details>

### Mouse

- [Tablua Muris](http://tabula-muris.ds.czbiohub.org/) - Mouse scRNA-seq of 100,605 cells from 20 organs and tissues. http://tabula-muris.ds.czbiohub.org/, [GitHub with gene matrix download](https://github.com/czbiohub/tabula-muris). [Raw data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109774). [Processed](https://figshare.com/articles/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells_v2_/5829687/1). <details>
    <summary>Paper</summary>
    Quake, Stephen R., Tony Wyss-Coray, and Spyros Darmanis. "Single-Cell Transcriptomic Characterization of 20 Organs and Tissues from Individual Mice Creates a Tabula Muris" https://doi.org/10.1101/237446  March 29, 2018. 
</details>

- [10X Genomics, The 1 million neuron data set from E18 mice](https://support.10xgenomics.com/single-cell-gene-expression/datasets). R packages for its analyses: 
    - [TENxGenomics](https://github.com/mtmorgan/TENxGenomics) - https://github.com/mtmorgan/TENxGenomics
    - [TENxBrainAnalysis](https://github.com/Bioconductor/TENxBrainAnalysis) - https://github.com/Bioconductor/TENxBrainAnalysis

- [DropViz](http://dropviz.org/) - Exploring the Mouse Brain through Single Cell Expression Profiles. Drop-seq to analyze 690,000 individual cells from nine different regions of the adult mouse brain. http://dropviz.org/

- [STARmap](https://www.starmapresources.com/data/) - in situ gene expression datasets of the mouse visual cortex (890 cells, 1,020 genes).

- [Mouse Cell Atlas](http://bis.zju.edu.cn/MCA/)

- scRNA-seq and scATAC-seq of normal mammary epithelial cells (MECs, mouse). 4 main clusters, their characteristics. Trajectory analysis, regulatory modules and TFs. Seurat/Signac, Monocle, Cicero, cisTopic, ChromVar, Homer. Processed data: [GSE157890](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157890). <details>
    <summary>Paper</summary>
    Pervolarakis, Nicholas, Quy H. Nguyen, Justice Williams, Yanwen Gong, Guadalupe Gutierrez, Peng Sun, Darisha Jhutty, et al. “Integrated Single-Cell Transcriptomics and Chromatin Accessibility Analysis Reveals Regulators of Mammary Epithelial Cell Identity.” Cell Reports 33, no. 3 (October 2020): 108273. https://doi.org/10.1016/j.celrep.2020.108273.
</details>

- [scRNA-seq of mouse gastrulation and early embryogenesis](https://github.com/MarioniLab/MouseGastrulationData).  116,312 cells, nine time points from 6.5 to 8.5 days post-fertilization. Methods, scran-based analysis, [tutorial](https://marionilab.cruk.cam.ac.uk/MouseGastrulation2018/), [download script](https://github.com/MarioniLab/EmbryoTimecourse2018), [Bioconductor data package](https://github.com/MarioniLab/MouseGastrulationData). <details>
    <summary>Paper</summary>
    Pijuan-Sala, Blanca, Jonathan A. Griffiths, Carolina Guibentif, Tom W. Hiscock, Wajid Jawaid, Fernando J. Calero-Nieto, Carla Mulas, et al. "A Single-Cell Molecular Map of Mouse Gastrulation and Early Organogenesis" https://doi.org/10.1038/s41586-019-0933-9  Nature, February 20, 2019.
</details>

- [scRNA-seq of aging](https://mca.research.calicolabs.com/). >50,000 cells from kidney, lung, and spleen in young (7 months) and aged (22-23 months) mice. Transcriptional variation using difference from the median. Cell-cell heterogeneity using the Euclidean distance from centroids. Aging trajectories derived from NMF embedding. Cell type identification by neural network trained on Tabula Muris. ln(CPM + 1) UMIs used for all analyses. Visualization, downloadable data, code, pre-trained network. <details>
    <summary>Paper</summary>
    Kimmel, Jacob C, Lolita Penland, Nimrod D Rubinstein, David G Hendrickson, David R Kelley, and Adam Z Rosenthal. "A Murine Aging Cell Atlas Reveals Cell Identity and Tissue-Specific Trajectories of Aging" https://doi.org/10.1101/657726  BioRxiv, January 1, 2019
</details>

- [Mouse spinal cord development scRNA-seq](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7320/). Temporal (embryonic day 9.5-13.5) and spatial (cervical and thoracic regions of the neural tube) profiling. 10X genomics protocol, Cell Ranger processing, filtering, combinatorial testing for differential expression, pseudotime reconstruction using Monocle2. UMI matrix (21465 cells), https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7320/, [code](https://github.com/juliendelile/MouseSpinalCordAtlas), [Table S1](http://www.biologists.com/DEV_Movies/DEV173807/TableS1.csv) - binarized matrix of gene markers for each neuronal subtype. <details>
    <summary>Paper</summary>
    Delile, Julien, Teresa Rayon, Manuela Melchionda, Amelia Edwards, James Briscoe, and Andreas Sagner. "Single Cell Transcriptomics Reveals Spatial and Temporal Dynamics of Gene Expression in the Developing Mouse Spinal Cord" https://doi.org/10.1242/dev.173807  Development, March 7, 2019
</details>

- scRNA-seq of mouse testes (2,500 cells from two 8-week-old C57Bl/6J mice, 10X, CellRanger 2.0). [GitHub](https://github.com/slukassen/SCS_testis/) with analysis scripts. Raw data at [GSE104556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104556), processed data at [FigShare](http://dx.doi.org/10.6084/m9.figshare.c.4119713). <details>
    <summary>Paper</summary>
    Lukassen, Soeren, Elisabeth Bosch, Arif B. Ekici, and Andreas Winterpacht. “Single-Cell RNA Sequencing of Adult Mouse Testes.” Scientific Data 5 (September 11, 2018): 180192. https://doi.org/10.1038/sdata.2018.192.
</details>

- [Sci-CAR, single-cell RNA- and ATAC-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117089). Two experiments: 1) Lung adenocarcinoma A549 cells, dexametasone treatment over 3 timepoints. 2) Mixture of HEK293T (human) and NIH3T3 (mouse) cells. Differential gene expression, accessibility analysis, clustering. Linking distal open chromatin to genes, 44% map to nearest, 21 to the second nearest. Gene expression counts and ATAC-seq peaks, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117089
    Cao, Junyue, Darren A. Cusanovich, Vijay Ramani, Delasa Aghamirzaie, Hannah A. Pliner, Andrew J. Hill, Riza M. Daza, et al. "Joint Profiling of Chromatin Accessibility and Gene Expression in Thousands of Single Cells" https://doi.org/10.1126/science.aau0730  Science 361, no. 6409 (September 28, 2018)

- [scRNA-seq (10X Genomics) of murine cerebellum](https://cellseek.stjude.org/cerebellum/). 39245 cells (after filtering) from 12 time points, 48 distinct clusters. [Cell Seek](https://cellseek.stjude.org/cerebellum/) for online exploration, and [GitHub]). [Data](https://www.ebi.ac.uk/ena/data/view/PRJEB23051) - Approx. 83,000 cells (BAM files per time point). No annotations. <details>
    <summary>Paper</summary>
    Carter, Robert A., Laure Bihannic, Celeste Rosencrance, Jennifer L. Hadley, Yiai Tong, Timothy N. Phoenix, Sivaraman Natarajan, John Easton, Paul A. Northcott, and Charles Gawad. "A Single-Cell Transcriptional Atlas of the Developing Murine Cerebellum" https://doi.org/10.1016/j.cub.2018.07.062  Current Biology, September 2018. 
</details>

- [Single-cell ATAC-seq](http://atlas.gs.washington.edu/mouse-atac/), approx. 100,000 single cells from 13 adult mouse tissues. Two sequence platforms, good concordance. Filtered data assigned into 85 clusters. Genes associated with the corresponding ATAC sites (Cicero for identification). Differential accessibility. Motif enrichment (Basset CNN). GWAS results enrichment. All data and metadata are available for download as text or rds format. <details>
    <summary>Paper</summary>
    Cusanovich, Darren A., Andrew J. Hill, Delasa Aghamirzaie, Riza M. Daza, Hannah A. Pliner, Joel B. Berletch, Galina N. Filippova, et al. "A Single-Cell Atlas of In Vivo Mammalian Chromatin Accessibility" https://doi.org/10.1016/j.cell.2018.06.052  Cell 174, no. 5 (August 2018)
</details>

- [Mouse pancreas scRNA-seq, 12,000 cells](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133). Deconvolution of bulk RNA-seq data using cell signatures derived from clusters of scRNA-seq cells, [GitHub](https://github.com/shenorrLab/bseqsc). <details>
    <summary>Paper</summary>
    Baron, Maayan, Adrian Veres, Samuel L. Wolock, Aubrey L. Faust, Renaud Gaujoux, Amedeo Vetere, Jennifer Hyoje Ryu, et al. "A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-Cell Population Structure" https://doi.org/10.1016/j.cels.2016.08.011  Cell Systems 3, no. 4 (26 2016)
</details>

- [The full-length SMART-Seq2 protocol with deep sequencing (6,558 genes/cell) to profile 765 multipotent mouse hematopoietic progenitors](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81682). <details>
    <summary>Paper</summary>
    Nestorowa, S., F. K. Hamey, B. Pijuan Sala, E. Diamanti, M. Shepherd, E. Laurenti, N. K. Wilson, D. G. Kent, and B. Gottgens. "A Single-Cell Resolution Map of Mouse Hematopoietic Stem and Progenitor Cell Differentiation" https://doi.org/10.1182/blood-2016-05-716480  Blood 128, no. 8 (August 25, 2016)
</details>

- Buettner, Florian, Kedar N Natarajan, F Paolo Casale, Valentina Proserpio, Antonio Scialdone, Fabian J Theis, Sarah A Teichmann, John C Marioni, and Oliver Stegle. "[Computational Analysis of Cell-to-Cell Heterogeneity in Single-Cell RNA-Sequencing Data Reveals Hidden Subpopulations of Cells](https://doi.org/10.1038/nbt.3102)"  Nature Biotechnology 33, no. 2 (March 2015).
    - [data/scLVM/nbt.3102-S7.xlsx] - Uncorrected and cell-cycle corrected expression values (81 cells x 7073 genes) for T-cell data. Includes cluster assignment to naive T cells vs. TH2 cells (GATA3 high marker). [Source](https://media.nature.com/original/nature-assets/nbt/journal/v33/n2/extref/nbt.3102-S7.xlsx)
    - [data/scLVM/nbt.3102-S8.xlsx] - Corrected and uncorrected expression values for the newly generated mouse ESC data. 182 samples x 9571 genes. [Source](https://media.nature.com/original/nature-assets/nbt/journal/v33/n2/extref/nbt.3102-S8.xlsx)

- Zeisel, A., Munoz-Manchado, A.B., Codeluppi, S., Lonnerberg, P., La Manno, G., Jureus, A., Marques, S., Munguba, H., He, L., Betsholtz, C., et al. (2015). Brain structure. [Cell types in the mouse cortex and hippocampus revealed by single-cell RNA- seq](https://science.sciencemag.org/content/347/6226/1138). Science 347, 1138–1142. - 3,005 single cells from the hippocampus and cerebral cortex of mice. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361), [Web](http://linnarssonlab.org/cortex/), and more on this site.
    - [data/Brain/Zeisel_2015_TableS1.xlsx] - Table S1 - gene signatures for Ependymal, Oligodendrocyte, Microglia, CA1 Pyramidal, Interneuron, Endothelial, S1 Pyramidal, Astrocyte, Mural cells. [Source](http://science.sciencemag.org/highwire/filestream/628248/field_highwire_adjunct_files/1/aaa1934_TableS1.xlsx)
    - [data/Brain/expression_mRNA_17-Aug-2014.txt] - 19,972 genes x 3005 cells. Additional rows with class annotations to interneurons, pyramidal SS, pyramidal CA1, oligodendrocytes, microglia, endothelial-mural, astrocytes_ependymal, further subdivided into 47 subclasses. [Source](https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt)



### Brain

- [HNOCA](https://github.com/theislab/neural_organoid_atlas) - human neural organoid cell atlas, scRNA-seq, over 1.7M cells. Three-step integration pipeline. <details>
    <summary>Paper</summary>
    He, Zhisong, Leander Dony, Jonas Simon Fleck, Artur Szalata, Katelyn X. Li, Irena Sliskovic, Hsiu-Chuan Lin, et al. “An Integrated Transcriptomic Cell Atlas of Human Neural Organoids.” Preprint. Developmental Biology, October 6, 2023. https://doi.org/10.1101/2023.10.05.561097.
</details>

- [Spatially resolved brain scRNA-seq isoform atlas](https://isoformatlas.com/). ScISOrSeq (single-cell isoform sequencing), long-read sequencing (ONT, PacBio), Visium Hippocampus and prefrontal cortex. [GEO GSE158450](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158450) - scRNA-seq data. [GitHub](https://github.com/noush-joglekar/scisorseqr). <details>
    <summary>Paper</summary>
    Joglekar, Anoushka. “A Spatially Resolved Brain Region- and Cell Type-Specific Isoform Atlas of the Postnatal Mouse Brain,” Nature Communications, 19 January 2021, https://doi.org/10.1038/s41467-020-20343-5
</details>

- [scRNA-seq and scATAC-seq integration, human forebrain development timecourse](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132403). Analysis, reporting, data. <details>
    <summary>Paper</summary>
    Trevino, Alexandro E., Nasa Sinnott-Armstrong, Jimena Andersen, Se-Jin Yoon, Nina Huber, Jonathan K. Pritchard, Howard Y. Chang, William J. Greenleaf, and Sergiu P. Pașca. "Chromatin Accessibility Dynamics in a Model of Human Forebrain Development" https://doi.org/10.1126/science.aay1645  Science 367, no. 6476 (January 24, 2020)
</details>

- scRNA-seq of mouse cerebellar cortex, cell types. Over 780K single nuclei from 16 cortex lobules. LIGER for data integration and clustering. Seurat processing. Gene expression changes are continuous across cell types. Data are available through [NeMO](https://portal.nemoarchive.org/) and [Single-cell portal](https://singlecell.broadinstitute.org/single_cell/study/SCP795/a-transcriptomic-atlas-of-the-mouse-cerebellum). Scripts on [GitHub](https://github.com/MacoskoLab/cerebellum-atlas-analysis). <details>
    <summary>Paper</summary>
    Kozareva, Velina. "A Transcriptomic Atlas of Mouse Cerebellar Cortex Comprehensively Defines Cell Types" https://doi.org/10.1038/s41586-021-03220-z  Nature, 06 October 2021
</details>

- scRNA-seq of aging mouse brain (>50K single cells, 2-3 month and 21-23 month old mice). 37K cells after filtering, clustered into 25 cell types Age-related genes, (bidirectional) signatures, pathways and ligand-receptor interactions. [Supplementary material](https://www.nature.com/articles/s41593-019-0491-3#Sec32)
    - [Supplementary Table 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0491-3/MediaObjects/41593_2019_491_MOESM3_ESM.xlsx) - List of abbreviations for all cell types and subtypes.
    - [Supplementary Table 2](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0491-3/MediaObjects/41593_2019_491_MOESM4_ESM.xlsx) - Transcriptional profile (TPMs) of each cell type and subtype by age.
    - [Supplementary Table 3](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0491-3/MediaObjects/41593_2019_491_MOESM5_ESM.xlsx) - List of most discriminating genes per cell type.
    - [Supplementary Table 4](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0491-3/MediaObjects/41593_2019_491_MOESM6_ESM.xlsx) - List of most discriminating genes per neuronal subtype.
    - [Supplementary Table 5](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0491-3/MediaObjects/41593_2019_491_MOESM7_ESM.xlsx) - Cell counts of each cell type and subtype by age and animal.
    - [Supplementary Table 6](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0491-3/MediaObjects/41593_2019_491_MOESM8_ESM.xlsx) - Differential gene expression data between young and old cell types.
    - [Supplementary Table 7](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0491-3/MediaObjects/41593_2019_491_MOESM9_ESM.xlsx) - Differential gene expression data between young and old neuronal subtypes.
    - [Supplementary Table 8](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0491-3/MediaObjects/41593_2019_491_MOESM10_ESM.xlsx) - Matrix of aging-related genes across cell types and neuronal subtypes.
    - [Supplementary Table 9](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0491-3/MediaObjects/41593_2019_491_MOESM11_ESM.xlsx) - List of aging-related pathways and processes by cell type and neuronal subtype.
    - [Supplementary Table 10](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-019-0491-3/MediaObjects/41593_2019_491_MOESM12_ESM.xlsx) - Matrix of aging-related pathways and processes across cell types and neuronal subtypes.
    Ximerakis, Methodios, Scott L. Lipnick, Brendan T. Innes, Sean K. Simmons, Xian Adiconis, Danielle Dionne, Brittany A. Mayweather, et al. "Single-Cell Transcriptomic Profiling of the Aging Mouse Brain" https://doi.org/10.1038/s41593-019-0491-3  Nature Neuroscience 22, no. 10 (October 2019)

- snRNA-seq of Alzheimer's disease. 24 disease and 24 healthy individuals. 80,660 cells, 75,060 after filtering. ROSMAP study. [News press release](https://news.mit.edu/2019/gene-map-alzheimers-brain-0501) [snRNA-seq data](https://www.radc.rush.edu/docs/omics.htm), analysis results in [Supplementary Material](https://www.nature.com/articles/s41586-019-1195-2#Sec33)
    - [Supplementary Table 2](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1195-2/MediaObjects/41586_2019_1195_MOESM4_ESM.xlsx) - cell type-specific differentially expressed genes
    - [Supplementary Table 3](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1195-2/MediaObjects/41586_2019_1195_MOESM4_ESM.xlsx) - clinical annotations
    - [Supplementary Table 6](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1195-2/MediaObjects/41586_2019_1195_MOESM8_ESM.xlsx) - Cell-type sub-population signature genes
    Mathys, Hansruedi, Jose Davila-Velderrain, Zhuyu Peng, Fan Gao, Shahin Mohammadi, Jennie Z. Young, Madhvi Menon, et al. "Single-Cell Transcriptomic Analysis of Alzheimer’s Disease" https://doi.org/10.1038/s41586-019-1195-2  Nature, May 1, 2019

- [scRNA-seq of the middle temporal gyrus of human cerebral cortex](http://celltypes.brain-map.org/rnaseq). Similarity with mouse cortex scRNA-seq, and differences (proportions, laminar distribution, gene expression, morphology). NeuN-labeled sorting to select for neuronal nuclei. SMART-Seq v4 library preparation. Downstream analysis includes PCA, clustering using Jaccard, Louvain community detection. 75 cell types. Gene expression matrix, clusters, more. [Interactive exploration (download)](https://viewer.cytosplore.org/). <details>
    <summary>Paper</summary>
    Hodge, Rebecca D, Trygve E Bakken, Jeremy A Miller, Kimberly A Smith, Eliza R Barkan, Lucas T Graybuck, Jennie L Close, et al. "Conserved Cell Types with Divergent Features between Human and Mouse Cortex" https://doi.org/10.1101/384826  Preprint. Neuroscience, August 5, 2018. 
</details>

- [Brain immune atlas scRNA-seq resource](http://www.brainimmuneatlas.org/). Border-associated macrophages from discrete mouse brain compartments, tissue-specific transcriptional signatures. [GitHub](https://github.com/saeyslab/brainimmuneatlas/). <details>
    <summary>Paper</summary>
    Van Hove, Hannah, Liesbet Martens, Isabelle Scheyltjens, Karen De Vlaminck, Ana Rita Pombo Antunes, Sofie De Prijck, Niels Vandamme, et al. "A Single-Cell Atlas of Mouse Brain Macrophages Reveals Unique Transcriptional Identities Shaped by Ontogeny and Tissue Environment" https://doi.org/10.1038/s41593-019-0393-4  Nature Neuroscience, May 6, 2019. 
</details>

- Cell-type specificity of schizophrenia SNPs judged by enrichment in expressed genes. scRNA-seq custom data collection. Difference between schizophrenia and neurological disorders.
    - [data/Brain_cell_type_gene_expression.xlsx](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0129-5/MediaObjects/41588_2018_129_MOESM3_ESM.xlsx) - Supplementary Table 4 - Specificity values for Karolinska scRNA-seq superset. Specificity represents the proportion of the total expression of a gene found in one cell type as compared to that in all cell types (i.e., the mean expression in one cell type divided by the mean expression in all cell types). Gene X cell type matrix. Level 1 (core cell types) and level 2 (extended collection of cell types) data. 
    Major Depressive Disorder Working Group of the Psychiatric Genomics Consortium et al., "Genetic Identification of Brain Cell Types Underlying Schizophrenia" https://doi.org/10.1038/s41588-018-0129-5  Nature Genetics 50, no. 6 (June 2018) 

- [Drop-seq scRNA-seq data of 690,000 cells from 9 regions of adult mouse brain](http://dropviz.org/). Independent Component Analysis (ICA). ICs grouped into 565 transcriptionally distinct clusters (323 neuronal) corresponding to biological signals using network-based clustering. [DropViz](http://dropviz.org/) visualization, data download in CSV, RData, with annotations. <details>
    <summary>Paper</summary>
    Saunders, Arpiar, Evan Z. Macosko, Alec Wysoker, Melissa Goldman, Fenna M. Krienen, Heather de Rivera, Elizabeth Bien, et al. "Molecular Diversity and Specializations among the Cells of the Adult Mouse Brain" https://doi.org/10.1016/j.cell.2018.07.028  Cell 174, no. 4 (09 2018)
</details>

- Single-cell methylation of human and mouse neuronal cells. Marker genes with cell type-specific methylation profiles - [Table S3](http://science.sciencemag.org/content/suppl/2017/08/09/357.6351.600.DC1). <details>
    <summary>Paper</summary>
    Luo, Chongyuan, Christopher L. Keown, Laurie Kurihara, Jingtian Zhou, Yupeng He, Junhao Li, Rosa Castanon, et al. "Single-Cell Methylomes Identify Neuronal Subtypes and Regulatory Elements in Mammalian Cortex" https://doi.org/10.1126/science.aan3351  Science (11 2017) 
</details>

- Single-cell RNA-seq of human neuronal cell types. Dimensionality reduction, clustering, WGCNA, defining cell type-specific signatures, comparison with other signatures (Zeng, Miller). [Supplementary material](http://science.sciencemag.org/content/suppl/2017/12/06/358.6368.1318.DC1). [Wired story](https://www.wired.com/story/neuroscientists-just-launched-an-atlas-of-the-developing-human-brain/). [Controlled access data on dbGAP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000989.v3.p1), [summarized matrix with annotations](https://cells.ucsc.edu/?ds=cortex-dev)
    - [data/Brain/Nowakowski_2017_Tables_S1-S11.xlsx] - Table S5 has brain region-specific gene signatures. [Source](http://science.sciencemag.org/highwire/filestream/703290/field_highwire_adjunct_files/1/aap8809_Nowakowski_SM-Tables-S1-S11.xlsx) 
    Nowakowski, Tomasz J., Aparna Bhaduri, Alex A. Pollen, Beatriz Alvarado, Mohammed A. Mostajo-Radji, Elizabeth Di Lullo, Maximilian Haeussler, et al. "Spatiotemporal Gene Expression Trajectories Reveal Developmental Hierarchies of the Human Cortex" https://doi.org/10.1126/science.aap8809  Science, (December 8, 2017) 

- Single-nucleus droplet-based sequencing (snDrop-seq) and single-cell transposome hypersensitive site sequencing (scTHS-seq) of >60K cells from various parts of human adult brain. Resolving subpopulations, integrating the datasets, predicting one modality from another using the GBM classifier, integration with GWAS signal. Detailed methods, Data processing using Pagoda2, Seurat, LIGER, Monocle, [Code](https://github.com/JEFworks/Supplementary-Code), [Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97942), [Supplementary Table 3](https://www.nature.com/articles/nbt.4038#Sec35) - neuronal subpopulation-specific gene lists. <details>
    <summary>Paper</summary>
    Lake, Blue B, Song Chen, Brandon C Sos, Jean Fan, Gwendolyn E Kaeser, Yun C Yung, Thu E Duong, et al. "Integrative Single-Cell Analysis of Transcriptional and Epigenetic States in the Human Adult Brain" https://doi.org/10.1038/nbt.4038  Nature Biotechnology 36, no. 1 (December 11, 2017)
</details>

- scRNA-seq of human and mouse embryonic midbrain during development. Comparison of cell types and their temporal patterns. [BackSPIN clustering algorithm](https://github.com/linnarsson-lab/BackSPIN). [Raw and processed data (approx. 2000 cells per organism)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76381). [Code to reproduce the results](https://github.com/linnarsson-lab/ipynb-lamanno2016). <details>
    <summary>Paper</summary>
    La Manno, Gioele, Daniel Gyllborg, Simone Codeluppi, Kaneyasu Nishimura, Carmen Salto, Amit Zeisel, Lars E. Borm, et al. "Molecular Diversity of Midbrain Development in Mouse, Human, and Stem Cells" https://doi.org/10.1016/j.cell.2016.09.027  Cell 167, no. 2 (October 6, 2016)
</details>

- Single cell brain transcriptomics, human. Fluidigm C1 platform. Healthy cortex cells (466 cells) containing: Astrocytes, oligodendrocytes, oligodendrocyte precursor cells (OPCs), neurons, microglia, and vascular cells. Single cells clustered into 10 clusters, their top 20 gene signatures are in Supplementary Table S3. [Raw data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835).
    - [data/Brain/TableS3.txt] - top 20 cell type-specific genes
    - [data/Brain/TableS3_matrix.txt] - genes vs. cell types with 0/1 indicator variables.
    Darmanis, S., Sloan, S.A., Zhang, Y., Enge, M., Caneda, C., Shuer, L.M., Hayden Gephart, M.G., Barres, B.A., and Quake, S.R. (2015). [A survey of human brain transcriptome diversity at the single cell level](https://www.pnas.org/content/112/23/7285). Proc. Natl. Acad. Sci. 

## Links

- [Single Cell Genomics Day](https://satijalab.org/scgd21/) - virtual practical workshops by Satija lab, 2021 edition, links to videos and material from previous years. [Single Cell Genomics Day 2020 - Overview](https://youtu.be/b9bqKulMLp8) - video, 38 min 

- [Seurat object presentation](https://osf.io/49q2u/) by Ming Tang, [Tweet](https://twitter.com/tangming2005/status/1224726014170787842?s=20)

- [Loompy](http://loompy.org/) - file format for large omics datasets. [Linnarson lab](http://loom.linnarssonlab.org/) website includes scRNA-seq datasets in `.loom` format. Example R script to read the data in, [loom.R](tools/loom.R)

- [loupeR](https://github.com/10XGenomics/loupeR) - an R package to Convert Seurat objects to 10x Genomics Loupe files.

- [EDS](https://github.com/COMBINE-lab/EDS) - A simple, intuitive and Efficient single cell binary Data Storage format. Converter between different formats. 

- [Figure depicting the breadth of multimodal scRNA-seq technologies](https://github.com/arnavm/multimodal-scRNA-seq)

- [Clustered Dot Plots in the ggverse](https://davemcg.github.io/post/lets-plot-scrna-dotplots/), by David McGaughey, [Tweet](https://twitter.com/David_McGaughey/status/1242088829655400448?s=20). Similar plot by Ming Tang, [clustered dotplot for single-cell RNAseq](https://divingintogeneticsandgenomics.rbind.io/post/clustered-dotplot-for-single-cell-rnaseq/), [Tweet](https://twitter.com/tangming2005/status/1369778269835395077?s=20)

### Papers

- Stuart, Tim, Andrew Butler, Paul Hoffman, Christoph Hafemeister, Efthymia Papalexi, William M. Mauck, Yuhan Hao, Marlon Stoeckius, Peter Smibert, and Rahul Satija. "[Comprehensive Integration of Single-Cell Data](https://doi.org/10.1016/j.cell.2019.05.031)"  Cell, (June 2019) - [Seurat v.3](https://satijalab.org/seurat/) paper. Integration of multiple scRNA-seq and other single-cell omics (spatial transcriptomics, scATAC-seq, immunophenotyping), including batch correction. Anchors as reference to harmonize multiple datasets. Canonical Correlation Analysis (CCA) coupled with Mutual Nearest Neighborhoors (MNN) to identify shared subpopulations across datasets. CCA to reduce dimensionality, search for MNN in the low-dimensional representation. Shared Nearest Neighbor (SNN) graphs to assess similarity between two cells. Outperforms scmap. Extensive validation on multiple datasets (Human Cell Atlas, STARmap mouse visual cortex spatial transcriptomics. Tabula Muris, 10X Genomics datasets, others in STAR methods). Data normalization, variable feature selection within- and between datasets, anchor identification using CCA (methods), their scoring, batch correction, label transfer, imputation. Methods correspond to details of each Seurat function. Preprocessing of real single-cell data. [GitHub with code for the paper](https://github.com/satijalab/Integration2019)

- [The single cell studies database](https://docs.google.com/spreadsheets/d/1En7-UV0k0laDiIfjFkdn7dggyR7jIk3WH8QgXaMOZF0/edit#gid=935747368), over 1000 studies. [Main database](https://www.nxn.se/single-cell-studies), [Tweet by Valentine Svensson](https://twitter.com/vallens/status/1263477436513206273?s=20)

- [scATACdb](https://docs.google.com/spreadsheets/d/1IJhYMysjt-gyFeG1P7BL0NbKkOl1QPKqljz7s6SLZ6U/edit#gid=1161067356) - list of scATAC-seq studies, Google Sheet by Caleb Lareau

- [Journal club on single-cell multimodal data technology and analysis](https://github.com/waldronlab/data-science-seminar/wiki/Single-cell-multimodal-data) - Data science seminar led by Levi Waldron

- Review of manifold learning-based methods for denoising scRNA-seq data, revealing gene interactions, extracting pseudotime progressions with model fitting, visualizing the cellular state space via dimensionality reduction, and clustering. Manifold = a mathematical construct that represents a locally-Euclidean smoothly varying space. Modeling single-cell data as a manifold: 1. a graph based on local affinities such as the minimal spanning tree or a k- nearest neighbors (nn); 2. data diffusion. Applications: 1. Denoising and gene interactions (Seurat, ZIFA, CIDR, PCoA, scIMPUTE, SAVER, MAGIC, DREMI); 2. Pseudotime (Monocle, Wishbone, Wanderlust, Diffusion Pseudotime); 3. Dimensionality reduction (t-SNE, diffusion map, PHATE); 4. Density estimation and clustering (k-nn graph, PhenoGraph). Future directions in manifold learning. <details>
    <summary>Paper</summary>
    Moon, Kevin R., Jay S. Stanley, Daniel Burkhardt, David van Dijk, Guy Wolf, and Smita Krishnaswamy. “Manifold Learning-Based Methods for Analyzing Single-Cell RNA-Sequencing Data.” Current Opinion in Systems Biology 7 (February 2018): 36–46. https://doi.org/10.1016/j.coisb.2017.12.008.
</details>

- Review of single-cell sequencing technologies, individual and combined, technical details of each. Combinatorial indexing. Genomic DNA, methylomes, histone modifications, open chromatin, 3D genomics, proteomics, spatial transcriptomics. Table 1 - multiomics technologies, summary. Areas of application, in cancer and cell atlases. Future development, e.g., single-cell metabolomics. <details>
    <summary>Paper</summary>
    Chappell, Lia, Andrew J. C. Russell, and Thierry Voet. "Single-Cell (Multi)Omics Technologies" https://doi.org/10.1146/annurev-genom-091416-035324  Annual Review of Genomics and Human Genetics 19 (31 2018)
</details>

- [sciRNA-seq](http://atlas.gs.washington.edu/hub/) - single-cell combinatorial indexing RNA-seq technology and sequencing of C. elegans, ~49,000 cells, 27 cell types. Data and R code. <details>
    <summary>Paper</summary>
    Cao, Junyue, Jonathan S. Packer, Vijay Ramani, Darren A. Cusanovich, Chau Huynh, Riza Daza, Xiaojie Qiu, et al. "Comprehensive Single-Cell Transcriptional Profiling of a Multicellular Organism" https://doi.org/10.1126/science.aam8940  Science 357, no. 6352 (August 18, 2017)
</details>
