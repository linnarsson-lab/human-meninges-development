# human-meninges-development
Supplementary code and data to our study of the developing human meninges by Vinsland et al. (manuscript in preparation).

## Preprint (bioRxiv)

TBA

## Browser

Visualizations of the complete dataset, as well as cell type specific subsets ("Superclasses"), are browsable at [CELLxGENE](https://cellxgene.cziscience.com/collections/7d66d871-091f-4602-9f42-85f86d2853e0).

## Data availability

#### Raw sequence reads

- BAM files will be available from the European Genome/Phenome Archive (https://ega-archive.org/) under accession number TBA. 

#### scRNA-seq expression matrices

- Complete count matrices (gene x cell counts) for the developing meninges and meningiomas are available as loom files [here](loom_files.md).
- The datasets can also be downloaded as .h5ad files from the browser: [CELLxGENE](https://cellxgene.cziscience.com/collections/7d66d871-091f-4602-9f42-85f86d2853e0). 

#### Xenium spatial data

- Raw Xenium data and images have been deposited at the [BioImage Archive](https://www.ebi.ac.uk/bioimage-archive/) under accession number S-BIAD1600.
- Complete count matrices (gene x cell counts) for meningiomas are available as loom files [here](loom_files.md).

## Code used for analysis and visualisation

- Clustering was performed using the cytograph-dev version of cytograph. This is the version used for our adult human brain project. Its installation and usage are described [here](https://github.com/linnarsson-lab/adult-human-brain/tree/main/cytograph). 
- [Jupyter](https://jupyter.org/) notebooks used to make figures are available [here](notebooks). The notebooks also import from cytograph-dev. (cytograph-shoji will *not* work).
- Jupyter notebooks used for Xenium data processing are found [here](https://storage.googleapis.com/linnarsson-lab-loom/meninges/Xenium_processing.zip).
