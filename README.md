# L2-UniFrac

L2-UniFrac is an improvement upon standard L1-UniFrac, otherwise known as EMDUniFrac, that pushes up the L2 normalization of the differential abundance vectors as opposed to the standard L1 normalization process (absolute values vs sum of squares). EMDUniFrac was first developed by Dr. Jason McClelland and Dr. David Koslicki and can be accessed via at either [arXiv](https://arxiv.org/abs/1611.04634) or [bioRxiv](https://www.biorxiv.org/content/10.1101/087171v2). This improvement was deemed necessary following issues regarding the presence of negative abundances in samples while testing data using the L1-UniFrac methodology. 

Later noted by Dr. McClelland, L2 metrics hold more biological significance due to the fact that their derivation directly yields a function that, for a given community, takes the sum of pairwise UniFrac distances times their relative abundances and subtracts the average spread of the same metric within a community and between communities. As a result, this metric allows researchers to quantify the evolutionary spread of a given community versus a group of communities, allowing for more differentiation among distinct species. His work related to Wasserstein β-Diversity metrics and more specifically for this repository, L2 normalization of such, can be found in the [Oregon State University Scholar's Archives](https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/76537620h) 

To further demonstrate the usefulness of L2-UniFrac in finding the biologically meaningful average sample of an environment, another set of results were published by Wei Wei, Andrew Millward and Dr. Koslicki in 2023, which can be found [here](https://academic.oup.com/bioinformatics/article/39/Supplement_1/i57/7210517). This repo now contains the supporting functions for this publication, which also serves as a prototype for users to explore the potential usages of L2-UniFrac.

## Installation ##

```
git clone -b user --single-branch https://github.com/KoslickiLab/L2-UniFrac.git
cd L2-UniFrac
conda create -n L2UniFrac python=3.10
conda activate L2UniFrac
python -m pip install -r requirements.txt
```

## Usage ##
### 1. Finding the average sample with respect to L2UniFrac - 16s data
- Required input:
  - Metagenomic samples in OTU table format. An example of such a table can be found in the data directory.
  - A metadata file that specifies the phenotype/environment each sample. The user is also required to specify the column of this file based on which the average samples will be generated.
  - Tree file. The default option is `gg_13_5_otus_99_annotated.tree` from Greengenes. 

- Output:
  - An OTU file containing the average sample for each of the environments/phenotypes under the specified column of the metadata file.

#### Example
```
python scripts/get_16s_L2UniFrac_average.py -i data/example_data/otu_table_body_sites.biom -m data/example_data/metadata_body_sites.txt -k sample_name -v body_site -o data/example_output/representative_otu_by_body_site.tsv 
```

The output after running the above command will be saved in `data/example_output`.

### 2. Finding the average sample with respect to L2UniFrac - WGS data
This extension is based on the method [WGSUniFrac](https://drops.dagstuhl.de/opus/volltexte/2022/17049/). The input and output will be in the format of CAMI profiles.

- Required input:
  - A directory containing WGS profiles, one for each sample. An example of such files can be found under `data/example_data/adenoma_266076/profiles`.
  - A metadata file that specifies the phenotype/environment each sample. The user is also required to specify the column of this file based on which the average samples will be generated.
  - Output format and path

- Output:
  - The user has two choices of output format, the CAMI profile format in which a each representative sample is represented as one single profile in the directory specified by the user. Or, an OTU table format in which the representative samples will be saved in one single OTU table under the name specified by the user.

### 3. Cluster in L2UniFrac space
- To be continued...