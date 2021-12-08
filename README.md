# L2-UniFrac

L2-UniFrac is an improvement upon standard L1-UniFrac, otherwise known as EMDUniFrac, that pushes up the L2 normalization of the differential abundance vectors as opposed to the standard L1 normalization process (absolute values vs sum of squares). EMDUniFrac was first developed by Dr. Jason McClelland and Dr. David Koslicki and can be accessed via at either [arXiv](https://arxiv.org/abs/1611.04634) or [bioRxiv](https://www.biorxiv.org/content/10.1101/087171v2). This improvement was deemed necessary following issues regarding the presence of negative abundances in samples while testing data using the L1-UniFrac methodology. 

Later noted by Dr. McClelland, L2 metrics hold more biological significance due to the fact that their derivation directly yields a function that, for a given community, takes the sum of pairwise UniFrac distances times their relative abundances and subtracts the average spread of the same metric within a community and between communities. As a result, this metric allows researchers to quantify the evolutionary spread of a given community versus a group of communities, allowing for more differentiation among distinct species. His work related to Wasserstein β-Diversity metrics and more specifically for this repository, L2 normalization of such, can be found in the [Oregon State University Scholar's Archives](https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/76537620h) 

To retrieve more biologiocally meaningful results that preserve fundamental characteristics such as non-negative sample abundances in averaged sets, L2-UniFrac set out to solve these aforementioned shortcomings with L1-UniFrac and apply the principles noted by Dr. McClelland to real-world data.

## Requirements ##
+ [dendropy](http://www.dendropy.org/)
+ numpy 
+ matplotlib - for plotting
+ scipy 
+ [scikit-bio](http://scikit-bio.org/) - for PCoA and Distance Matrix
+ [pandas](https://pandas.pydata.org/) - for metadata frame
+ [biom-format](https://biom-format.org/) - for chosen dataset
+ h5py

## Setup ##

Miniconda Setup:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +777 Miniconda3-latest-Linux-x86_64.sh
bash
./Miniconda3-latest-Linux-x86_64.sh
conda create -n env1
conda activate env1
conda install -c conda-forge xorg-makedepend
```

Git Setup
```
sudo apt-get install aptitude
sudo aptitude install git
```

Krona Setup:
```
git clone https://github.com/marbl/Krona.git
```
