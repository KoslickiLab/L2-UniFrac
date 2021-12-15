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

Git and Curl Setup
```
sudo apt update && sudo apt upgrade
sudo apt-get install aptitude
sudo aptitude install git
sudo aptitude install curl
sudo aptitude install build essential
```

Krona Setup:
```
git clone https://github.com/marbl/Krona.git
cd Krona/KronaTools/
sudo ./install.pl
mkdir taxonomy/
sudo ./updateTaxonomy.sh #If fails to download taxdump.tar.gz, download directly from https://ftp.ncbi.nih.gov/pub/taxonomy/ and move to taxonomy/ folder.
sudo ./updateAccessions.sh #NOT REQUIRED. If fails, mkdir accession2taxid/ in taxonomy, and download dead_nucl.accession2taxid.gz, dead_prot.accession2taxid.gz, dead_wgs.accession2taxid.gz, nucl_gb.accession2taxid.gz, nucl_wgs.accession2taxid.gz, and prot.accession2taxid.gz from https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid
```

L2-UniFrac Install
```
git clone https://github.com/KoslickiLab/L2-UniFrac.git
```

Python Dependencies
```
sudo apt-get update
sudo apt-key adv --refresh-keys --keyserver keyserver.ubuntu.com
sudo apt-get -y install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt-get update
sudo apt-get -y install python3.9 python3.9-distutils
sudo apt remove python3.5-minimal -y
sudo ln -sf /usr/bin/python3.9 /usr/bin/python3
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3.9 get-pip.py
rm get-pip.py
conda install -c conda-forge scikit-bio
cd L2-UniFrac/
python -m pip install -r requirements.txt
```
