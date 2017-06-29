#!/bin/bash


# Install clustalw

sudo apt-get install clustalw

# Download Silva databases

mkdir SILVA_db
cd SILVA_db
wget https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_SSURef_tax_silva.fasta.gz
gunpzip SILVA_128_SSURef_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/taxonomy/tax_slv_ssu_128.txt

# Add bin folder to PATH

cd ../bin

path_to_NG_tax=$(pwd)
echo -e "\n#Added by NG-Tax installer" >>  ~/.bashrc
echo -e 'export PATH="$PATH:'$path_to_NG_tax'"\n' >>  ~/.bashrc
source  ~/.bashrc

cd ../test_set

# Update changes

bash
