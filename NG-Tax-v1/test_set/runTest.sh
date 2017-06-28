#!/bin/bash

# Check for usearch file
if [[ ! -s ../bin/usearch ]]
	then
		echo -e "\nERROR: usearch script not present or empty\n"
		exit 1
fi

# Remove biom file if exist
rm -f mapping_file_test.biom

# Create customized 16 rRNA gene databases.
customized_database_generator.sh -d test_db.fa -t test_taxonomy.txt -k GTGCCAGC[AC]GCCGCGGTAA -p GGACTAC[ACT][ACG]GGGT[AT]TCTAAT -q ATTAGA[AT]ACCC[TCG][ATG]GTAGTCC -f primer_for_test_db -r primer_rev_test_db -x test_taxonomy_NG_Tax -o 71 -e 70

# Filter pair-end fastq libraries.
library_filtering.sh -a library_test_1.fastq -b library_test_2.fastq -p library_test -n 01 -f GTGCCAGC[AC]GCCGCGGTAA -r GGACTAC[ACT][ACG]GGGT[AT]TCTAAT -l 8

# OTU picking.
OTU_picking_pair_end.sh -m mapping_file_test.txt -p library_test -a 0.1 -c 0.985 -f primer_for_test_db -r primer_rev_test_db -o 71 -e 70 -t test_taxonomy_NG_Tax -q 2 -k 100 -n 4

# Check for tax files
if [[ ! -s tax_files/Sample.talarrubias_tax_file ]]
	then
		echo -e "\nERROR: tax_files not present or empty\n"
		exit 1
fi

# Examples of tax files:
echo -e "Example of tax files."
echo -e "\nSample.talarrubias:"
cat tax_files/Sample.talarrubias_tax_file
echo -e "\nSample.trigueros:"
cat tax_files/Sample.trigueros_tax_file

# Make biom file.
make_biom_file.sh -t tax_files -m mapping_file_test.txt

# Check for biom file
if [[ ! -s mapping_file_test.biom ]]
	then
		echo -e "\nERROR: mapping_file_test.biom not present or empty\n"
		exit 1
else
    echo -e "\nTest passed succesfully\n"
fi

cd ..
