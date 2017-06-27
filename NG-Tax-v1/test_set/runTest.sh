###### NG-Tax test ######
### Pair-end fastq libraries: library_test_1.fastq, library_test_2.fastq
### Mapping file: mapping_file_test.txt
### 16S rRNA gene database: test_db.fa
### 16S rRNA gene taxonomy: test_taxonomy.txt

### Test code

# Remove biom file if exist
rm -f mapping_file_test.biom

# Create customized 16 rRNA gene databases.
customized_16S_database_generator.sh -d test_db.fa -t test_taxonomy.txt -k GTGCCAGC[AC]GCCGCGGTAA -p GGACTAC[ACT][ACG]GGGT[AT]TCTAAT -q ATTAGA[AT]ACCC[TCG][ATG]GTAGTCC -f primer_for_test_db -r primer_rev_test_db -x test_NG_Tax -o 71 -e 70

# Filter pair-end fastq libraries.
library_filtering.sh -a library_test_1.fastq -b library_test_2.fastq -p library_test -n 01 -f GTGCCAGC[AC]GCCGCGGTAA -r GGACTAC[ACT][ACG]GGGT[AT]TCTAAT -l 8

# OTU picking.
OTU_picking_pair_end.sh -m mapping_file_test.txt -p library_test -a 0.1 -c 0.985 -f primer_for_test_db -r primer_rev_test_db -o 71 -e 70 -t test_NG_Tax -q 2 -k 100 -n 4

# Make biom file.
make_biom_file.sh -t tax_files -m mapping_file_test.txt

# Check for biom file
if [[ ! -s mapping_file_test.biom ]]
	then
		echo -e "\nERROR: mapping_file_test.biom not present or empty\n"
		exit 1
else
  then
    echo -e "\nTest passed succesfully\n"
fi
