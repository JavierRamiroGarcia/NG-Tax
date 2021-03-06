#!/bin/bash

# Usage message

usage="
\n---USAGE---\n
\n
    Description:\n\n
This script downloads the Silva database and the taxonomy (release_128) named\n
SILVA_128_SSURef_tax_silva.fasta and tax_slv_ssu_128.txt. A different database\n
or taxonomy files can be provided by entering the link to them as optional input.\n\n
It also installs clustalw and add the bin folder to PATH\n
It should be run from the folder containing the script\n
\n
    Optional inputs\n
\n
-d Link to reference 16S rRNA database file [fasta]\n
-t Link to SILVA taxonomy\n
\n
    Output:\n
\n
A reference 16S rRNA database in fasta format \n
A taxonomy file customized for NG-Tax\n\n
Example of usage:\n\n
Default -> Silva release_128:\n
$0
Providing links to a different Silva reference 16S database or to a different Silva taxonomy:\n
$0 -d link_to/SILVA_128_SSURef_tax_silva.fasta -t link_to/tax_slv_ssu_128.txt\n\n
output:\n
SILVA_128_SSURef_tax_silva.fasta\n
tax_slv_ssu_128.txt\n\n
--- End of USAGE ---\n
"

# Check if NG-Tax_installer_mac.sh is run from a different folder

if [[ "$0" != "./NG-Tax_installer_unix.sh" ]]
  then
    echo -e "ERROR: Move to the folder containing NG-Tax_installer_unix.sh\nTry ./NG-Tax_installer_unix.sh -h for help"
    exit 1
fi


# Help option, check number of argument.

if [[ "$1" == "-h" ]]
  then
    echo -e $usage
    exit 0
elif [[ $# -ne 0 && $# -ne 2 && $# -ne 4 ]]
  then
    echo -e "ERROR: Invalid number of args \nTry $0 -h for help"
    exit 1
fi

# Assign arguments to variables, check for invalid arguments.

while getopts d:t: opt
  do
  case "$opt" in
    d)
      database_link=${OPTARG}
      ;;
    t)
      taxonomy_link=${OPTARG}
      ;;
    *)
      echo -e "\nERROR: invalid arguments used\n";
      exit 1
      ;;
    esac
done


# Install wget gunzip or clustalw packages if they are not already installed

packages_to_install=( wget gzip gunzip clustalw )
for i in "${packages_to_install[@]}"
	do
		package_to_test=$i
		if ! type "$package_to_test" &> /dev/null
			then
				sudo apt-get install $package_to_test
		fi
done

# Download Silva databases using realease 128 as default. Other releases could be downloaded by providing the link as an optional input

mkdir SILVA_db
cd SILVA_db

if [[ "$database_link" != "" ]]

  then
    wget SILVA_db $database_link
    bname_db=$(basename "$database_link")
    gunzip $bname_db
else
    wget SILVA_db https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_SSURef_tax_silva.fasta.gz
    gunzip SILVA_128_SSURef_tax_silva.fasta.gz
fi

if [[ "$taxonomy_link" != "" ]]

  then
    wget SILVA_db $taxonomy_link
else
    wget SILVA_db https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/taxonomy/tax_slv_ssu_128.txt
fi

cd ..

# Add bin folder to PATH

path_to_NG_tax=$(pwd)"/bin"
echo -e "\n#Added by NG-Tax installer" >>  ~/.bashrc
echo -e 'export PATH="$PATH:'$path_to_NG_tax'"\n' >>  ~/.bashrc
source  ~/.bashrc

cd test_set

# Update changes

bash
