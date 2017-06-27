#!/bin/bash

# Usage message

usage="
\n---USAGE---\n
\n
    Description:\n\n
This script reassigns OTUs to alternative taxonomies introduced by the user. \n
The new taxonomy should be recorded in the 3rd column of the alternative \n
reassignment file. If the word remove is recorded instead of a new taxonomy \n
the OTU will be removed. /n
\n
    Required inputs\n\n
-a alternative taxonomy file \n
-t tax folder \n
\n
    Output:\n
\n
A folder containing reassigned tax files \n\n
Example of usage:\n
$0 -a alternative_taxonomy_file -t tax_files \n\n
  output:\n
alternative_reassigned_tax_files\n\n
--- End of USAGE ---\n
"

# Help option, check number of argument.

if [ "$1" == "-h" ]
  then
    echo -e $usage
    exit 0
elif [[ $# -ne 4 ]]
  then
    echo -e "ERROR: Invalid number of args \nTry $0 -h for help"
    exit 1
fi

# Assign arguments to variables, check for invalid arguments.

while getopts a:t: opt
  do
    case "$opt" in
    a)
      alternative_taxonomy=${OPTARG}
      ;;
    t)
      tax_folder=${OPTARG}
      ;;
    *)
      echo -e "\nERROR: invalid arguments used\n";
      exit 1
      ;;
    esac
done

# Check argument files presence.

if [[ ! -s $alternative_taxonomy  ]]
  then
  	echo -e "\nERROR: $alternative_taxonomy not present or empty\n"
  	exit 1
elif [[ ! -d $tax_folder ]]
  then
    echo -e "\nERROR: $tax_folder not present or not directory \n"
    exit 1
fi

# Create NG-tax folder.

mkdir -p "alternative_reassigned_"$tax_folder

# Create alternative reassigned files.

for file in $tax_folder/*
  do
    awk -v v_taxonomy=$alternative_taxonomy 'BEGIN{ \
      while(( getline line<v_taxonomy) > 0 ) { \
        split(line,name,"\t"); \
        sequence=name[1]; \
        array_remove[sequence]=name[3]; \
        array_taxonomy[sequence]=name[3]"\t"name[4]"\t"name[5]"\t"name[6] \
      } \
    } \
    { \
      if(substr($1,1,1)==">"){ \
        otu=$1; \
        tax=$2"\t"$3"\t"$4"\t"$5 \
      } \
      else{ \
        if(array_taxonomy[$1]!=""){ \
          if(array_remove[$1]!="remove"){ \
            print otu"\t"array_taxonomy[$1]; \
            print $1 \
          } \
        } \
        else{ \
          print otu"\t"tax; \
          print $1 \
        } \
      } \
    }' $file > "alternative_reassigned_"$file
done
