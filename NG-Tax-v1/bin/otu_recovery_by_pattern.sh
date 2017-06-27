#!/bin/bash

# Usage message

usage="
\n---USAGE---\n
\n
    Description:\n\n
This script retrieves OTUs by searching a pattern in the sequence or in the taxonomic assignment.\n
\n
    Required inputs\n\n
-t folder containing the tax files \n
-n pattern name \n
-p pattern to be retrieved \n
-s search by taxonomy or sequence \n
\n
    Output:\n
\n
A folder containing: \n
\n
A file with all the OTUs matching the pattern \n
A file with all the unique OTUs matching the pattern [fasta] \n
A file to be used for alternative reassignment \n\n
Example of usage:\n
$0 -t tax_files -n non_assigned -p NA -s taxonomy \n\n
  output:\n
otu_retrievement_files/non_assigned_otus_file \n
otu_retrievement_files/non_assigned_unique_otus_file.fasta \n
otu_retrievement_files/non_assigned_alternative_taxonomy_file \n\n
--- End of USAGE ---\n
"

# Help option, check number of argument.

if [ "$1" == "-h" ]
  then
    echo -e $usage
    exit 0
elif [[ $# -ne 8 ]]
  then
    echo -e "ERROR: Invalid number of args \nTry $0 -h for help"
    exit 1
fi

# Assign arguments to variables, check for invalid arguments.

while getopts t:n:p:s: opt
 do
  case "$opt" in
  t)
    tax_files=${OPTARG}
    ;;
  n)
    pattern_name=${OPTARG}
    ;;
  p)
    pattern_to_search=${OPTARG}
    ;;
  s)
    taxonomy_or_sequence=${OPTARG}
    ;;
  *)
    echo -e "\nERROR: invalid arguments used\n";
    exit 1
    ;;
  esac
done

# Check argument files presence.

if [[ ! -d $tax_files  ]]
  then
    echo -e "\nERROR: $tax_files not present or not directory \n"
    exit 1
elif [[ $taxonomy_or_sequence != "taxonomy" && $taxonomy_or_sequence != "sequence" ]]
  then
    echo -e "\nERROR: $taxonomy_or_sequence value not valid. Valid values are taxonomy or sequence \n"
    exit 1
fi

# Create NG-tax folder

mkdir -p "otu_retrievement_files"

# Search pattern by taxonomy

if [[ "$taxonomy_or_sequence" == "taxonomy" ]]
  then
    cat $tax_files"/"* | \
    grep $pattern_to_search -A1 --no-group-separator | \
    grep -v \^\> --no-group-separator | \
    sort | \
    uniq | \
    awk '{\
      i+=1;\
      print ">seq_"i;\
      print $0 \
    }' > "otu_retrievement_files/"$pattern_name"_unique_otus_file.fasta"

    cat $tax_files"/"* | \
    grep $pattern_to_search -A1 --no-group-separator | \
    grep -v \^\> --no-group-separator | \
    sort | \
    uniq | \
    awk '{ \
      i+=1; \
      print $0"\t>seq_"i"\tNA;__p;__c;__o;__f;__g\tNA\tNA\tNA\tExplanation" \
    }' > "otu_retrievement_files/"$pattern_name"_alternative_taxonomy_file"

    cat $tax_files"/"* | \
    grep $pattern_to_search -A1 --no-group-separator | \
    awk '{ \
      if(substr($1,1,1)==">"){\
        name=$0 \
      } \
      else{ \
        print $1"\t"name \
      } \
    }' | \
    sort | \
    awk 'BEGIN{\
      i=0; \
      sequence="" \
    } \
    {\
      if($1!=sequence){ \
        i+=1; \
        sequence=$1 \
      }; \
      print $1"\t>seq_"i"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 \
    }' > "otu_retrievement_files/"$pattern_name"_otus_file"
fi

# Search pattern by sequence.

if [[ "$taxonomy_or_sequence" == "sequence" ]]
  then
    cat $tax_files"/"* | \
    grep $pattern_to_search -B1 --no-group-separator | \
    grep -v \^\> --no-group-separator | \
    sort | \
    uniq | \
    awk '{ \
      i+=1; \
      print ">seq_"i; \
      print $0 \
    }' > "otu_retrievement_files/"$pattern_name"_unique_otus_file.fasta"

    cat $tax_files"/"* | \
    grep $pattern_to_search -B1 --no-group-separator | \
    grep -v \^\> --no-group-separator | \
    sort | \
    uniq | \
    awk '{ \
      i+=1; \
      print $0"\t>seq_"i"\tNA;__p;__c;__o;__f;__g\tNA\tNA\tNA\tExplanation" \
    }' > "otu_retrievement_files/"$pattern_name"_alternative_taxonomy_file"

    cat $tax_files"/"* | \
    grep $pattern_to_search -B1 --no-group-separator | \
    awk '{ \
      if(substr($1,1,1)==">"){ \
        name=$0 \
      } \
      else{ \
        print $1"\t"name \
      } \
    }' | \
    sort | \
    awk '{ \
      if($1!=sequence){\
        i+=1; \
        sequence=$1 \
      }; \
      print $1"\t>seq_"i"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 \
    }' > "otu_retrievement_files/"$pattern_name"_otus_file"
fi
