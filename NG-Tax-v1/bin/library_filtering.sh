#!/bin/bash


# Usage message

usage="
\n---USAGE---\n
\n
    Description:\n\n
This script creates a library suitable for NG-Tax by filtering those reads with not  \n
matching barcodes from a pair-end library. For every sample, forward and reverse reads \n
should have the same barcode. Last nucleotide of every read is also removed for quality \n
reasons.\n
\n
    Required inputs\n\n
-a sequencing library file 1 [fastq]\n
-b sequencing library file 2 [fastq]\n
-p name given to the project\n
-n number given to the library (as indicated in the third column of the mapping file)\n
-f forward primer sequence (degenerate positions between brackets)\n
-r reverse primer sequence (degenerate positions between brackets)\n
-l barcode length\n
\n
    Output:\n
\n
A filtered library \n\n
Example of usage:\n
$0 -a lib1_1.fastq -b lib1_2.fastq -p Mock_project -n 01 -f GTGCCAGC[AC]GCCGCGGTAA -r GGACTAC[ACT][ACG]GGGT[AT]TCTAAT -l 8 \n\n
  output:\n
Mock_project_01\n\n
--- End of USAGE ---\n
"

# Help option, check number of argument.

if [ "$1" == "-h" ]
  then
    echo -e $usage
    exit 0
elif [[ $# -ne 14 ]]
  then
    echo -e "ERROR: Invalid number of args \nTry $0 -h for help"
    exit 1
fi

# Assign arguments to variables, check for invalid arguments.

while getopts a:b:p:n:f:r:l: opt
  do
  case "$opt" in
    a)
      forward_library=${OPTARG}
      ;;
    b)
      reverse_library=${OPTARG}
      ;;
    p)
      project_name=${OPTARG}
      ;;
    n)
      library_number=${OPTARG}
      ;;
    f)
      forward_primer=${OPTARG}
      ;;
    r)
      reverse_primer=${OPTARG}
      ;;
    l)
      barcode_length=${OPTARG}
      ;;
    *)
    echo -e "\nERROR: invalid arguments used\n";
    exit 1
    ;;
    esac
done

# Check library files presence.

if [[ ! -s $forward_library ]]
	then
		echo -e "\nERROR: $forward_library not present or empty\n"
		exit 1

elif [[ ! -s $reverse_library  ]]
	then
		echo -e "\nERROR: $reverse_library not present or empty\n"
		exit 1
fi

# Create filtered library

paste  $forward_library $reverse_library | \
awk -v v_forward_primer=$forward_primer -v v_reverse_primer=$reverse_primer -v v_barcode_length=$barcode_length -v v_library_number=$library_number \
	'BEGIN{ \
    i=0 \
  } \
  { \
    i++; \
    if(i==1){ \
      name=$1;allow_print=0 \
    }; \
    if(i==2){ \
      if(substr($1,1,v_barcode_length)==substr($2,1,v_barcode_length)){ \
        if(match($1,v_forward_primer)>0 && match($2,v_reverse_primer)>0){ \
          forward_seq=$1; \
					reverse_seq=$2; \
					allow_print=1 \
        }; \
        if(match($2,v_forward_primer)>0 && match($1,v_reverse_primer)>0){ \
          forward_seq=$2; \
          reverse_seq=$1; \
					allow_print=1 \
        } \
      } \
    }; \
    if(i==4){ \
      i=0; \
      if(allow_print==1){ \
        print name"\t"v_library_number""substr(forward_seq,1,v_barcode_length)"\t"substr(forward_seq,1,length(forward_seq)-1)"\t"substr(reverse_seq,1,length(reverse_seq)-1) \
      } \
    } \
  }' | \
sed "s/$forward_primer/\\t/" | \
sed "s/$reverse_primer/\\t/" | \
awk '{print $1"\t"$2"\t"$4"\t"$6; \
}' > $project_name"_"$library_number
