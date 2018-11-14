#!/bin/bash

# Usage message

usage="
\n---USAGE---\n
\n
    Description:\n\n
This script allows for comparison of two 16S rRNA regions by generating tax files with OTUs that have both, \n
the amplified region and the region you want to compare with. The folders both_regions_tax_files generated \n
for each of the regions should be merged into a new folder and then this folder can be used for generating \n
the biom file.\n
\n
    Required inputs\n
\n
-m mapping file containing metadata [txt]\n
-t folder containing tax files \n
-f NG-Tax customized forward-primer database file for the region you want to predict\n
-r NG-Tax customized reverse-primer database file for the region you want to predict\n
-p position (1 to attach the predicted region after the OTU sequence and 2 to attach the predicted region before) \n
\n
    Output:\n
\n
A folder containing region prediction files \n
A folder containing both regions tax files \n\n
Example of usage:\n
$0 -m map_Mocks.txt -t tax_files -f primer_F515_71_nt_1mm_db -r primer_R806_70_nt_1mm_db -p 1 \n\n
output:\n
region_prediction_files \n
both_regions_tax_files \n\n
--- End of USAGE ---\n
"

# Help option, check number of argument.

if [ "$1" == "-h" ]
  then
    echo -e $usage
    exit 0
elif [[ $# -ne 10 ]]
  then
    echo -e "ERROR: Invalid number of args \nTry $0 -h for help"
    exit 1
fi

# Assign arguments to variables, check for invalid arguments.

while getopts m:t:f:r:p: opt
  do
    case "$opt" in
    m)
      map=${OPTARG}
      ;;
    t)
      tax_folder=${OPTARG}
      ;;
    f)
      forward_primer_db=${OPTARG}
      ;;
    r)
      reverse_primer_db=${OPTARG}
      ;;
    p)
      position=${OPTARG}
      ;;
    *)
      echo -e "\nERROR: invalid arguments used\n";
      exit 1
      ;;
    esac
done

# Check argument files presence.

if [[ ! -s $map ]]
	then
  	echo -e "\nERROR: $map not present or empty\n"
		exit 1
elif [[ ! -s $forward_primer_db  ]]
  then
  	echo -e "\nERROR: $forward_primer_db not present or empty\n"
  	exit 1
elif [[ ! -s $reverse_primer_db  ]]
  then
    echo -e "\nERROR: $reverse_primer_db not present or empty\n"
  	exit 1
elif [[ ! -d $tax_folder  ]]
  then
    echo -e "\nERROR: $tax_folder not present or not a directory\n"
    exit 1
fi

# Create NG-tax folders

mkdir -p  region_prediction_files both_regions_tax_files

# Create region prediction files

awk  '{ \
  if(FNR==1){ \
    file_number++ \
  }; \
  if(file_number==1){ \
    if(substr($1,1,1)==">"){ \
      id_db=substr($1,2,length($1)-1) \
    } \
    else{ \
      array_seq_db_forward[id_db]=substr($1,1,length($1)-10) \
    } \
  }; \
  if(file_number==2){ \
    if(substr($1,1,1)==">"){ \
      id_db=substr($1,2,length($1)-1) \
    } \
    else{ \
      array_seq_db_reverse[id_db]=substr($1,1,length($1)-10) \
    } \
  }; \
  if(file_number==3){ \
    otu=substr($1,2,length($1)-1); \
    array_otu_identity_threshold[otu]=$3 \
  }; \
	if(file_number==4){ \
    if($1=="H" && array_otu_identity_threshold[$9]<=$4){ \
      array_otu_fr_hit[$9$10]=$9$10 \
    } \
  }; \
	if(file_number==5){ \
    if($1=="H" && array_otu_fr_hit[$9$10]==$9$10 && array_otu_identity_threshold[$9]<=$4){ \
      print $9"\t"array_seq_db_forward[$10] > "region_prediction_files/database_otus_region_prediction_fr_file_tmp" ; \
      print $9"\t"array_seq_db_reverse[$10] > "region_prediction_files/database_otus_region_prediction_rr_file_tmp" \
    } \
  } \
}' $forward_primer_db $reverse_primer_db "complementary_tax_files/database_tax_file" "clustering_results_files/database_results_otus_file_fr.uc" "clustering_results_files/database_results_otus_file_rr.uc"

cat <(sort "region_prediction_files/database_otus_region_prediction_fr_file_tmp" | uniq -c | sort -k2,2 -nk1,1) <( echo "" ) > "region_prediction_files/database_otus_region_prediction_fr_file"
cat <(sort "region_prediction_files/database_otus_region_prediction_rr_file_tmp" | uniq -c | sort -k2,2 -nk1,1) <( echo "" ) > "region_prediction_files/database_otus_region_prediction_rr_file"

rm "region_prediction_files/database_otus_region_prediction_fr_file_tmp" "region_prediction_files/database_otus_region_prediction_rr_file_tmp"

# Create both regions tax files

for sample_name in $(awk '{if(NR>1 && $1!=""){print $1}}' $map)
  do
    awk -v v_sample_name=$sample_name -v v_position=$position '{ \
      if(FNR==1){ \
        file_number++ \
      }; \
			if(file_number==1){ \
				if(otu!=$2 && otu!=""){ \
          print otu"\t"seq"\t"number"\t"number/count > "region_prediction_files/database_otus_seq_predicted_fr_file" ; \
          array_predicted_seq_fr[otu]=seq; \
          otu=$2; \
          seq=$3; \
          number=$1; \
          count=$1 \
        } \
        else{ \
          otu=$2; \
          seq=$3; \
          number=$1; \
          count+=$1 \
        } \
      }; \
			if(file_number==2){ \
				if(otu!=$2 && otu!=""){ \
          print otu"\t"seq"\t"number"\t"number/count > "region_prediction_files/database_otus_seq_predicted_rr_file" ; \
          array_predicted_seq_rr[otu]=seq; \
          otu=$2; \
          seq=$3; \
          number=$1; \
          count=$1 \
        } \
        else{ \
          otu=$2; \
          seq=$3; \
          number=$1; \
          count+=$1 \
        } \
      }; \
			if(file_number==3){ \
        if(substr($1,1,1)==">"){ \
          otu=substr($1,2,length($1)-1) \
        } \
        else{ \
          seq_database=$0; \
          array_otu_seq[seq_database]=otu \
        } \
      }; \
			if(file_number==4){ \
				if(substr($1,1,1)==">"){ \
          print $0 \
        } \
        else{ \
          otu_database=array_otu_seq[$1]; \
          if(v_position==1){ \
            print $1array_predicted_seq_fr[otu_database]array_predicted_seq_rr[otu_database] \
          }; \
          if(v_position==2){ \
            print array_predicted_seq_fr[otu_database]array_predicted_seq_rr[otu_database]$1 \
          } \
        } \
      } \
    }'  "region_prediction_files/database_otus_region_prediction_fr_file" "region_prediction_files/database_otus_region_prediction_rr_file" "all_otus_files/database_otus" $tax_folder"/"$sample_name"_tax_file" > "both_regions_tax_files/"$sample_name"_tax_file"

done
