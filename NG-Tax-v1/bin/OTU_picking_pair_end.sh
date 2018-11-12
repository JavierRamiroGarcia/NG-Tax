#!/bin/bash

# Usage message

usage="
\n---USAGE---\n
\n
    Description:\n\n
This script demultiplexes the raw data into samples using the information contained in the mapping file, \n
it also generates an OTU table per each sample after removing chimeras and assigns taxonomy to the OTUs. \n
NG-Tax is designed for short reads, 70-100 nucleotides is the recommended read length. Reads can be trimmed \n
to this length by the script. Longer length can be selected by the user but comparison with 70-100 nucleotides \n
results is advisable. \n
\n
    Required inputs\n
\n
-m mapping file containing metadata [txt]\n
-p NG-tax project name.\n
-a minimum threshold detectable, expressed in percentage (0.1 recommended) \n
-c error correction clustering percentage (only one mismatch recommended, 0.985 for ~70 nt)\n
-f NG-Tax customized forward-primer database file\n
-r NG-Tax customized reverse-primer database file\n
-o forward-read length\n
-e reverse-read length\n
-t taxonomic table file\n
-q ratio otu_parent_abundance/otu_chimera_abundance (recommended 2, both otu parents must be two times more abundant than the otu chimera) \n
-k identity level between parents and chimera (recommended 100, no error allowed, chimera as perfect combination of two otus) \n
-n number of threads (for parallel computing) \n
\n
    Output:\n
\n
A folder containing files with demultiplexed reads \n
A folder containing files with OTUs before chimera checking \n
A folder containing files with chimera clustering \n
A folder containing files with final OTUs \n
A folder containing files with final OTUs for all the samples \n
A folder containing the clustering results for all the OTUs against the 16S database \n
A folder containing the complementary taxonomic assignment for all the otus \n
A folder containing the final taxonomic assignment for every sample \n\n
Example of usage:\n
$0 -m map_Mocks.txt -p Mock_project -a 0.1 -c 0.985 -f primer_F515_71_nt_1mm_db -r primer_R806_70_nt_1mm_db -o 71 -e 70 -t Silva_111_taxa_map_RDP_6_levels_full.txt -q 2 -k 100 -n 24 \n\n
output:\n
total_sample_files\n
pre_set_otu_files\n
chimera_clustering_files\n
otus_files\n
all_otus_file\n
clustering_results_files\n
complementary_tax_files\n
tax_files \n\n
--- End of USAGE ---\n
"

# Help option, check number of argument.

if [ "$1" == "-h" ]
  then
    echo -e $usage
    exit 0
elif [[ $# -ne 24 ]]
  then
    echo -e "ERROR: Invalid number of args \nTry $0 -h for help"
    exit 1
fi

# Assign arguments to variables, check for invalid arguments.

while getopts m:p:a:c:f:r:o:e:t:q:k:n: opt
  do
    case "$opt" in
    m)
      map=${OPTARG}
      ;;
    p)
      project_name=${OPTARG}
      ;;
    a)
      abundance_threshold=${OPTARG}
      ;;
    c)
      clustering_percentage=${OPTARG}
      ;;
    f)
      forward_primer_db=${OPTARG}
      ;;
    r)
      reverse_primer_db=${OPTARG}
      ;;
    o)
      length_forward_read=${OPTARG}
      ;;
    e)
      length_reverse_read=${OPTARG}
      ;;
    t)
      taxonomy=${OPTARG}
      ;;
    q)
      chimera_ratio=${OPTARG}
      ;;
    k)
      chimera_identity=${OPTARG}
      ;;
    n)
      number_threads=${OPTARG}
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
elif [[ ! -s $taxonomy  ]]
  then
    echo -e "\nERROR: $taxonomy not present or empty\n"
    exit 1
fi

# Convert mapping file to unix format

strings $map > $map"_tmp"
mv $map"_tmp" $map

for number in $(awk '{if(NR>1 && $1!=""){$3}}' $map);
  do
    library_number=$number
    if [[ ! -s $project_name"_"$library_number ]]
      then
        echo -e "\nERROR: $project_name"_"$library_number not present or empty\n"
        echo -e "Add path to libraries if they are not in the current directory\n"
        echo -e "Example: Path_to/libraries_folder/$project_name"
    		exit 1
    fi
done


# Create NG-tax directories

mkdir -p  total_sample_files pre_set_otu_files chimera_clustering_files otus_files all_otus_files clustering_results_files total_sample_distribution_files complementary_tax_files tax_files tmp_databases

# Demultiplexing, chimera checking and OTU picking

for line in $(awk '{if(NR>1 && $1!=""){print $1"separator"$2"separator"$3}}' $map);
  do

    sample_name=$( echo $line | sed 's/separator/\t/g' | awk '{print $1}')
    barcode_sequence=$(echo $line | sed 's/separator/\t/g' | awk '{print $2}')
    library_number=$(echo $line | sed 's/separator/\t/g' | awk '{print $3}')

## Demultiplexing

    awk -v v_library_number_barcode_sequence=$library_number$barcode_sequence '{ \
      if($2==v_library_number_barcode_sequence && length($3)>='$length_forward_read' && length($4)>='$length_reverse_read' ){ \
        print substr($3,1,'$length_forward_read')substr($4,1,'$length_reverse_read')"\t"$1 \
      } \
    }' $project_name"_"$library_number | sort | awk '{print $2;print $1}' | sed 's/@/>/' > "total_sample_files/"$sample_name"_total_sample_file"

## calculating threshold sample.

    total_reads=$( \
      grep -v \^\> "total_sample_files/"$sample_name"_total_sample_file" |  uniq -c | sort -n -r | \
      awk -v v_abundance_threshold=$abundance_threshold '{ \
        if(NR==1){ \
          counter+=$1 \
        } \
        else{ \
          if($1/(counter+$1)>0.01*v_abundance_threshold){ \
            counter+=$1 \
          } \
          else{ \
            exit \
          } \
        } \
      } \
      END{ \
        print counter \
      }' \
      )

## Create a pre-set of OTUs

    grep -v \^\> "total_sample_files/"$sample_name"_total_sample_file" | \
    uniq -c | \
    sort -n -r | \
    awk -v v_sample_name=$sample_name -v v_total_reads=$total_reads '{ \
      counter+=$1; \
      if(counter<=v_total_reads){ \
        print ">"NR"_position_"v_sample_name"_"$1"_"$1/v_total_reads"\n"$2 \
      } \
      else{ \
        exit \
      } \
    }' > "pre_set_otu_files/"$sample_name"_pre_set_OTU_file"

## Generate forward and reverse reads files to check for chimeras

    awk -v v_sample_name=$sample_name '{ \
      if(substr($1,1,1)==">"){ \
        print $1 >"pre_set_otu_files/"v_sample_name"_pre_set_OTU_file_fr"; \
        print $1 >"pre_set_otu_files/"v_sample_name"_pre_set_OTU_file_rr" \
      } \
      else{ \
        print substr($1,1,'$length_forward_read')>"pre_set_otu_files/"v_sample_name"_pre_set_OTU_file_fr"; \
        print substr($1,'$length_forward_read'+1,'$length_reverse_read')>"pre_set_otu_files/"v_sample_name"_pre_set_OTU_file_rr" \
      } \
    }' "pre_set_otu_files/"$sample_name"_pre_set_OTU_file"

## Chimera checking

    $(dirname $0)"/"usearch -usearch_global "pre_set_otu_files/"$sample_name"_pre_set_OTU_file_fr" -strand plus -maxaccepts 0 -maxrejects 0 -db "pre_set_otu_files/"$sample_name"_pre_set_OTU_file_fr" -id 0.97 -uc "chimera_clustering_files/"$sample_name"_chimera_fr" -quiet 2> "chimera_clustering_files/"$sample_name"_chimera_fr_log" & wait $!
    $(dirname $0)"/"usearch -usearch_global "pre_set_otu_files/"$sample_name"_pre_set_OTU_file_rr" -strand plus -maxaccepts 0 -maxrejects 0 -db "pre_set_otu_files/"$sample_name"_pre_set_OTU_file_rr" -id 0.97 -uc "chimera_clustering_files/"$sample_name"_chimera_rr" -quiet 2> "chimera_clustering_files/"$sample_name"_chimera_rr_log" & wait $!

    sed 's/_/\t/g' "chimera_clustering_files/"$sample_name"_chimera_fr" | \
    sort -n -k9 -k14 | \
    awk -v v_chimera_ratio=$chimera_ratio -v v_chimera_identity=$chimera_identity 'BEGIN{ \
      print "chimeric_candidate\torigin_seq\tchimeric_candidate_abundance\torigin_seq_abundance\tratio" \
    } \
    { \
      if($1=="H" && $4>=v_chimera_identity && $9>$14 && $17/$12>v_chimera_ratio && $9>n){ \
        n=$9; \
        print $9"\t"$14"\t"$12"\t"$17"\t"$17/$12 \
      } \
    }' > "chimera_clustering_files/"$sample_name"_chimera_results_fr"

    sed 's/_/\t/g' "chimera_clustering_files/"$sample_name"_chimera_rr" | \
    sort -n -k9 -k14 | \
    awk -v v_chimera_ratio=$chimera_ratio -v v_chimera_identity=$chimera_identity 'BEGIN{ \
      print "chimeric_candidate\torigin_seq\tchimeric_candidate_abundance\torigin_seq_abundance\tratio"
    } \
    { \
      if($1=="H" && $4>=v_chimera_identity && $9>$14 && $17/$12>v_chimera_ratio && $9>n){ \
        n=$9; \
        print $9"\t"$14"\t"$12"\t"$17"\t"$17/$12 \
      } \
    }' > "chimera_clustering_files/"$sample_name"_chimera_results_rr"

    awk '{ \
      if (FNR==1){ \
        number_file++ \
      } \
    } \
    { \
      if (number_file==1){ \
        chimera_fr_array[$1*2-1]=$1; \
        chimera_fr_array[$1*2]=$1 \
      } \
      if (number_file==2){ \
        chimera_rr_array[$1*2-1]=$1; \
        chimera_rr_array[$1*2]=$1 \
      } \
      if (number_file==3){ \
        if ((FNR in chimera_fr_array && FNR in chimera_rr_array) !=1){ \
        print $0 \
        } \
      } \
    }' "chimera_clustering_files/"$sample_name"_chimera_results_fr" "chimera_clustering_files/"$sample_name"_chimera_results_rr" "pre_set_otu_files/"$sample_name"_pre_set_OTU_file" > "otus_files/"$sample_name"_otus"

## Error correction through clustering


    $(dirname $0)/usearch -usearch_global "total_sample_files/"$sample_name"_total_sample_file" -db "pre_set_otu_files/"$sample_name"_pre_set_OTU_file" -strand plus -maxaccepts 0 -maxrejects 0 -top_hit_only -id $clustering_percentage -uc "clustering_results_files/"$sample_name"_result_clustering_missmatch_correction.uc" -quiet 2> "clustering_results_files/"$sample_name"_results_missmatch_correction_log" & wait $!

    awk '{ \
      if (FNR==1){ \
        x++ \
      } \
    } \
    { \
      if (x==1){ \
        valid_otu[$1]="valid" \
      } \
    } \
    { \
      if (x==2){ \
        if($1=="H" && valid_otu[">"$10]=="valid"){ \
          print $10; \
          count+=1 \
        } \
      } \
    } \
    END{ \
      print "0"count \
    }' "otus_files/"$sample_name"_otus" "clustering_results_files/"$sample_name"_result_clustering_missmatch_correction.uc" | \
    sort | \
    uniq -c | \
    sort -n | \
    awk '{ \
      if(NR==1){ \
        count=substr($2,2,length($2)-1) \
      } \
      else{ \
        print $1/count"\t"$1"\t"$2"\t"$2 \
      } \
    }' | \
    sed 's/_position_/\t/' | \
    sed 's/_/\t/' | \
    sed 's/_/\t/' | \
    awk '{ \
      print $1"\t"$2"\t"$6"\t"$5"\t"$3"\t"$4"\t"$7 \
    }' > "total_sample_distribution_files/"$sample_name"_total_sample_distribution"

done

# Generate OTU database

awk '{ \
  if(substr($1,1,1)!=">"){ \
    print $1 \
  } \
}' "otus_files/"*"_otus" | \
sort | \
uniq | \
awk '{ \
  counter+=1; \
  print ">"counter"_OTU\n"$1 \
}' >  "all_otus_files/database_otus"

awk '{ \
  if(substr($1,1,1)==">"){ \
    name=$1 \
  } \
  else{ \
    print name"\n"substr($1,1,'$length_forward_read') > "all_otus_files/database_otus_fr" ;  \
    print name"\n"substr($1,'$length_forward_read'+1,'$length_reverse_read') > "all_otus_files/database_otus_rr" \
  } \
}' "all_otus_files/database_otus"

# Assign taxonomy through clustering.

split -l $(( ($( grep \^\> -c $forward_primer_db ) / 200 + 1)*2 )) -d -a 3 $forward_primer_db "tmp_databases/forward_primer_db_"
split -l $(( ($( grep \^\> -c $reverse_primer_db ) / 200 + 1)*2 )) -d -a 3 $reverse_primer_db "tmp_databases/reverse_primer_db_"

for i in $(seq -w 0 1 199)
  do
    $(dirname $0)/usearch -usearch_global all_otus_files/database_otus_fr -maxaccepts 0 -maxrejects 0 -strand plus -db tmp_databases/forward_primer_db_$i -id 0.90 -uc clustering_results_files/database_results_otus_file_fr_$i.uc -quiet 2> clustering_results_files/database_log_otus_file_fr_$i &
    let counter++
    pids="$pids $!"
    if [ "$counter" -ge "$number_threads" ]
      then
        echo $counter
        counter=0
        echo $pids
        wait $pids
        pids=""
      fi
    $(dirname $0)/usearch -usearch_global all_otus_files/database_otus_rr -maxaccepts 0 -maxrejects 0 -strand plus -db tmp_databases/reverse_primer_db_$i -id 0.90 -uc clustering_results_files/database_results_otus_file_rr_$i.uc -quiet 2> clustering_results_files/database_log_otus_file_rr_$i &
    let counter++
    pids="$pids $!"
    if [ "$counter" -ge "$number_threads" ]
      then
        echo $counter
        counter=0
        echo $pids
        wait $pids
        pids=""
      fi
    if [[ $i == 199 ]]
      then
        echo "hola"
        wait $pids
      fi
done

# Join clustering results

cat "clustering_results_files/database_results_otus_file_fr_"*  > "clustering_results_files/database_results_otus_file_fr.uc"

cat "clustering_results_files/database_results_otus_file_rr_"*  > "clustering_results_files/database_results_otus_file_rr.uc"

rm "clustering_results_files/database_results_otus_file_fr_"* "clustering_results_files/database_results_otus_file_rr_"*

# Generate taxonomic files.

for identity in 90 92 95 97 98 100
  do
    awk -v v_identity=$identity '{ \
      if (FNR==1){ \
        x++ \
      } \
    } \
    { \
      if (x==1){ \
        tax[$1]=$0 \
      } \
    } \
    { \
      if (x==2){ \
        if ($1=="H" && $4>=v_identity){ \
          arr[$9"_"$10]=$10 \
        } \
      } \
    } \
    { \
      if (x==3){ \
        if ($1=="H" && $4>=v_identity){ \
          if ($9"_"$10 in arr != 0){ \
            print $9"\t"tax[$10] \
          } \
        } \
      } \
    }' $taxonomy  "clustering_results_files/database_results_otus_file_fr.uc" "clustering_results_files/database_results_otus_file_rr.uc" | \
    awk '{print $1"\t"$3$4$5$6"\t"$7$8}' | \
    sort | \
    uniq -c | \
    sort -k2 -nk1 | \
    awk -v v_identity=$identity '{ \
      if(NR==1){ \
        m=$1; \
        n=$1; \
        read=$2; \
        tax=$3$4; \
        print_final_line=1 \
      } \
      else{ \
        if(read!=$2){ \
          print ">"read"\t"tax"@"v_identity"@"n"@"n/m; \
          n=$1; \
          m=$1; \
          read=$2; \
          tax=$3$4 \
        } \
        else{ \
          n=$1; \
          m+=$1; \
          tax=$3$4 \
        } \
      } \
    } \
    END{ \
      if(print_final_line==1){ \
        print ">"read"\t"tax"@"v_identity"@"n"@"n/m \
      } \
    }' | \
    LANG=en_EN sort  > "complementary_tax_files/database_otu_genera_"$identity"_tax"

    awk -v v_identity=$identity '{ \
      if (FNR==1){ \
        x++ \
      } \
    } \
    { \
      if (x==1){ \
        tax[$1]=$0 \
      } \
    } \
    { \
      if (x==2){ \
        if ($1=="H" && $4>=v_identity){ \
          arr[$9"_"$10]=$10 \
        } \
      } \
    } \
    { \
      if (x==3){ \
        if ($1=="H" && $4>=v_identity){ \
          if ($9"_"$10 in arr != 0){ \
            print $9"\t"tax[$10] \
          } \
        } \
      } \
    }' $taxonomy  "clustering_results_files/database_results_otus_file_fr.uc" "clustering_results_files/database_results_otus_file_rr.uc" | \
    awk '{print $1"\t"$3$4$5"\t"$6$7}' | \
    sort | \
    uniq -c | \
    sort -k2 -nk1 | \
    awk -v v_identity=$identity '{ \
      if(NR==1){ \
        m=$1; \
        n=$1; \
        read=$2; \
        tax=$3$4; \
        print_final_line=1 \
      } \
      else{ \
        if(read!=$2){ \
          print ">"read"\t"tax"@"v_identity"@"n"@"n/m; \
          n=$1; \
          m=$1; \
          read=$2; \
          tax=$3$4 \
        } \
        else{ \
          n=$1; \
          m+=$1; \
          tax=$3$4 \
        } \
      } \
    } \
    END{ \
      if(print_final_line==1){ \
        print ">"read"\t"tax"@"v_identity"@"n"@"n/m \
      } \
    }' | \
    LANG=en_EN sort > "complementary_tax_files/database_otu_family_"$identity"_tax"
done

awk '{if(substr($1,1,1)==">"){name=$1}else{print name"\t"$1}}' "all_otus_files/database_otus" | LANG=en_EN sort > "complementary_tax_files/database_otu_combined_family_tax"
awk '{if(substr($1,1,1)==">"){name=$1}else{print name"\t"$1}}' "all_otus_files/database_otus" | LANG=en_EN sort > "complementary_tax_files/database_otu_combined_genera_tax"

for identity in 100 98 97 95 92 90; do

LANG=en_EN join -a1  "complementary_tax_files/database_otu_combined_genera_tax" "complementary_tax_files/database_otu_genera_"$identity"_tax" > "complementary_tax_files/database_otu_combined_genera_tax_tmp"

mv "complementary_tax_files/database_otu_combined_genera_tax_tmp" "complementary_tax_files/database_otu_combined_genera_tax"

LANG=en_EN join -a1 "complementary_tax_files/database_otu_combined_family_tax" "complementary_tax_files/database_otu_family_"$identity"_tax" > "complementary_tax_files/database_otu_combined_family_tax_tmp"

mv "complementary_tax_files/database_otu_combined_family_tax_tmp" "complementary_tax_files/database_otu_combined_family_tax"

done

paste <( \
  awk '{ \
    for(i=4;i>2;i--){ \
      if($i==""){ \
        $i="NA;__p;__c;__o;__f;__g@NA@NA@NA" \
      } \
    }; \
    print $1"\t"$2"\t"$3"\t"$4 \
  }' "complementary_tax_files/database_otu_combined_genera_tax" \
) \
<( \
  awk '{ \
    for(i=4;i>2;i--){ \
      if($i==""){ \
        $i="NA;__p;__c;__o;__f;@NA@NA@NA" \
      } \
    }; \
    print $1"\t"$3"\t"$4 \
  }' "complementary_tax_files/database_otu_combined_family_tax" \
  ) | \
  sed 's/;@/;__g@/g' | \
  sed 's/@/\t/g' | \
  sed 's/;/\t/g' | \
  awk '{ \
    if($39<0.5){ \
      $35="__f" \
    }; \
    if($30<0.5){ \
      $26="__f" \
    }; \
    if($9==97){ \
      $17="__g" \
    }; \
    if($9==95){ \
      $8="__g";$16="__f";$17="__g";$35="__f" \
    }; \
    if($9==92){ \
      $7="__f";$8="__g";$16="__f";$17="__g";$26="__f";$35="__f" \
    }; \
    if($9==90){ \
      $7="__f";$8="__g";$26="__f" \
    }; \
      print $1"\t"$2"\t"$3";"$4";"$5";"$6";"$7";"$8"\t"$9"\t"$10"\t"$11"\t"$12";"$13";"$14";"$15";"$16";"$17"\t"$18"\t"$19"\t"$20"\t"$22";"$23";"$24";"$25";"$26";"$27"\t"$28"\t"$29"\t"$30"\t"$31";"$32";"$33";"$34";"$35";"$36"\t"$37"\t"$38"\t"$39 \
    }' | \
    awk '{ \
      if($5>1 || $5=="NA" || $5==$9 || $9=="NA"){ \
        n=0; \
        if($6>=0.5){ \
          m=0 \
        } \
        else{ \
          m=8 \
        } \
      } \
      else{ \
        n=4; \
        if($10>=0.5){ \
          m=0 \
        } \
        else{ \
          m=8 \
        } \
      }; \
      print $1"\t"$(1+n+m+2)"\t"$(2+n+m+2)"\t"$(3+n+m+2)"\t"$(4+n+m+2)"\t"$2 \
    }' > "complementary_tax_files/database_tax_file"


for sample_name in $(awk '{if(NR>1 && $1!=""){print $1}}' $map)
  do
    awk '{ \
      if (FNR==1){ \
        x++ \
      } \
    } \
    { \
      if (x==1){ \
        tax[$6]=$2"\t"$3"\t"$4"\t"$5 \
      } \
    } \
    { \
      if (x==2){ \
        abundance_after_error_correction[">"$7]=$2"_"$1 \
      } \
    } \
    { \
      if (x==3){ \
        if (substr($1,1,1)==">"){ \
          otu=$1"_"abundance_after_error_correction[$1] \
        } \
        else{ \
          print otu"\t"tax[$1]"\n"$1 \
        } \
      }
    }' "complementary_tax_files/database_tax_file" "total_sample_distribution_files/"$sample_name"_total_sample_distribution" "otus_files/"$sample_name"_otus" > "tax_files/"$sample_name"_tax_file"
done

rm -r total_sample_distribution_files tmp_databases
