#!/bin/bash

# Usage message

usage="
\n---USAGE---\n
\n
    Description:\n\n
This script assigns taxonomy to a fasta file, containing concatenated pair-end reads. The headers should start \n
with a number, for example, >1_sequence_name. \n
\n
    Required inputs\n
\n
-i File to be assigned [fasta]\n
-f NG-Tax customized forward-primer database name\n
-r NG-Tax customized reverse-primer database name\n
-o forward-read length\n
-e reverse-read length\n
-t taxonomic tablefile\n
-n number of threads\n
\n
    Output:\n
A folder containing the OTUs file \n
A folder containing the clustering result files \n
A folder containing the complementary taxonomic assignment files\n
A folder containing the final taxonomic assignment file \n
\n
Example of usage:\n
$0 -i Mock_theoretical -f primer_F515_71_nt_1mm_db -r primer_R806_70_nt_1mm_db -o 71 -e 70 -t Silva_111_taxa_map_RDP_6_levels_full.txt -n 24 \n\n
output:\n
otus_files_Mock_theoretical \n
clustering_results_files_Mock_theoretical \n
complementary_tax_file_Mock_theoretical \n
tax_files_Mock_theoretical \n\n
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

while getopts i:f:r:o:e:t:n: opt
  do
    case "$opt" in
    i)
      sample_name=${OPTARG}
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

if [[ ! -s $sample_name ]]
	then
		echo -e "\nERROR: $sample_name not present or empty\n"
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

# Create NG-tax directories

mkdir -p  "otus_files_"$sample_name "clustering_results_files_"$sample_name  "complementary_tax_files_"$sample_name "tax_files_"$sample_name tmp_databases

# Create OTU files

cat $sample_name | \
awk -v v_sample_name=$sample_name '{ \
  if(substr($1,1,1)==">"){ \
    name=$1 \
  } \
  else{ \
    print name"\n"substr($1,1,'$length_forward_read') > "otus_files_"v_sample_name"/"v_sample_name"_otus_fr" ; \
    print name"\n"substr($1,'$length_forward_read'+1,'$length_reverse_read') > "otus_files_"v_sample_name"/"v_sample_name"_otus_rr" \
  } \
}'


# Assign taxonomy through clustering.

split -l $(( ($( grep \^\> -c $forward_primer_db ) / 100 + 1)*2 )) -d -a 2 $forward_primer_db "tmp_databases/forward_primer_db_"
split -l $(( ($( grep \^\> -c $reverse_primer_db ) / 100 + 1)*2 )) -d -a 2 $reverse_primer_db "tmp_databases/reverse_primer_db_"

rm -f clustering_commands

for i in $(seq -w 0 1 99)
  do
    echo "/home/jramirogarcia/Work/NG-Tax/usearch -usearch_global otus_files_"$sample_name"/"$sample_name"_otus_fr -maxaccepts 0 -maxrejects 0 -strand plus -db tmp_databases/forward_primer_db_"$i" -id 0.90 -uc clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_fr_"$i".uc -quiet 2> clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_fr_log_"$i" &" >> clustering_commands
    let "counter++"
    if [ "$counter" -ge "$number_threads" ]
      then
        counter=0
        echo "wait" >> clustering_commands
      fi
    echo "/home/jramirogarcia/Work/NG-Tax/usearch -usearch_global otus_files_"$sample_name"/"$sample_name"_otus_rr -maxaccepts 0 -maxrejects 0 -strand plus -db tmp_databases/reverse_primer_db_"$i" -id 0.90 -uc clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_rr_"$i".uc -quiet 2> clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_rr_log_"$i" &" >> clustering_commands
    let "counter++"
    if [ "$counter" -ge "$number_threads" ]
      then
        counter=0
        echo "wait" >> clustering_commands
      fi
done

echo "wait" >> clustering_commands

sh clustering_commands

rm clustering_commands

cat "clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_fr_"*  > "clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_fr.uc"

cat "clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_rr_"*  > "clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_rr.uc"

rm "clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_fr_"* "clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_rr_"*


# Generate taxonomic files

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
    }' $taxonomy  "clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_fr.uc" "clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_rr.uc" | \
    awk '{ \
      print $1"\t"$3$4$5$6"\t"$7$8 \
    }' | \
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
      }
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
    }' | LANG=en_EN sort  > "complementary_tax_files_"$sample_name"/"$sample_name"_otu_genera_"$identity"_tax"

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
}' $taxonomy  "clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_fr.uc" "clustering_results_files_"$sample_name"/"$sample_name"_results_otus_file_rr.uc" | \
awk '{\
  print $1"\t"$3$4$5"\t"$6$7 \
}' | \
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
  }\
}' | LANG=en_EN sort > "complementary_tax_files_"$sample_name"/"$sample_name"_otu_family_"$identity"_tax"


done

awk '{ \
  if(substr($1,1,1)==">"){ \
    name=$1 \
  } \
  else{ \
    print name"\t"$1 \
  } \
}' $sample_name | \
LANG=en_EN sort > "complementary_tax_files_"$sample_name"/"$sample_name"_otu_combined_family_tax"

awk '{ \
  if(substr($1,1,1)==">"){ \
    name=$1 \
  } \
  else{ \
    print name"\t"$1 \
  } \
}' $sample_name | \
LANG=en_EN sort > "complementary_tax_files_"$sample_name"/"$sample_name"_otu_combined_genera_tax"

for identity in 100 98 97 95 92 90
  do
    LANG=en_EN join -a1  "complementary_tax_files_"$sample_name"/"$sample_name"_otu_combined_genera_tax" "complementary_tax_files_"$sample_name"/"$sample_name"_otu_genera_"$identity"_tax" > "complementary_tax_files_"$sample_name"/"$sample_name"_otu_combined_genera_tax_tmp"
    mv "complementary_tax_files_"$sample_name"/"$sample_name"_otu_combined_genera_tax_tmp" "complementary_tax_files_"$sample_name"/"$sample_name"_otu_combined_genera_tax"
    LANG=en_EN join -a1 "complementary_tax_files_"$sample_name"/"$sample_name"_otu_combined_family_tax" "complementary_tax_files_"$sample_name"/"$sample_name"_otu_family_"$identity"_tax" > "complementary_tax_files_"$sample_name"/"$sample_name"_otu_combined_family_tax_tmp"
    mv "complementary_tax_files_"$sample_name"/"$sample_name"_otu_combined_family_tax_tmp" "complementary_tax_files_"$sample_name"/"$sample_name"_otu_combined_family_tax"

done

paste <( \
  awk '{ \
    for(i=4;i>2;i--){ \
      if($i==""){ \
        $i="NA;__p;__c;__o;__f;__g@NA@NA@NA" \
      } \
    }; \
    print $1"\t"$2"\t"$3"\t"$4 \
  }' "complementary_tax_files_"$sample_name"/"$sample_name"_otu_combined_genera_tax" \
) \
<( \
awk '{ \
  for(i=4;i>2;i--){ \
    if($i==""){ \
      $i="NA;__p;__c;__o;__f;@NA@NA@NA" \
    } \
  }; \
  print $1"\t"$3"\t"$4 \
}' "complementary_tax_files_"$sample_name"/"$sample_name"_otu_combined_family_tax" \
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
    $8="__g"; \
    $16="__f"; \
    $17="__g"; \
    $35="__f" \
  }; \
  if($9==92){ \
    $7="__f"; \
    $8="__g"; \
    $16="__f"; \
    $17="__g"; \
    $26="__f"; \
    $35="__f" \
  }; \
  if($9==90){ \
    $7="__f"; \
    $8="__g"; \
    $26="__f" \
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
  print $1"\t"$(1+n+m+2)"\t"$(2+n+m+2)"\t"$(3+n+m+2)"\t"$(4+n+m+2)"\n"$2 \
}' > "tax_files_"$sample_name"/"$sample_name"_tax_file"

rm -r -f tmp_databases
