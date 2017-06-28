#!/bin/bash

# Usage message

usage="
\n---USAGE---\n
\n
    Description:\n\n
This script generates NG-Tax customized databases using the primer sequences and the read lengths introduced by de user.\n
Degenerate positions in the primers should be included between brackets. A list of allowed primer sequences containing \n
mismatches can be provided as optional input.\n
\n
    Required inputs\n
\n
-d Reference 16S rRNA database file [fasta]\n
-t SILVA taxonomy\n
-k forward-primer sequence, degenerate position between brackets\n
-p reverse-primer sequence, degenerate position between brackets\n
-q complementary reversed reverse-primer sequence, degenerate position between brackets\n
-f NG-Tax customized forward-primer database name\n
-r NG-Tax customized reverse-primer database name\n
-x NG-Tax customized taxonomy table name\n
-o forward-read length\n
-e reverse-read length\n
\n
    Optional input:\n
\n
-y allowed mismatching forward-primer sequences list\n
-z allowed mismatching reverse-primer sequences list\n
\n
    Output:\n
\n
A database with all sequences that have matching (or allowed mismatching) forward primer \n
A database with all sequences that have matching (or allowed mismatching) reverse primer \n
A taxonomy file customized for NG-Tax\n\n
Example of usage:\n
$0 -d SILVA_128_SSURef_tax_silva.fasta -t tax_slv_ssu_128.txt -k GTGCCAGC[AC]GCCGCGGTAA -p GGACTAC[ACT][ACG]GGGT[AT]TCTAAT -q ATTAGA[AT]ACCC[TCG][ATG]GTAGTCC -f primer_F515_71_nt_1mm_db -r primer_R806_70_nt_1mm_db -x Silva_taxonomic_table -o 71 -e 70 -y primer_F515_1mm -z primer_R806_1mm\n\n
output:\n
primer_F515_71_nt_1mm_db\n
primer_R806_70_nt_1mm_db\n
Silva_taxonomic_table\n\n
--- End of USAGE ---\n
"

# Help option, check number of argument.

if [[ "$1" == "-h" ]]
  then
    echo -e $usage
    exit 0
elif [[ $# -ne 20 && $# -ne 24 ]]
  then
    echo -e "ERROR: Invalid number of args \nTry $0 -h for help"
    exit 1
fi

# Assign arguments to variables, check for invalid arguments.

while getopts d:t:k:p:q:f:r:x:o:e:y:z: opt
  do
  case "$opt" in
    d)
      database_address=${OPTARG}
      ;;
    t)
      taxonomy_address=${OPTARG}
      ;;
    k)
      forward_primer_sequence=${OPTARG}
      ;;
    p)
      reverse_primer_sequence=${OPTARG}
      ;;
    q)
      complementary_reversed_reverse_primer_sequence=${OPTARG}
      ;;
    f)
      forward_primer_database_name=${OPTARG}
      ;;
    r)
      reverse_primer_database_name=${OPTARG}
      ;;
    x)
      taxonomy_name=${OPTARG}
      ;;
    o)
      length_forward_read=${OPTARG}
      ;;
    e)
      length_reverse_read=${OPTARG}
      ;;
    y)
      allowed_forward_primers=${OPTARG}
      ;;
    z)
      allowed_reverse_primers=${OPTARG}
      ;;
    *)
      echo -e "\nERROR: invalid arguments used\n";
      exit 1
      ;;
    esac
done

# Check reference 16S rRNA database file presence. If provided, check also allowed primer files presence.

if [[ ! -s $database_address ]]
	then
		echo -e "\nERROR: $database_address not present or empty\n"
		exit 1

elif [[ -n $allowed_forward_primers && ! -s $allowed_forward_primers  ]]
	then
		echo -e "\nERROR: $allowed_forward_primers not present or empty\n"
		exit 1

elif [[ -n $allowed_reverse_primers && ! -s $allowed_reverse_primers  ]]
	then
		echo -e "\nERROR: $allowed_reverse_primers not present or empty\n"
		exit 1
fi


# Generate a temporary database containing complementary reverse sequences.

awk '{ \
  if(substr($1,1,1)==">"){ \
    if(NR==1){ \
      print $0 \
    } \
    else{ \
      print "\n"$0 \
    } \
  } \
  else{ \
    printf("%s",$0) \
  } \
} \
END{ \
  printf("\n") \
}' $database_address | \
awk '{ \
  if(substr($1,1,1)==">"){ \
    print $1 \
  } \
  else{ \
    gsub("U","T",$1); \
    for(i=length($1);i!=0;i--){ \
      complementary_reversed_sequence=complementary_reversed_sequence substr($1,i,1) \
    }; \
    gsub("T","1",complementary_reversed_sequence); \
    gsub("G","2",complementary_reversed_sequence); \
    gsub("C","3",complementary_reversed_sequence); \
    gsub("A","4",complementary_reversed_sequence); \
    gsub("1","A",complementary_reversed_sequence); \
    gsub("2","C",complementary_reversed_sequence); \
    gsub("3","G",complementary_reversed_sequence); \
    gsub("4","T",complementary_reversed_sequence); \
    print $1"\t"complementary_reversed_sequence; \
    complementary_reversed_sequence="" \
  } \
}' > "tmp_"$forward_primer_database_name"_"$reverse_primer_database_name


# Assign length to the NG-Tax customized databases, 10 nucleotides are added to account for insertions.

extra_length_forward_read=$(echo $(($length_forward_read + 10)))
extra_length_reverse_read=$(echo $(($length_reverse_read + 10)))


# Retrieve sequences with matching forward primer from the reference 16S rRNA database to built the NG-Tax customized forward primer database.

grep -e $forward_primer_sequence'[A-Z]\{'$extra_length_forward_read'\}' -e \>.*  -o "tmp_"$forward_primer_database_name"_"$reverse_primer_database_name | \
awk '{ \
  if(substr($1,1,1)==">"){ \
    p=1;name=$1 \
  } \
  else{ \
    if(p==1){ \
      print name"\n"$1; \
		  p=0 \
    } \
  } \
}' | \
sed 's/'$forward_primer_sequence'//' > $forward_primer_database_name &

# Retrieve sequences with matching reverse primer from the reference 16S rRNA database to built the NG-Tax customized reverse primer database.

grep -e $reverse_primer_sequence'[A-Z]\{'$extra_length_reverse_read'\}' -e \>.*  -o "tmp_"$forward_primer_database_name"_"$reverse_primer_database_name | \
awk '{ \
  if(substr($1,1,1)==">"){ \
    p=1; \
    name=$1 \
  } \
  else{ \
    if(p==1){ \
      print name"\n"$1; \
      p=0 \
    } \
  } \
}' | \
sed 's/'$reverse_primer_sequence'//' > $reverse_primer_database_name &

wait

# If a list of allowed mismatching forward primers is provided: Add to the NG-Tax customized forward primer database those sequences with allowed mismatching forward primers.

if [[ "$allowed_forward_primers" != "" ]]

	then

    bname_allowed_forward_primers=$(basename "$allowed_forward_primers")
		strings $allowed_forward_primers > $bname_allowed_forward_primers"_non_return_carriage"

		awk  '{ \
      if (FNR==1){ \
        file_number++ \
      } \
    } \
		{ \
      if(file_number==1){ \
        if (substr($1,1,1)==">"){ \
          name[$1]="present" \
        } \
      } \
    } \
    { \
      if(file_number==2){ \
        if (substr($1,1,1)==">"){ \
          if(name[$1]!="present"){ \
            print $1; \
            p=1 \
          } \
        } \
        else{ \
          if(p==1){ \
            print; \
            p=0 \
          } \
        } \
      } \
    }' \
		$forward_primer_database_name \
		"tmp_"$forward_primer_database_name"_"$reverse_primer_database_name > "db_no_present_in_"$forward_primer_database_name

		rm -f "sequences_allowed_by_mismatch_"$bname_allowed_forward_primers"_"$forward_primer_database_name

		for forward_primer_to_check in  `awk '{if($1!=""){print $1}}' $bname_allowed_forward_primers"_non_return_carriage"`
    do \
      grep -e $forward_primer_to_check'[A-Z]\{'$extra_length_forward_read'\}' -e \>.*  -o  "db_no_present_in_"$forward_primer_database_name | \
			awk -v v_extra_length_forward_read=$extra_length_forward_read '{ \
        if(substr($1,1,1)==">"){ \
          p=1; \
          name=$1 \
        } \
        else{ \
          if(p=1){ \
            print name"\t"substr($1,1 + length($1) - v_extra_length_forward_read,v_extra_length_forward_read); \
            p=0 \
          } \
        } \
      }'  >> "sequences_allowed_by_mismatch_"$bname_allowed_forward_primers"_"$forward_primer_database_name
    done

	  cat $forward_primer_database_name <( \
      LANG=en_EN sort "sequences_allowed_by_mismatch_"$bname_allowed_forward_primers"_"$forward_primer_database_name | \
		  awk '{ \
        if(name!=$1){ \
          print $1"\n"$2; \
          name=$1 \
        } \
      }' \
    ) > "allowed_mismatches_"$forward_primer_database_name

    cat "allowed_mismatches_"$forward_primer_database_name > $forward_primer_database_name
    rm -f $bname_allowed_forward_primers"_non_return_carriage"
    rm -f "allowed_mismatches_"$forward_primer_database_name
    rm -f "sequences_allowed_by_mismatch_"$bname_allowed_forward_primers"_"$forward_primer_database_name
    rm -f "db_no_present_in_"$forward_primer_database_name
fi

# If a list of allowed mismatching reverse primers is provided: Add to the NG-Tax customized reverse primer database those sequences with allowed mismatching reverse primers.

if [[ "$allowed_reverse_primers" != "" ]]

	then

    bname_allowed_reverse_primers=$(basename "$allowed_reverse_primers")
		strings $allowed_reverse_primers > $bname_allowed_reverse_primers"_non_return_carriage"

		awk  '{ \
      if (FNR==1){ \
        file_number++ \
      } \
    } \
    { \
      if(file_number==1){ \
        if (substr($1,1,1)==">"){ \
          name[$1]="present" \
        } \
      } \
    } \
    { \
      if(file_number==2){ \
        if (substr($1,1,1)==">"){ \
          if(name[$1]!="present"){ \
            print $1; \
            p=1 \
          } \
        } \
        else{ \
          if(p==1){ \
            print; \
            p=0 \
          } \
        } \
      } \
    }' $reverse_primer_database_name "tmp_"$forward_primer_database_name"_"$reverse_primer_database_name > "db_no_present_in_"$reverse_primer_database_name

		rm -f "sequences_allowed_by_mismatch_"$bname_allowed_reverse_primers"_"$reverse_primer_database_name

		for reverse_primer_to_check in `awk '{if($1!=""){print $1}}' $bname_allowed_reverse_primers"_non_return_carriage"`
    do \
			grep -e $reverse_primer_to_check'[A-Z]\{'$extra_length_reverse_read'\}' -e \>.*  -o  "db_no_present_in_"$reverse_primer_database_name | \
			awk -v v_extra_length_reverse_read=$extra_length_reverse_read '{ \
        if(substr($1,1,1)==">"){ \
          p=1; \
          name=$0 \
        } \
        else{ \
          if(p==1){ \
            print name"\t"substr($1,1 + length($1) - v_extra_length_reverse_read,v_extra_length_reverse_read); \
            p=0 \
          } \
        } \
      }' >> "sequences_allowed_by_mismatch_"$bname_allowed_reverse_primers"_"$reverse_primer_database_name
		done

    cat $reverse_primer_database_name <( \
      LANG=en_EN sort "sequences_allowed_by_mismatch_"$bname_allowed_reverse_primers"_"$reverse_primer_database_name | \
      awk '{ \
        if(name!=$1){ \
          print $1"\n"$2; \
          name=$1 \
        } \
      }' \
    ) > "allowed_mismatches_"$reverse_primer_database_name

    cat "allowed_mismatches_"$reverse_primer_database_name > $reverse_primer_database_name
    rm -f $bname_allowed_reverse_primers"_non_return_carriage"
    rm -f "allowed_mismatches_"$reverse_primer_database_name
    rm -f "sequences_allowed_by_mismatch_"$bname_allowed_reverse_primers"_"$reverse_primer_database_name
    rm -f "db_no_present_in_"$reverse_primer_database_name
fi

wait

# Remove temporary database containing complementary reverse sequences and extra files

rm "tmp_"$forward_primer_database_name"_"$reverse_primer_database_name

# Create NG-Tax taxonomy file

awk -F "\t" '{ \
  if(FNR==1){ \
    x++ \
  } \
} \
{ \
  if(x==1){ \
    if(substr($1,1,7)=="Archaea" || substr($1,1,8)=="Bacteria"){ \
      gsub(" ","_",$1); \
      array_classification[$1]=$1 \
    } \
  } \
} \
{ \
 if(x==2){ \
  if(substr($1,1,1)==">"){ \
    sub(">","",$1); \
    sub(" ","\t"); \
    gsub(" ","_"); \
    n=split($2,array_names,";"); \
    name=array_names[1]";"; \
    for(i=2;i<n;++i){ \
      name=name""array_names[i]";" \
    }; \
    last_name=name""array_names[n]";"; \
    if(array_classification[last_name]!=""){ \
      classification=array_classification[last_name] \
    } \
    else{ \
      classification=array_classification[name] \
    }; \
    if(classification!=""){ \
      split(classification,classification_levels,";"); \
      if(classification_levels[2]==""){ \
        classification_levels[2]="p" \
      }; \
      if(classification_levels[3]==""){ \
        classification_levels[3]="c" \
      }; \
      if(classification_levels[4]==""){ \
        classification_levels[4]="o" \
      }; \
      if(classification_levels[5]==""){ \
        classification_levels[5]="f" \
      }; \
      if(classification_levels[6]==""){ \
        classification_levels[6]="g" \
      }; \
      print $1"\t"classification_levels[1]";\t__"classification_levels[2]";\t__"classification_levels[3]";\t__"classification_levels[4]";\t__"classification_levels[5]";\t__"classification_levels[6]}}}}' $taxonomy_address $database_address > $taxonomy_name
