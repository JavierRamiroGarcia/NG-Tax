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
-k forward-primer sequence, degenerate position between brackets\n
-p reverse-primer sequence, degenerate position between brackets\n
-q complementary reversed reverse-primer sequence, degenerate position between brackets\n
-f NG-Tax customized forward-primer database name\n
-r NG-Tax customized reverse-primer database name\n
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
A database with all sequences that have matching (or allowed mismatching) reverse primer \n\n
Example of usage:\n
$0 -d Silva_111_full_unique.fasta -k GTGCCAGC[AC]GCCGCGGTAA -p GGACTAC[ACT][ACG]GGGT[AT]TCTAAT -q ATTAGA[AT]ACCC[TCG][ATG]GTAGTCC -f primer_F515_71_nt_1mm_db -r primer_R806_70_nt_1mm_db -o 71 -e 70 -y primer_F515_1mm -z primer_R806_1mm\n\n
output:\n
primer_F515_71_nt_1mm_db\n
primer_R806_70_nt_1mm_db\n\n
--- End of USAGE ---\n
"

# Help option, check number of argument.

if [[ "$1" == "-h" ]]
  then
    echo -e $usage
    exit 0
elif [[ $# -ne 16 && $# -ne 20 ]]
  then
    echo -e "ERROR: Invalid number of args \nTry $0 -h for help"
    exit 1
fi

# Assign arguments to variables, check for invalid arguments.

while getopts d:k:p:q:f:r:o:e:y:z: opt
  do
  case "$opt" in
    d)
      database_address=${OPTARG}
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
    print $1 \
  } \
  else{ \
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
}' \
$database_address > "tmp_"$forward_primer_database_name"_"$reverse_primer_database_name


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

		strings $allowed_forward_primers > "tmp_"$allowed_forward_primers
    mv "tmp_"$allowed_forward_primers $allowed_forward_primers

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

		rm -f "sequences_allowed_by_mismatch_"$allowed_forward_primers"_"$forward_primer_database_name

		for forward_primer_to_check in `awk '{if($1!=""){print $1}}' $allowed_forward_primers`
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
      }'  >> "sequences_allowed_by_mismatch_"$allowed_forward_primers"_"$forward_primer_database_name
    done

	  cat $forward_primer_database_name <( \
      LANG=en_EN sort "sequences_allowed_by_mismatch_"$allowed_forward_primers"_"$forward_primer_database_name | \
		  awk '{ \
        if(name!=$1){ \
          print $1"\n"$2; \
          name=$1 \
        } \
      }' \
    ) > "allowed_mismatches_"$forward_primer_database_name

    cat "allowed_mismatches_"$forward_primer_database_name > $forward_primer_database_name
    rm -f "allowed_mismatches_"$forward_primer_database_name
    rm -f "sequences_allowed_by_mismatch_"$allowed_forward_primers"_"$forward_primer_database_name
    rm -f "db_no_present_in_"$forward_primer_database_name
fi

# If a list of allowed mismatching reverse primers is provided: Add to the NG-Tax customized reverse primer database those sequences with allowed mismatching reverse primers.

if [[ "$allowed_reverse_primers" != "" ]]

	then

		strings $allowed_reverse_primers > "tmp_"$allowed_reverse_primers
    mv "tmp_"$allowed_reverse_primers $allowed_reverse_primers

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

		rm -f "sequences_allowed_by_mismatch_"$allowed_reverse_primers"_"$reverse_primer_database_name

		for reverse_primer_to_check in `awk '{if($1!=""){print $1}}' $allowed_reverse_primers`
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
      }' >> "sequences_allowed_by_mismatch_"$allowed_reverse_primers"_"$reverse_primer_database_name
		done

    cat $reverse_primer_database_name <( \
      LANG=en_EN sort "sequences_allowed_by_mismatch_"$allowed_reverse_primers"_"$reverse_primer_database_name | \
      awk '{ \
        if(name!=$1){ \
          print $1"\n"$2; \
          name=$1 \
        } \
      }' \
    ) > "allowed_mismatches_"$reverse_primer_database_name

    cat "allowed_mismatches_"$reverse_primer_database_name > $reverse_primer_database_name
    rm -f "allowed_mismatches_"$reverse_primer_database_name
    rm -f "sequences_allowed_by_mismatch_"$allowed_reverse_primers"_"$reverse_primer_database_name
    rm -f "db_no_present_in_"$reverse_primer_database_name
fi

wait

# Remove temporary database containing complementary reverse sequences and extra files

rm "tmp_"$forward_primer_database_name"_"$reverse_primer_database_name
