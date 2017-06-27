#!/bin/bash

# Usage message

usage="
\n---USAGE---\n
\n
    Description:\n\n
This script make a biom file using a mapping file and a folder containing tax files. \n
\n
    Required inputs\n
\n
-m mapping file containing metadata [txt]\n
-t folder contianing tax files \n
\n
    Output:\n
\n
An OTU database \n
An OTU database [fasta] \n
A tax database \n
A biom file \n\n
Example of usage:\n
$0 -m map_Mocks.txt -t tax_files \n\n
output:\n
map_Mocks_otu_database \n
map_Mocks_otu_database.fa \n
map_Mocks_tax_database \n
map_Mocks.biom \n\n
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

while getopts m:t: opt
  do
    case "$opt" in
    m)
      map=${OPTARG}
      ;;
    t)
	    folder_tax_files=${OPTARG}
      ;;
    *)
      echo -e "\nERROR: invalid arguments used\n";
      exit 1
      ;;
    esac
done

# Check library files presence.

if [[ ! -s $map ]]
	then
		echo -e "\nERROR: $map not present or empty\n"
		exit 1

elif [[ ! -d $folder_tax_files  ]]
	then
		echo -e "\nERROR: $folder_tax_files not present or not a directory\n"
		exit 1
fi

map_name=$(basename $map)

rm -f $map_name"_otu_database" $map_name"_tax_database" "tmp_sparse_otu_table_"$map_name

cat $folder_tax_files"/"*"_tax_file" | awk '{if(substr($1,1,1)==">"){tax=$2}else{print $1"\t"tax}}' | \
LANG=en_EN sort | \
uniq | \
awk '{i+=1;print $0"\t"i-1}' | \
sed 's/;__/\tp__/' | \
sed 's/;__/\tc__/' | \
sed 's/;__/\to__/' | \
sed 's/;__/\tf__/' | \
sed 's/;__/\tg__/' | \
awk -v v_map_name=$map_name '{print $1"\t"$8 > v_map_name"_otu_database" ; print $8"\tk__"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7 > v_map_name"_tax_database"}'


for sample in $(awk '{if(NR>1 && $1!=""){print $1}}' $map)
  do
    sample_position=$((sample_position+1))
    awk  -v v_sample_position=$sample_position '{ \
      if (FNR==1){ \
        x++ \
      }\
    } \
    { \
      if (x==1){ \
        otu_number[$1]=$2 \
      } \
    } \
    { \
      if (x==2){ \
        if (substr($1,1,1)==">"){ \
          otu=$1 \
        } \
        else{ \
          print otu_number[$1]"\t"v_sample_position-1"\t"otu \
        } \
      } \
    }' $map_name"_otu_database" $folder_tax_files"/"$sample"_tax_file" >> "tmp_sparse_otu_table_"$map_name
done

sed 's/_/\t/g' "tmp_sparse_otu_table_"$map_name | \
awk '{ \
	print $1"\t"$2"\t"$8".0" \
}' | \
sort -n -k1,1 -k2,2 > "sparse_otu_table_"$map_name

current_date=$(date | sed 's/ /-/g')

awk -v v_current_date=$current_date 'BEGIN{ \
  print "@@@@{!!!!id!!!!: !!!!None!!!!";
  print "!!!!format!!!!: !!!!Biological Observation Matrix 1.0.0!!!!";
  print "!!!!format_url!!!!: !!!!http://biom-format.org!!!!";
  print "!!!!generated_by!!!!: !!!!NG-Tax v1.1!!!!";
  print "!!!!matrix_element_type!!!!: !!!!float!!!!";
  print "!!!!type!!!!: null";
  print "!!!!matrix_type!!!!: !!!!sparse!!!!";
  print "!!!!date!!!!: !!!!"v_current_date"!!!!"
} \
{ \
  if(FNR==1){ \
    x++ \
  } \
} \
{ \
  if(x==1){ \
    if(NR>1 && $1!=""){ \
      s_number+=1;
      sample_name[s_number]="{!!!!id!!!!: !!!!"$1"!!!!, !!!!metadata!!!!: null}"
    } \
  }\
} \
{ \
  if(x==2){
    o_number+=1;
    taxonomy[$1]="{!!!!id!!!!: !!!!"$1"!!!!, !!!!metadata!!!!: {!!!!taxonomy!!!!: [!!!!"$2"!!!!, !!!!"$3"!!!!, !!!!"$4"!!!!, !!!!"$5"!!!!, !!!!"$6"!!!!, !!!!"$7"!!!!]}}"
  } \
} \
{ \
  if(x==3){
    if(FNR==1){print "!!!!shape!!!!: ["o_number", "s_number"]";
      print "!!!!data!!!!: [["$1","$2","$3"]"
    } \
  else{ \
      print "["$1","$2","$3"]"
    } \
  } \
} \
END{
  print "@@@@], !!!!rows!!!!: [@@@@";
  for(i=0;i<o_number;i++){
    print taxonomy[i]
  };
  print "@@@@],!!!!columns!!!!: [@@@@";
  for(e=1;e<=s_number;e++){
    print sample_name[e]
  };
  print "@@@@]}"
}' $map $map_name"_tax_database" "sparse_otu_table_"$map_name  | \
awk '{
  line=line","$0
}
END{
  print line
}' | \
sed 's/!!!!/\"/g' | \
sed 's/@@@@,//g' | \
sed 's/,@@@@//g' > $map_name".biom"

rm "tmp_sparse_otu_table_"$map_name "sparse_otu_table_"$map_name

awk '{ \
  print ">"$2; \
  print $1 \
}' $map_name"_otu_database" > $map_name"_otu_database.fa"
