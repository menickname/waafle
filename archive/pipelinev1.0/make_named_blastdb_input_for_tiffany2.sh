#! /usr/bin/env bash

pangenomes_path="/n/huttenhower_lab/data/repophlan_chocophlan_pangenomes/repophlan_31122013_speciescentroids_ffn"
microbes_list="/n/huttenhower_lab/data/repophlan_chocophlan_pangenomes/single_celled_orgs.tax.txt"
for pangenome in $pangenomes_path/*.ffn
do
    # get species name
    species=$(echo "$pangenome" | perl -pe "s/.*s__(.*?)\..*/\1/")
    fulltaxa=$(grep --max-count=1 "$species" $microbes_list)
    taxaname=$(echo "$fulltaxa" | sed s/'|t__.*$'/''/g | sed s/'|'/'.'/g)
    # read fasta, appending species name
    # *** this syntax guarantees that you read lines, not words ***
    while IFS= read -r line
    do
        if [[ ${line:0:1} == ">" ]]
        then
            echo "$line|$taxaname"
        else
            echo "$line"
        fi
    done < "$pangenome"
done
