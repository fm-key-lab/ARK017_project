#!/bin/bash

# this script gathers all orinial annotations which have a corresponding prokka annotation

cd /ptmp/iclight/pestis_work/prokka_pestis_pseudotb

ls -1 | while read prokka_gff
    do 
    gff_basename=$(echo $prokka_gff | cut -d'.' -f1)
    grep "ID" ${prokka_gff} | while read line
        do
        query=$(echo ${line} | cut -d' ' -f1,4,5 | sed 's/\t/\.*/g')
        grep -E -m1 "${query}" ../fastas/${gff_basename}*.gff >> ../correponsing_line_gffs/${gff_basename}_correponding_lines.gff
    done
done
