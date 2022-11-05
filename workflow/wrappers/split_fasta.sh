#!/bin/bash
set +e 
input=$1
output=$2

echo $1

cd $output
echo $?
echo $(ls $output)
echo $?
echo lol
echo $?
zcat $input | head | awk '{print $6}'
echo $?
zcat $input | awk '{print > "lol.txt"}'
echo $?
zcat $input | awk '{print > $6}'
echo $?
#echo $(ls $output)
