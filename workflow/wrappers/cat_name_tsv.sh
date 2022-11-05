#OUTPUT=${@:$#} # last parameter 
#INPUTS=${*%${!#}} # all parameters except the last
INPUTS=$@

#BEGIN {OFS = ","} 
awk 'FNR==1{if (NR==1) print $0, "filename"; next} {print $0, FILENAME}' $INPUTS #> $OUTPUT