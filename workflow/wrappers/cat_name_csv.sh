OUTPUT=${@:$#} # last parameter 
INPUTS=${*%${!#}} # all parameters except the last

echo "last" $OUTPUT
echo "other" $INPUTS

awk 'BEGIN {OFS = ","} FNR==1{if (NR==1) print "filename", $0; next} {print FILENAME, $0}' $INPUTS > $OUTPUT