for file in prom*assocTSS_*
do
	awk -F "\t" '{print $4}' $file  | uniq -c > count_$file
done
