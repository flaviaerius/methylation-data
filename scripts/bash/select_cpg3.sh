for file in count*
do
	awk '{if ($1 >=3) print $0;}' $file > cpg3_$file
	mv cpg3_$file cpg3/
done
