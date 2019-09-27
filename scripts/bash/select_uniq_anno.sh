for file in anno*
do
	uniq -u $file > tmp$file
	mv tmp$file $file
done
