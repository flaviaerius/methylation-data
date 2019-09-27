for file in anno*
do
	sort -d $file > srt_$file
	mv srt_$file $file
done
