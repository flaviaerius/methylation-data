# importante que os awks aqui estao com field separator espaco (default) ao inves de tab, pq eh esse q ta no arquivo assocTSS e q eu mudei qdo juntei as listas assocTSS e gprofiler, que eh a lista "complete_list..", no caso esta esta com field separator tab.

# select promoters

for file in assocTSS_*
do
    awk '{ if ($3 <= 1000 && $3 >= -1000) print $0;}' $file > prom1000_$file
    awk '{ if ($3 <= 500 && $3 >= -500) print $0;}' $file > prom500_$file
done


