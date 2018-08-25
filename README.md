# mitotic-bookmarking-CAR

 These in-house scrtipts are associted with mitotic-bookmarking-CAR paper .

 filtering_distance.R reads genecount data genereated by iontorrent RNA-seq plugin, and filter low abundant genes by keeping genes that rpm>10 in at least 4 samples. Then sample clustering were performed by R package "cluster" and "factoextra".

main.R is the major R script containing calling CARs from A549 and Hela cells; Functinal annotation;coverage analysis; chip-seq peak intergration; gglpot2 for a variaty of plots used in manuscript. 
anno_from_HGNC.pl is a perl script to annotate gene based on HGNC and genecode V19 datasets.
enrichmentTest.R is for statistic test on differenct group of CAR/genes enrichment in varient dataset.
filt_MARGI_distal.pl is for keeping only distal and inter-chromosal RNA-DNA interaction link data.
