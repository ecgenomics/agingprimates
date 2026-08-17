mkdir genelists

cat owm.analysis.allcaas.filtered.tsv | awk '{print $1}'| grep -v Gene | sort | uniq > genelists/owm.analysis.allcaas.filtered.tsv.genes
cat pla.analysis.allcaas.filtered.tsv | awk '{print $1}'| grep -v Gene | sort | uniq > genelists/pla.analysis.allcaas.filtered.tsv.genes
cat str.analysis.allcaas.filtered.tsv | awk '{print $1}'| grep -v Gene | sort | uniq > genelists/str.analysis.allcaas.filtered.tsv.genes
