zcat ../../examples/YRI.chr22.vcf.gz | ../000.scripts/getAF | gzip > YRI.af.txt.gz

#you need to replace the path to homo_sapiens_ancestor_GRCh37_e71
path=/picb/humpopg-bigdata2/zhangrui/AS2_refData/ 
../000.scripts/annoAnc YRI.af.txt.gz ${path}/human_ancestor/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor YRI.anc.txt
gzip YRI.anc.txt
