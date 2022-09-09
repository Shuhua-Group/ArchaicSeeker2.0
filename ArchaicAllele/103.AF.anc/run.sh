zcat ../../examples/YRI.chr22.vcf.gz | ../000.scripts/getAF | gzip > YRI.af.txt.gz

../000.scripts/annoAnc YRI.af.txt.gz ../../examples/homo_sapiens_ancestor 22 22 YRI.anc.txt
gzip YRI.anc.txt
