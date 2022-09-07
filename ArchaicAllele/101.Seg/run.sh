../000.scripts/getSeg ind.txt Han_chr22.seg chr22
cat chr22_Neanderthal.txt | ../000.scripts/segFiltMerge | ../000.scripts/con2gvcf | bgzip > Han_Neanderthal.vcf.gz
