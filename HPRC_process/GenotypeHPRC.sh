#!/bin/bash

export REF_CACHE=/liulab/lsong/tmp
for prefix in `cat sample_list.out`
do
	minimap2 ${prefix}.f1v2g.fa.gz kir_genome_seq.fa -x asm5 -t 8 > ${prefix}_mm2.paf
	python3 WholeGenomeGenotypeBwaExon.py ${prefix}_mm2.paf  <(zcat ${prefix}.f1v2g.fa.gz) kir_genome_seq.fa -o $prefix > ${prefix}_bwa_exon_genotype.out
done
