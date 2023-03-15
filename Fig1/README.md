### Running commands for T1K, arcasHLA and OptiType on HLA genotyping
T1K:
`
run-t1k -o -1 read1 -2 read2 -f hlaidx/hla_{rna,dna}_seq.fa -t 8 --stage 0 --preset hla
`

arcasHLA:
```
s=sample
arcashla extract sample_starAligned.sortedByCoord.out.bam -o $s -t 8
arcashla genotype $s/${s}.extracted.1.fq.gz $s/${s}.extracted.2.fq.gz -t 8
```

OptiType:
`
OptiTypePipeline.py --input read1 read2 --{rna,dna} --config optitype_config.ini
`
### Creating the whitelist for T1K based on the list extracted from OptiType's code.

1. Run `grep "freq_alleles =" /OptiType/OptiTypePipeline.py | cut -f4 -d\' | tr " " "\n" > optitype_freq_alleles.list` to extract the frequent allele list in OptiType.

2. Run `perl filterSeq t1k/hlaidx/hla_dna_seq.fa optitype_freq_alleles.list > t1k_whitelist.out`. The file t1k_whitelist.out can be used in the "--alleleWhitelist" option in T1K. 

### Creating the T1K's reference filtered by OptiType's reference sequence.

1. Run `grep ">" OptiType/data/hla_reference_dna.fasta  | cut -f2 -d' ' > optitype_dna.list` to get the alleles in OptiType.

2. Run `perl filterSeq t1k/hlaidx/hla_dna_seq.fa optitype_dna.list > hla_dna_seq_filt.fa` which only keeps the alleles whose first four digits are in OptiType database.
