This is a README for creating the T1K's reference filtered by OptiType's reference sequence.

1. Run `grep ">" OptiType/data/hla_reference_dna.fasta  | cut -f2 -d' ' > optitype_dna.list` to get the alleles in OptiType.

2. Run `perl filterSeq t1k/hlaidx/hla_dna_seq.fa optitype_dna.list > hla_dna_seq_filt.fa` which only keeps the alleles whose first four digits are in OptiType database.
