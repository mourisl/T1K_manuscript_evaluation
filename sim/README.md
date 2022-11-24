This part is for how to generate the simulated data in T1K. We used mason v0.1.2 to generate the simulated reads. The kir_rna_seq.fa is the reference file for T1K's genotpying on RNA-seq data, which can be found in the kiridx_2_10_0.tar.gz file. The overall procedure is as following:

1.Randomly select KIR genes and alleles for 100 simulated samples:

	perl GeneratePairedTestHaplotypes.pl kir_rna_seq.fa 100 > random_paired_haplotype_100.tsv

2. Generate the simulated sequences:

	perl BuildExploreSim.pl kir_rna_seq.fa random_paired_haplotype_100.tsv

It will generate files like sample_[0..99]_{1,2}.fq.

3. Run T1K with commands like ```run-t1k -1 sample_0_1.fq -2 sample_0_2.fq -o sample_0 -t 8 -f kir_rnaseq.fa```
4. Compare T1K's results with simulated data in the T1K folder

	perl ../EvaluateSim.pl ../random_paired_haplotype_100.tsv
