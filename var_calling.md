# Table of contents

- [SH1000](#sh1000)
 * [Check the reference](#check-the-reference)
 * [Map and call variants](#map-and-call-variants)
 * [SH1000 variants](#sh1000-variants)
- [ATCC 6538](#atcc-6538)
 * [Assemble with `a5miseq` and annotate with `prokka`](#assemble-with-a5miseq-and-annotate-with-prokka)
 * [Map and call variants](#map-and-call-variants)
- [CC398](#cc398)
 * [Assemble with `a5miseq` and annotate with `prokka`](#assemble-with-a5miseq-and-annotate-with-prokka)
 * [Map and call variants](#map-and-call-variants)

# SH1000

## Check the reference

In theory there should not be any SNPs between the ancestor and the reference since they are both derived from the same stock of the same strain but check anyway:

Align reads from ancestor and calls SNPs:

```sh
snippy --cpus 12 --outdir ref_fix_sh1000 --ref sh1000_pol_6.gbk --rgid 53_S3 --prefix 53_S3 --pe1 ../53_S3_L001_R1_001.fastq.gz --pe2 ../53_S3_L001_R2_001.fastq.gz &> 53_S3.log
```

So there are actually a few SNPs. The usual ribosomal operons etc.:

|TYPE|NT_POS  |AA_POS | EFFECT                                |GENE  |PRODUCT                       |
|----|--------|-------|---------------------------------------|------|------------------------------|
|snp |100/882 |34/293 | missense_variant c.100G>T p.Ala34Ser  |*nanA*|N-acetylneuraminate lyase     |
|snp |697/1389|233/462| missense_variant c.697G>A p.Val233Ile |*tcyP*|L-cystine uptake protein TcyP |
|snp |        |       | non_coding_transcript_variant         |      |16S ribosomal RNA             |
|snp |184/282 | 62/93 | missense_variant c.184T>A p.Ser62Thr  |      |hypothetical protein          |

Table 1. SNPs between the SH1000 reference genome and the particular SH1000 strain used here.

Alter the reference to reflect the differences and re-annotate the new reference:

```sh
bcftools consensus -f   53_S3/53_S3.vcf.gz  -o sh1000_dueppel.fna
prokka --cpus 12 sh1000_dueppel.fna --genus Staphylococcus --prefix sh1000_dueppel &> sh1000_dueppel.log
```

## Map and call variants

```sh
for i in `cut -f 5 sh1000.tsv | grep -v suffix`
do
  snippy --cpus 12 --outdir $i --ref sh1000_dueppel/sh1000_dueppel.gbk --rgid $i --prefix $i --pe1 ..//${i}_L001_R1_001.fastq.gz --pe2 ../${i}_L001_R2_001.fastq.gz &> $i.log
done
```

Merge the vcf files:

```sh
for i in `cut -f 5 sh1000.tsv | grep -v suffix`
do
  cp sh1000/${i}/${i}.vcf.gz* sh1000/vcf
done
vcf-merge *vcf.gz > sh1000.vcf
```

## SH1000 variants

You can download this as a csv [here](https://github.com/Perugolate/biocides/blob/master/sh1000/sh1000_variants.csv).

|Treatment|  Type        | Label   |  Mutation                                                                                       |  Locus tag   | Annotation|  Function                                        |
|---------|--------------|----------|-------------------------------------------------------------------------------------------------|--------------|-----------|--------------------------------------------------|
|  PG     |  population  |    4     |  frameshift variant c.117delG p.Gly41fs                                                         |  PROKKA_01659|   *hemY*  |  Protoporphyrinogen oxidase                      | 
|  PG     |  colony      | 4-s-c1   |  frameshift variant c.117delG p.Gly41fs                                                         |  PROKKA_01659|   *hemY*  |  Protoporphyrinogen oxidase                      | 
|  PG     |  population  |   32     |  frameshift variant c.117delG p.Gly41fs                                                         |  PROKKA_01659|   *hemY*  |  Protoporphyrinogen oxidase                      | 
|  PG     |  colony      | 32-s-c1  |  frameshift variant c.117delG p.Gly41fs                                                         |  PROKKA_01659|   *hemY*  |  Protoporphyrinogen oxidase                      | 
|  PG     |  population  |    4     |  stop gained c.991A>T p.Lys331\*                                                                |  PROKKA_01871|   *cls2*  |  Cardiolipin synthase                            | 
|  PG     |  colony      | 4-s-c1   |  stop gained c.991A>T p.Lys331\*                                                                |  PROKKA_01871|   *cls2*  |  Cardiolipin synthase                            | 
|  PG     |  population  |   32     |  frameshift variant c.786delA p.Lys262fs                                                        |  PROKKA_01871|   *cls2*  |  Cardiolipin synthase                            | 
|  PG     |  colony      | 32-s-c1  |  frameshift variant c.786delA p.Lys262fs                                                        |  PROKKA_01871|   *cls2*  |  Cardiolipin synthase                            | 
|  PG     |  colony      | 32-s-c1  |  missense variant c.515C>T p.Pro172Leu                                                          |  PROKKA_01192|   *parE*  |  DNA topoisomerase 4 subunit B                   | 
|  PG     |  population  |   33     |  stop lost c.1485A>T p.Ter495Tyrext\*?                                                          |  PROKKA_01871|   *cls2*  |  Cardiolipin synthase                            | 
|  PG     |  colony      | 33-s-c1  |  stop lost c.1485A>T p.Ter495Tyrext\*?                                                          |  PROKKA_01871|   *cls2*  |  Cardiolipin synthase                            | 
|  PG     |  population  |   33     |  stop gained c.703C>T p.Gln235\*                                                                |  PROKKA_01660|   *hemH*  |  Ferrochelatase                                  | 
|  PG     |  colony      | 33-s-c1  |  stop gained c.703C>T p.Gln235\*                                                                |  PROKKA_01660|   *hemH*  |  Ferrochelatase                                  | 
|  PG     |  population  |   33     |  g.2185374C>T                                                                                   |  intergenic  |    -      |  -                                               | 
|  PG     |  colony      | 33-s-c1  |  g.2185374C>T                                                                                   |  intergenic  |    -      |  -                                               | 
|  BAC    |  population  |   14     |  g.687271T>A                                                                                    |  intergenic  |    -      |  -                                               | 
|  BAC    |  colony      | 14-s-c2  |  g.687271T>A                                                                                    |  intergenic  |    -      |  -                                               | 
|  BAC    |  population  |    7     |  g.687271T>A                                                                                    |  intergenic  |    -      |  -                                               | 
|  BAC    |  colony      | 7-s-c1   |  g.687271T>A                                                                                    |  intergenic  |    -      |  -                                               | 
|  BAC    |  population  |   14     |  missense variant c.477C>G p.Asp159Glu                                                          |  PROKKA_02369|    -      |  Baeyer-Villiger flavin-containing monooxygenase | 
|  BAC    |  colony      | 14-s-c2  |  missense variant c.477C>G p.Asp159Glu                                                          |  PROKKA_02369|    -      |  Baeyer-Villiger flavin-containing monooxygenase | 
|  BAC    |  population  |    7     |  missense variant c.477C>G p.Asp159Glu                                                          |  PROKKA_02369|    -      |  Baeyer-Villiger flavin-containing monooxygenase | 
|  BAC    |  colony      | 7-s-c1   |  missense variant c.477C>G p.Asp159Glu                                                          |  PROKKA_02369|    -      |  Baeyer-Villiger flavin-containing monooxygenase | 
|  BAC    |  colony      | 14-s-c2  |  missense variant c.636G>T p.Gln212His                                                          |  PROKKA_01209|   *trpB*  |  Tryptophan synthase beta chain                  | 
|  BAC    |  colony      | 14-s-c2  |  frameshift variant c.618dupT p.Arg207fs                                                        |  PROKKA_01818|   *agrA*  |  Accessory gene regulator protein A              |
|  BAC    |  colony      | 14-s-c2  |  missense variant & inframe deletion c.622_630delCATAATATTinsTCTTTCp.His208_Ile210delinsSerPhe  |  PROKKA_01818|   *agrA*  |  Accessory gene regulator protein A              |
|  BAC    |  colony      | 14-s-c2  |  missense variant c.778G>A p.Ala260Thr                                                          |  PROKKA_01862|   *kdpD*  |  Sensor protein KdpD                             | 
|  BAC    |  colony      | 7-s-c1   |  synonymous variant c.591G>A p.Gln197Gln                                                        |  PROKKA_01039|    -      |  hypothetical protein                            | 
|  BAC    |  colony      | 7-s-c1   |  missense variant c.684A>T p.Glu228Asp                                                          |  PROKKA_01817|   *dpiB*  |  Sensor histidine kinase DpiB                    | 
|  BAC    |  colony      | 7-s-c1   |  synonymous variant c.525T>A p.Thr175Thr                                                        |  PROKKA_01929|   *czcD*  |  Cadmium, cobalt and zinc/H( )-K( ) antiporter   | 
|  BAC    |  population  |   35     |  missense variant c.277G>A p.Ala93Thr                                                           |  PROKKA_01081|   *topA*  |  DNA topoisomerase 1                             | 
|  BAC    |  colony      | 35-s-c2  |  missense variant c.277G>A p.Ala93Thr                                                           |  PROKKA_01081|   *topA*  |  DNA topoisomerase 1                             | 
|  BAC    |  population  |   35     |  synonymous variant c.804T>A p.Ile268Ile                                                        |  PROKKA_01457|    -      |  putative AAA domain-containing protein          | 
|  BAC    |  colony      | 35-s-c2  |  synonymous variant c.804T>A p.Ile268Ile                                                        |  PROKKA_01457|    -      |  putative AAA domain-containing protein          | 
|  BAC    |  population  |   35     |  stop gained c.2035C>T p.Gln679\*                                                               |  PROKKA_01463|   *relA*  |  GTP pyrophosphokinase                           | 
|  BAC    |  colony      | 35-s-c2  |  stop gained c.2035C>T p.Gln679\*                                                               |  PROKKA_01463|   *relA*  |  GTP pyrophosphokinase                           | 
|  CON    |  population  |   11     |  g.1960175G>A                                                                                   |  intergenic  |    -      |  -                                               | 

Table 2. Summary of all mutations. PG, pexiganan; BAC, benzalkonium chloride; CON, control.

# ATCC 6538

We will assemble the PacBio reads (when available) and polish with the illumina reads. But in the meantime we will work with an assembly of only the illumina reads.

## Assemble with `a5miseq` and annotate with `prokka`:

```sh
a5_pipeline.pl --threads 12 51_S1_L001_R1_001.fastq.gz 51_S1_L001_R2_001.fastq.gz atcc6538_a5
prokka --cpus 12 atcc6538_a5.contigs.fasta --genus Staphylococcus --prefix atcc6538_a5 --force &> atcc6538_a5.log
```

|File Name  | Contigs | Scaffolds | Genome Size | Longest Scaffold | N50  | Raw reads | EC Reads | % reads passing EC | Raw nt  | EC nt     | % nt passing EC | Raw cov | EC cov | Median cov | 10th percentile cov | bases >= Q40 | % GC |
|-----------|---------|-----------|-------------|------------------|------|-----------|----------|--------------------|---------|-----------|-----------------|---------|--------|------------|---------------------|--------------|------| 
|atcc6538_a5|   25    |  25       | 2775939     |  662992          |306853|  1090946  | 1080446  |  99.04             |308243493| 281315666 |     91.26       | 111.04  | 101.34 |  107       |   82                | 2775503      | 32.7 |

Table 3. Summary of ATCC6538 assembly.

The assembly ends up pretty good. I'm sure PacBio will get this is to 1 contig (plus any plasmids).

## Map and call variants

```sh
```

|Treatment|  Type        | Label    |  Mutation                                                                                       |  Locus tag   | Annotation|  Function                                        |
|---------|--------------|----------|-------------------------------------------------------------------------------------------------|--------------|-----------|--------------------------------------------------|
|  BAC    |  population  |    17    |  g.332893G>C                                                                                    |  Intergenic  |     -     |                      -                           | 
|  BAC    |  colony      | 17-s-c1  |  g.332893G>C                                                                                    |  Intergenic  |     -     |                      -                           |

Table 4. Summary of all mutations in ATCC6538 lines. BAC, benzalkonium chloride.

Can't find any variants in the control lines. Will check investigate the SNP a bit more ad also check for CNVs.

# CC398

## Assemble with `a5miseq` and annotate with `prokka`:

```sh
a5_pipeline.pl --threads 12 52_S2_L001_R1_001.fastq.gz 52_S2_L001_R2_001.fastq.gz cc398_a5
prokka --cpus 12 cc398_a5.contigs.fasta --genus Staphylococcus --prefix cc398_a5 --force &> cc398_a5.log
```

|File Name|Contigs|Scaffolds|Genome Size|Longest Scaffold|N50|Raw reads|EC Reads|% reads passing EC|Raw nt|EC nt|% nt passing EC|Raw cov|EC cov|Median cov|10th percentile cov|bases >= Q40|% GC|
|---------|-------|---------|-----------|----------------|---|---------|--------|------------------|------|-----|---------------|-------|------|----------|-------------------|------------|----|
|cc398    |   28  |    28   |  2781805  |     643841  |308974| 1036352 |1025876 |98.99       94411942| 266521751|      90.53    |105.83 | 95.81|   103    | 81                | 2781376    |32.8|

Table 5. Summary of CC398 assembly.

Also a pretty good assembly.

## Map and call variants

```sh
```

|Treatment|  Type        | Label    |  Mutation                                                                                       |  Locus tag   | Annotation|  Function                                        |
|---------|--------------|----------|-------------------------------------------------------------------------------------------------|--------------|-----------|--------------------------------------------------|
|  BAC    |  population  |    21    |  g.550439.G>A                                                                                   |  Intergenic  |     -     |                      -                           | 
|  BAC    |  colony      | 21-s-c1  |  g.550439.G>A                                                                                   |  Intergenic  |     -     |                      -                           |
|  BAC    |  population  |    22    |  g.550439.G>A                                                                                   |  Intergenic  |     -     |                      -                           | 
|  BAC    |  colony      | 22-s-c1  |  g.550439.G>A                                                                                   |  Intergenic  |     -     |                      -                           |
|  BAC    |  colony      | 22-s-c1  |  g.550439.G>A                                                                                   |  Intergenic  |     -     |                      -                           |
|  BAC    |  colony      | 21-s-c1  |  frameshift variant c.88delC p.Gln30fs                                                          |  PROKKA_00233|  *pbpX_1* |     Putative penicillin-binding protein PbpX     |
|  BAC    |  colony      | 22-s-c1  |  missense variant c.841C>A p.Gln281Lys                                                          |  PROKKA_00233|  *pbpX_1* |     Putative penicillin-binding protein PbpX     |
|  BAC    |  colony      | 40-s-c2  |  missense variant c.808G>A p.Glu270Lys                                                          |  PROKKA_00848|  *sigA*   |     RNA polymerase sigma factor SigA             |

Table 6. Summary of all mutations in CC398 lines. BAC, benzalkonium chloride.

Can't find any variants in the control lines. Will check investigate the SNP a bit more ad also check for CNVs.


