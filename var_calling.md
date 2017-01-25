# Table of contents

- [SH1000](#sh1000)
 * [Check the reference](#check-the-reference)
 * [Map and call variants](#map-and-call-variants)
 * [SH1000 variants](#sh1000-variants)
 * [SH1000 CNVs](#sh1000-cnvs)
- [ATCC 6538](#atcc-6538)
 * [Assemble ATCC6538 PacBio reads with `canu`, correct with `pilon`, and annotate with `prokka`](#assemble-atcc6538-pacbio-reads-with-canu-correct-with-pilon-and-annotate-with-prokka)
 * [Map and call ATCC6538 variants](#map-and-call-atcc6538-variants)
 * [ATCC 6538 CNVs](#atcc-6538-cnvs)
- [CC398](#cc398)
 * [Assemble CC398 PacBio reads with `canu`, correct with `pilon`, and annotate with `prokka`](#assemble-cc398-pacbio-reads-with-canu-correct-with-pilon-and-annotate-with-prokka)
 * [Map and call CC398 variants](#map-and-call-cc398-variants)
 * [CC398 CNVs](#cc398-cnvs)

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

## SH1000 CNVs

```sh
for i in `cut -f5 ../sh1000.tsv | grep -v suffix`
do
  bwa mem -t 12 -R "@RG\tID:${i}\tSM:${i}" ../sh1000_dueppel/sh1000_dueppel.fna ../../../${i}_L001_R1_001.fastq.
gz ../../../${i}_L001_R2_001.fastq.gz > $i.sam
  picard-tools SortSam INPUT=$i.sam OUTPUT=$i.bam SORT_ORDER=coordinate
  picard-tools BuildBamIndex INPUT=$i.bam
done
```

```r
library("cn.mops")
library("magrittr")
BAMFiles <- list.files(pattern=".bam$")
bamDataRanges <- getReadCountsFromBAM(BAMFiles, mode="paired")
res <- haplocn.mops(bamDataRanges)
res <- calcIntegerCopyNumbers(res)
plot(res, which=1)
```

No CNVs. Check coverage plots:

```sh
for i in *bam; do /home/paul/opt/IGVTools/igvtools count $i $i.tdf ../sh1000_dueppel/sh1000_dueppel.genome; done
```

![](https://github.com/Perugolate/biocides/blob/master/figs/sh1000_count.png)

Nothing obvious here. The high-coverage region is a transposon that we always see in SH1000.

# ATCC 6538

## Assemble ATCC6538 PacBio reads with `canu`, correct with `pilon`, and annotate with `prokka`

```sh
# assemble the pacbio reads
canu -p pb1 -d pb1 genomeSize=2.7m -pacbio-raw 11829_1.filtered_subreads.fastq.gz
cd pb1
# this is annotation is not necessary - just for seeing where the errors were
prokka --cpus 12 pb2.contigs.fasta --genus Staphylococcus --prefix pb1 &> pb1.log
cd pb1
# align cognate illumina library to correct the errors
bwa index pb1.fna
bwa mem -t 12 pb1.fna  ../../../../51_S1_L001_R1_001.fastq.gz ../../../../51_S1_L001_R2_001.fastq.gz  > 51.sam
picard-tools SortSam INPUT=51.sam OUTPUT=51.bam SORT_ORDER=coordinate
picard-tools BuildBamIndex INPUT=51.bam
# use the alignment to correct the assembly
java -jar ~/opt/pilon-1.21.jar --genome pb1.fna --frags 51.bam
# annotate the corrected assembly
prokka --cpus 12 pilon.fasta --genus Staphylococcus --prefix pb1.pi &> pb1.pi.log
```

This produces a complete assembly.

## Map and call ATCC6538 variants

```sh
cd pb1.pi
for i in `cat atcc6538.tsv`
do
  snippy --cpus 12 --outdir $i --reference pb1.pi.gbk --rgid $i --prefix $i --pe1 ../../../../../${i}_L001_R1_001.fastq.gz --pe2 ../../../../../${i}_L001_R2_001.fastq.gz &> $i.log
done
```

|Treatment|  Type        | Label    |  Mutation                                                                                       |  Locus tag   | Annotation|  Function                                        |
|---------|--------------|----------|-------------------------------------------------------------------------------------------------|--------------|-----------|--------------------------------------------------|
|  BAC    |  population  |    17    |  g.1724366C>G                                                                                    |  Intergenic  |     -     |                      -                           | 
|  BAC    |  colony      | 17-s-c1  |  g.1724366C>G                                                                                    |  Intergenic  |     -     |                      -                           |

Table 3. Summary of all mutations in ATCC6538 lines. BAC, benzalkonium chloride.

Can't find any variants in the control lines.

## ATCC 6538 CNVs

```sh
bwa index pb2.pi.fna
for i in `cat atcc6538.tsv`
do
  bwa mem -t 12 -R "@RG\tID:${i}\tSM:${i}" pb1.pi.fna ../../../${i}_L001_R1_001.fastq.gz ../
../../${i}_L001_R2_001.fastq.gz > $i.sam
  picard-tools SortSam INPUT=$i.sam OUTPUT=$i.bam SORT_ORDER=coordinate
  picard-tools BuildBamIndex INPUT=$i.bam
done
```

```r
library("cn.mops")
library("magrittr")
BAMFiles <- list.files(pattern=".bam$")
bamDataRanges <- getReadCountsFromBAM(BAMFiles, mode="paired")
res <- haplocn.mops(bamDataRanges)
res <- calcIntegerCopyNumbers(res)
plot(res, which=1)
```

No CNVs.

```sh
for i in *bam; do /home/paul/opt/IGVTools/igvtools count $i $i.tdf pb1.pi.fna; done
```

![](https://github.com/Perugolate/biocides/blob/master/figs/atcc6538_count.png)

# CC398

## Assemble CC398 PacBio reads with `canu`, correct with `pilon`, and annotate with `prokka`

```sh

canu -p pb2 -d pb2 genomeSize=2.7m -pacbio-raw 11829_2.filtered_subreads.fastq.gz
cd pb2
prokka --cpus 12 pb2.contigs.fasta --genus Staphylococcus --prefix pb2 &> pb2.log
cd pb2
bwa index pb2.fna
bwa mem -t 12 pb2.fna  ../../../../52_S2_L001_R1_001.fastq.gz ../../../../52_S2_L001_R2_001.fastq.gz  > 52.sam
picard-tools SortSam INPUT=52.sam OUTPUT=52.bam SORT_ORDER=coordinate
picard-tools BuildBamIndex INPUT=52.bam
java -jar ~/opt/pilon-1.21.jar --genome pb2.fna --frags 52.bam
prokka --cpus 12 pilon.fasta --genus Staphylococcus --prefix pb2.pi &> pb2.pi.log
```

As with ATCC6538, this results in what appears to be a complete assembly.

## Map and call CC398 variants

```sh
cd pb2.pi
for i in `cat /media/tera5/pigs/biocides/cc398/cc398.tsv`
do
  snippy --cpus 12 --outdir $i --reference pb2.pi.gbk --rgid $i --prefix $i --pe1 ../../../../../${i}_L001_R1_001.fastq.gz --pe2 ../../../../../${i}_L001_R2_001.fastq.gz &> $i.log
done
```

|Treatment|  Type        | Label    |  Mutation                                                                                       |  Locus tag   | Annotation|  Function                                        |
|---------|--------------|----------|-------------------------------------------------------------------------------------------------|--------------|-----------|--------------------------------------------------|
|  BAC    |  population  |    21    |  g.1051564.C>T                                                                                   |  Intergenic  |     -     |                      -                           | 
|  BAC    |  colony      | 21-s-c1  |  g.1051564.C>T                                                                                   |  Intergenic  |     -     |                      -                           |
|  BAC    |  population  |    22    |  g.1051564.C>T                                                                                   |  Intergenic  |     -     |                      -                           | 
|  BAC    |  colony      | 22-s-c1  |  g.1051564.C>T                                                                                   |  Intergenic  |     -     |                      -                           |
|  BAC    |  colony      | 21-s-c1  |  frameshift variant c.88delC p.Gln30fs                                                          |  PROKKA_01239|  *pbpX* |     Putative penicillin-binding protein PbpX     |
|  BAC    |  colony      | 22-s-c1  |  missense variant c.841C>A p.Gln281Lys                                                          |  PROKKA_01239|  *pbpX* |     Putative penicillin-binding protein PbpX     |
|  BAC    |  colony      | 40-s-c2  |  missense variant c.808G>A p.Glu270Lys                                                          |  PROKKA_01802|  *sigA*   |     RNA polymerase sigma factor SigA             |

Table 4. Summary of all mutations in CC398 lines. BAC, benzalkonium chloride.

The population sample 40 and the cognate colony sample 40-s-c2 are a good example of why it can be useful to have both types of sample: The mutation g.1955940.C>T (which cause the substitution c.808G>A p.Glu270Lys) is called in 40-s-c2 but not in 40. However when you inspect the alignments (see the figure below) you can see that the mutation is present in the population but not at a high enough frequency for it to be identified by the variant caller.

![](https://github.com/Perugolate/biocides/blob/master/figs/cc398_snp1955940.png)

The same is true for frameshift variant c.88delC p.Gln30fs in samples 21 and 21-s-c1:

![](https://github.com/Perugolate/biocides/blob/master/figs/cc398_1350970.png)

And also for missense variant c.808G>A p.Glu270Lys in samples 22 and 22-s-c1:

![](https://github.com/Perugolate/biocides/blob/master/figs/cc398_1351724.png)

## CC398 CNVs

```sh
bwa index pb2.pi.fna
for i in `cat cc398.tsv`   
do
  bwa mem -t 12 -R "@RG\tID:${i}\tSM:${i}" pb2.pi.fna ../../../../../${i}_L001_R1_001.fastq.gz ../../../../../${i}_L001_R2_001.fastq.gz > $i.sam       
  picard-tools SortSam INPUT=$i.sam OUTPUT=$i.bam SORT_ORDER=coordinate
  picard-tools BuildBamIndex INPUT=$i.bam
done
```

```r
library("cn.mops")
library("magrittr")
BAMFiles <- list.files(pattern=".bam$")
bamDataRanges <- getReadCountsFromBAM(BAMFiles, mode="paired")
res <- haplocn.mops(bamDataRanges)
res <- calcIntegerCopyNumbers(res)
plot(res, which=1)
```

No CNVs.

```sh
for i in *bam; do /home/paul/opt/IGVTools/igvtools count $i $i.tdf pb2.pi.fna; done
```

![](https://github.com/Perugolate/biocides/blob/master/figs/cc398_count.png)

