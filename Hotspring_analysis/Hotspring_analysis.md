---
editor_options: 
  markdown: 
    wrap: 72
---

# HOTSPRING ANALYSIS

There is an interesting effect in a hot spring located in Iceland:

Coinciding with the activity of a nearby volcano, the hot spring
undergoes events of very high temperature. After such episodes of high
temperature (close to 90 degrees), a bloom of algae living in the same
environment happens. It is proposed to perform an in depth genomic and
metagenomic exploration of this singular hot spring ecosystem, which
involved eight work packages:

1.  Metagenomic analysis
2.  Genome analysis (basic checks)
3.  Genome analysis (read mapping)
4.  Genome analysis (variant calling)
5.  Differential expression analysis
6.  Functional analysis of the responsible genes
7.  Phylogenetic analysis of the responsible genes
8.  Conclusions

## 1. Metagenomics

As a first step, it was run shotgun metagenomic sequencing of the
prokaryotic microbiome in two conditions obtained at different times: 1)
one sample taken during the high temperature episodes, and 2) another
sample taken right after the episodes, when the temperature is back to
normal and the bloom of algae has started. The DNA present in each
sample was extracted, and sent for Illumina shotgun sequencing. The raw
read files (reverse and forward) were produced by Illumina pair-end
sequencing:

-   Forward and reverse reads from the high temperature sample.
    metagenomics-hotspring-hightemp.1.fq.gz
    metagenomics-hotspring-hightemp.2.fq.gz

-   Forward and reverse reads from the normal temperature sample.
    metagenomics-hotspring-normaltemp.1.fq.gz
    metagenomics-hotspring-normaltemp.2.fq.gz

A taxonomic profiling of the samples was performed using the tool mOTUs:

``` bash

cd ~/Hotspring_analysis/metagenomics/

motus profile -f ~/Hotspring_analysis/sequencing_files/metagenomics-hotspring-hightemp.1.fq.gz -r ~/Hotspring_analysis/sequencing_files/metagenomics-hotspring-hightemp.2.fq.gz -o HighT.motus -t 20

motus profile -f ~/Hotspring_analysis/sequencing_files/metagenomics-hotspring-normaltemp.1.fq.gz -r ~/Hotspring_analysis/sequencing_files/metagenomics-hotspring-normaltemp.2.fq.gz -o NormalT.motus -t 20

awk '{print $NF,$0}' HighT.motus | sort -n | cut -f2- -d' ' | grep -v '#' | tail -1
```

The most abundant organism in high-temperature is *Aquifex aeolicus*
with a relative abundance of 0.9462243550 In NCBI we can find a
reference genome and the article of origin: Deckert G et al. The
complete genome of the hyperthermophilic bacterium Aquifex aeolicus.
Nature. 1998 Mar 26;392(6674):353-8. doi: 10.1038/32831. *A.aeolicus* is
a Gram negative bacteria with the shape of cocci. It is aerobic,
chemolithoautotrophic and hyperthermophilic (optimum temperature 96° C).
It has a genome of 1.59079 Mb that encodes 1724 proteins.

`awk '{print $NF,$0}' NormalT.motus | sort -n | cut -f2- -d' ' | grep -v '#' | tail -1`

In normal temperature, the most abundant organism is also *Aquifex
aeolicus* with 0.0432142543 relative abundance.

`perl -F"\t" -lane 'print if $F[1]>0' HighT.motus | wc -l`

`perl -F"\t" -lane 'print if $F[1]>0' NormalT.motus | wc -l`

Nevertheless, The alpha biodiversity is greater in normal temperature,
with 230 detected organisms. In high temperature there are only 16
organisms detected.

## 2. Genome Analysis (basic checks)

To further characterize *A.aeolicus*, it is performed a RNA-Seq analysis
of samples from both cultures at normal-temperature and high-temperature
conditions, two biological replicates each. The genome is already
assembled and a quality checking for the reads was performed in all
samples, providing only high quality reads in fasta format. There are 4
samples, two replicates of each temperature in forward and reverse
orientation:

-   hightemp01.r1.fq (sample 1, forward)
-   hightemp01.r2.fq (sample 1, reverse)
-   hightemp02.r1.fq (sample 2, forward)
-   hightemp02.r2.fq (sample 2, reverse)
-   normal01.r1.fq (sample 1, forward)
-   normal01.r2.fq (sample 1, reverse)
-   normal02.r1.fq (sample 2, forward)
-   normal02.r2.fq (sample 2, reverse)

``` bash

cd ~/Hotspring_analysis/sequencing_files/RNAseq/

for f in *.gz ; do gunzip -c "$f" > ~/Hotspring_analysis/genome_analysis/"${f%.*}" ; done

cd ~/Hotspring_analysis/genome_analysis/

cat normal02.r1.fq | grep -E '^[ATCG]' | awk '{ print length }' | uniq -c

cat hightemp02.r1.fq | grep -E '^[ATCG]' | awk '{ print length }' | uniq -c

cat hightemp01.r2.fq | grep -vE '^[ATCG]' | grep -v '@' | grep -v '+' | uniq -c

cat normal01.r2.fq | grep -vE '^[ATCG]' | grep -v '@' | grep -v '+' | uniq -c
```

There are 318730 reads in each of the high-temperature files and 288742
in each of the normal-temperature files. They are paired-end reads and
have a length of 100 bases, so they were produced by Illumina paired-end
sequencing. All the reads in all the files have the maximum quality of
Illumina, 40, represented as the letter I.

## 3. Genome Analysis (read mapping)

Several downstream analyses were performed, including variant calling,
expression analysis, etc. First the reads were mapped to the assembled
genome.

``` bash

bwa index ~/Hotspring_analysis/genome_analysis/genome.fasta

bwa mem genome.fasta hightemp01.r1.fq hightemp01.r2.fq > hightemp01.sam

bwa mem genome.fasta hightemp02.r1.fq hightemp02.r2.fq > hightemp02.sam

bwa mem genome.fasta normal01.r1.fq normal01.r2.fq > normal01.sam

bwa mem genome.fasta normal02.r1.fq normal02.r2.fq > normal02.sam

conda activate samtools

samtools view -b -h hightemp01.sam > hightemp01.bam

samtools view -b -h hightemp02.sam > hightemp02.bam

samtools view -b -h normal01.sam > normal01.bam

samtools view -b -h normal02.sam > normal02.bam

rm hightemp01.sam hightemp02.sam normal01.sam normal02.sam

samtools sort hightemp01.bam > hightemp01.sorted.bam

samtools sort hightemp02.bam > hightemp02.sorted.bam

samtools sort normal01.bam > normal01.sorted.bam

samtools sort normal02.bam > normal02.sorted.bam

samtools flagstat hightemp01.sorted.bam

samtools flagstat normal01.sorted.bam
```

There are 637995 records in the first duplicate of high-temperature
sample, and 638010 records in the second duplicate. On the other hand,
there are 577504 records in the first duplicate of normal-temperature
sample and 577497 records in the second duplicate. In the
high-temperature samples, there are 318730 different reads for forward
and reverse; and there are 288742 different reads in both
normal-temperature samples. If we sum the number of reads in forward and
reverse orientation, we almost get the number of records. The reads that
are missing are supplementary reads and they don't align.

## 4. Genome Analysis (variant calling)

Using the mappings, a variant calling analysis was carried out. Maybe
some mutation is related to the sudden proliferation of these organisms.

``` bash

samtools merge -h normal01.sorted.bam merged.bam normal01.sorted.bam normal02.sorted.bam hightemp01.sorted.bam hightemp02.sorted.bam

samtools faidx genome.fasta

conda activate bcftools

bcftools mpileup -f genome.fasta merged.bam > merged.vcf

bcftools call --threads 1 -mv -Ob -o calls.bcf merged.vcf

bcftools mpileup -f genome.fasta normaltemp01.sorted.bam normaltemp02.sorted.bam hightemp01.sorted.bam hightemp02.sorted.bam > separate.vcf

bcftools call -mv -Ob -o variant_calling_separate.vcf separate.vcf

bcftools view calls.bcf | grep -v "^##" | head

bcftools stats calls.bcf | more

bcftools view -H calls.bcf | sort -k6,6gr
```

3 variants were obtained, all of them are SNPs (2 T\>A and 1 T\>G). Only
one variant has quality and depth of coverage greater than 100. With the
samples merged, the SNP T\>G has 213 of quality and 249 in coverage; and
with the samples separated, T\>G has a quality of 943 and a coverage of
682.

``` cp ~/Hotspring_analysis/sequencing_files/genome.gff ~/``Hotspring_analysis``/genome_analysis/genome.gff ```

`cat genome.gff | awk -F "\t" '{ if ((1264965 <= $5 ) && (1264965 >= $4 )) print $0 }'`

The variant with the best quality, the SNP T\>G is in position 1264965
of the *Aquifex aeolicus* genome, so the gene that receives the variant
is nifA. Using Integrative Genomics Viewer (IGV), it can be shown the
mutation site (Figure 1).

![**Figure 1.** nifA mutation. IGV
visualization](Images/igv_snapshot.png "nifA mutation")

## 5. Differential expression analysis

Given the expression data for each sample, it can be compared the
expression differences between the samples grown under normal and
high-temperature conditions, which could provide additional information
about important genes involved.

First of all, it is done a ?read-count? using `htseq-count` and the
following files were generated:

-   hightemp01.count
-   hightemp02.count
-   normal01.count
-   normal01.count

``` bash

conda activate samtools; for file in *.sorted.bam; do samtools index $file; done; conda deactivate

for file in *.sorted.bam; do htseq-count -i locus_tag -t CDS $file genome.gff > "${HOME}/Hotspring_analysis/DEA/${file%.sorted.bam}.counts"; done
```

The Differential expression analysis process is described in the
corresponding Rmarkdown file, DEA.Rmd.

Five genes showed a statistical (p-adj \< 0.01) differential expression.
They were annotated using the genome.gff file and table 1 was generated
showing the altered genes, their description, name, p-val, p-adj and
fold change:

![**Table 1.** Differential
expression](Images/tabla.png "Differential expression")

AQUIFEX_01423, AQUIFEX_01759 and AQUIFEX_01761 are over expressed under
high-heat conditions, while AQUIFEX_01754 and AQUIFEX_01760 are under
expressed.

## 6. Functional analysis

The differential expression results gave an idea about potentially
important overexpressed genes.

`cd ~/Hotspring_analysis/functional_analysis/`

`awk 'BEGIN {RS=">"} /AQUIFEX_01423/ {print ">"$0}' proteome.faa > AQUIFEX_01423.fasta`

`awk 'BEGIN {RS=">"} /AQUIFEX_01759/ {print ">"$0}' proteome.faa > AQUIFEX_01759.fasta`

`awk 'BEGIN {RS=">"} /AQUIFEX_01761/ {print ">"$0}' proteome.faa > AQUIFEX_01761.fasta`

First, I performed a BLASTp with the sequences of the 3 proteins
overexpressed. Then I search for the function in UniprotKB. For
AQUIFEX_01423, there?s a 100% match with sigma-54-dependent Fis family
transcriptional regulator of *Aquifex aeolicus* and *Aquifex sp*. When
its function is searched in Uniprot (accession WP_010881164) it is
obtained that it is in charge of ATP binding, sequence-specific DNA
binding and regulation of DNA-templated transcription. Kegg indicates
that it is a transcriptional regulator of the NifA subfamily, also
related to ATPase AAA+. This protein has three domains: GAF (named after
some of the proteins it is found in, including cGMP-specific
phosphodiesterases, adenylyl cyclases and FhlA); Sigma-54 interaction
domain; Bacterial regulatory protein, Fis family HTH_8 domain. For
AQUIFEX_01759, BLASTp gives 98.67% confidence to FeMo cofactor
biosynthesis protein NifB, and 100% to radical SAM protein from
*Methanococcus maripaludis*. Checking these proteins in Uniprot and
kegg, we see that, indeed, the two results state that it is a cofactor
of the maturase NifB related to the apical meristem assembler protein.
The functions are catalytic activity, iron-sulfur cluster binding and
metal ion binding. For AQUIFEX_01759 there is a single domain, radical
SAM domain. For AQUIFEX_01761, we have a 100% coverage and identity with
nitrogenase iron protein from *Methanococcus maripaludis*. In Uniprot
the function that appears suggests that it is part of the nitrogenase
complex which has 2 components: the iron protein and the molybdenum-iron
protein and it catalyzes the nitrogen fixation. Also there are some
other functional annotations: 4 iron, 4 sulfur cluster binding, ATP
binding, carbonyl sulfide nitrogenase activity, metal ion binding and
nitrogenase activity. There is also a single domain in this protein,
Fer4_NifH (4Fe-4S iron sulfur cluster binding proteins, NifH/frxC
family).

The genes that are overexpressed are Nif genes. Their proteins are
responsible for the fixation of atmospheric nitrogen to ammonium under
conditions of low dissolved N~2~. So, when there is a lack of N~2~,
AQUIFEX_01423 is responsible for activating the fixation, which involves
the activation of the nitrogenase complex. In STRING database, it says
that: The Nitrogen fixation protein nifB is involved in the biosynthesis
of the iron-molybdenum cofactor found in the dinitrogenase enzyme of the
nitrogenase complex in nitrogen-fixing microorganisms. And we know that
AQUIFEX_01759 is the cofactor of NifB. That nitrogenase complex has
AQUIFEX_01761 as a component so it is to be expected that its production
will also increase.

As an hypothesis, when the hot spring reaches high temperatures,
*Aquifex aeolicus* start to grow in its ideal temperature and due to the
high temperatures, the gases dissolved in water become less soluble
generating conditions of low N~2~ concentration which makes the genes
related to nitrogen fixation to activate. When the hot spring reaches
normal temperatures, the population of *A. aeolicus* descends and the
nitrogen levels are high enough to create an ideal condition for algae
to grow and create a bloom.

## 7. Phylogenetic analysis

The functional analysis of the overexpressed genes gave an idea about
the biological processes happening in the hotspring during the high
temperature episodes. To refine the functional inferences and
investigate the evolutionary origin of the overexpressed genes in the
isolated genome, it is performed an in depth phylogenetic analysis
comparing each overexpressed gene against their homologs in other
prokaryotic genomes. The reference set of organisms includes the public
genome of the same isolated species , and 6 other bacteria and archaea
that are known to be related to the biological processes identified
previously:

-   CLOPA - *Clostridium pasteurianum*
-   9AQUI - *Hydrogenivirga caldilitoris*
-   METV3 - *Methanococcus voltae*
-   NOSS1 - *Nostoc sp.*
-   AQUAE - *Aquifex aeolicus* (strain VF5)
-   METMP - *Methanococcus maripaludis*
-   RHOCB - *Rhodobacter capsulatus*

The complete proteome of all 7 species are downloaded from Uniprot and
put them together into a FASTA file (all_reference_proteomes.faa).

For each over expressed gene, it is performed a standard phylogenetic
workflow:

1.  Run a blast search for each over expressed protein against all
    reference proteomes
2.  Extract hits with e-value \<= 0.001
3.  Create a FASTA file with all the sequences of selected hits
4.  Build a phylogenetic tree out of the FASTA file
5.  Visualize the result

``` bash

cd ~/Hotspring_analysis/phylogenetic_analysis/

makeblastdb -dbtype prot -in all_reference_proteomes.faa -out proteomes.blastdb

blastp -task blastp -query AQUIFEX_01423.fasta -db proteomes.blastdb -outfmt 6 -evalue 0.001 > A_01423_homologs.blastout

blastp -task blastp -query AQUIFEX_01759.fasta -db proteomes.blastdb -outfmt 6 -evalue 0.001 > A_01759_homologs.blastout

blastp -task blastp -query AQUIFEX_01761.fasta -db proteomes.blastdb -outfmt 6 -evalue 0.001 > A_01761_homologs.blastout

python /home/Hotspring_analysis/phylogenetic_analysis/scripts/extract_seqs_from_blast_result.py A_01423_homologs.blastout all_reference_proteomes.faa > A_01423_homologs.faa

python /home/Hotspring_analysis/phylogenetic_analysis/scripts/extract_seqs_from_blast_result.py A_01759_homologs.blastout all_reference_proteomes.faa > A_01759_homologs.faa

python /home/Hotspring_analysis/phylogenetic_analysis/scripts/extract_seqs_from_blast_result.py A_01761_homologs.blastout all_reference_proteomes.faa > A_01761_homologs.faa

cat AQUIFEX_01423.fasta >> A_01423_homologs.faa

cat AQUIFEX_01759.fasta >> A_01759_homologs.faa

cat AQUIFEX_01761.fasta >> A_01761_homologs.faa

mafft A_01423_homologs.faa > A_01423_homologs.alg

mafft A_01759_homologs.faa > A_01759_homologs.alg

mafft A_01761_homologs.faa > A_01761_homologs.alg

iqtree -s A_01423_homologs.alg -m LG

iqtree -s A_01759_homologs.alg -m LG

iqtree -s A_01761_homologs.alg -m LG

python /home/Hotspring_analysis/phylogenetic_analysis/scripts/midpoint_rooting.py A_01423_homologs.alg.treefile | ete3 view --text

python /home/Hotspring_analysis/phylogenetic_analysis/scripts/midpoint_rooting.py A_01759_homologs.alg.treefile | ete3 view --text

python /home/Hotspring_analysis/phylogenetic_analysis/scripts/midpoint_rooting.py A_01761_homologs.alg.treefile | ete3 view --text
```

![**Figure
2.**](Images/AQUIFEX_01423.PNG "AQUIFEX_01423 phylogenetic tree")

![**Figure
3.**](Images/AQUIFEX_01759.PNG "AQUIFEX_01759 phylogenetic tree")

![**Figure
4.**](Images/AQUIFEX_01761.PNG "AQUIFEX_01761 phylogenetic tree")

The closes ortholog for AQUIFEX_01423 is
tr\|A0A497XW95\|A0A497XW95_9AQUI from *Hydrogenivirga caldilitoris*
(Figure 2); the closest one to AQUIFEX_1759 is tr\|Q6LZH0\|Q6LZH0_METMP
from *Methanococcus maripaludis* (Figure 3); and the closest ortholog
for AQUIFEX_1761 is sp\|P0CW57\|NIFH_METMP from *Methanococcus
maripaludis* (Figure 4).

The ortholog of AQUIFEX_01423 has the same functions (ATP binding and
DNA regulation) and domains (GAF_2, HTH_8 and Sigma54_activat); the
ortholog of AQUIFEX_01759 has similar functions (catalytic activity,
iron-sulfur cluster binding and metal ion binding), participate in the
nitrogen fixation process and has the same domain (Radical_SAM); and the
ortholog of AQUIFEX_01761 is a nitrogenase iron protein, part from the
nitrogenase complex as the A. aeolicus protein with the same function
and domain (Fer4_NifH); so these orthology assignments support the
previous functional annotations.

Since only AQUIFEX_01423 produced a hit in the public genome of A.
aeolicus, neither AQUIFEX_01759 nor AQUIFEX_01761 are present on the
public genome.

AQUIFEX_01423 is native from *A. aeolicus* and its origin is possibly
speciation as it is related to Nif gene from *Hydrogenivirga
caldilitoris* present in the same branch. For both AQUIFEX_01759 and
AQUIFEX_01761, it seems probable that they appeared via vertical or
horizontal transfer from *Methanococcus maripaludis* as it was one of
the most abundant organisms in the normal-temperature condition.

## 8. Conclusions

The hotspring in normal temperature has a high alpha biodiversity, but
when the temperature reaches about 90 degrees, only the hyperthermophile
organisms like *Aquifex aeolicus* and perhaps, other organisms in
sporulated form can survive. In normal conditions, the gases like O~2~
and N loose solubility in water due to the high temperatures and the
nutrients in the hotspring descend. Thanks to lateral transfer, *Aquifex
aeolicus* can obtain two genes from the second bigger population
(*Methanococcus maripaludis*) that can help it to obtain more nutrients
by nitrogen fixation in its optimal temperature. This event causes the
bacteria to overgrow in high temperature and fix a lot of nitrogen. When
the temperatures descend and the alpha biodiversity restores the
previous levels, *Aquifex aeolicus* population descends because of
temperature and all the nitrogen fixed causes algae to overgrow and
create blooms.
