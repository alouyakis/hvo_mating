project name: Hvo mating expression  
author: Artemis Louyakis  
date: 2017.12.01  

---
seq info: Illumina MiSeq 2x75  
samples: 24
treatments: time (0,2,4,8,24 hours) vs no mating control (shaking)
date of sequencing: 2014.03.29 [note: fastqs regenerated on 2017.12.01 due to missing sample entry in seq sample sheet]
seqs owner: Andrea Makkay
---
programs used:
  Scythe
  Sickle
  #bowtie2 v2.3.1
  Salmon v0.14.0 / v0.9.1
  samtools v1.3.1
  cufflinks v2.2.1 (http://cole-trapnell-lab.github.io/cufflinks/)
---
organism: Haloferx volcanii DS2
ref_seq: GCF_000025685.1_ASM2568v1
  length: 4.0129 Mb
  %GC: 65.46
  gene count: 4023
---
```shell
wget -r ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/Haloferax_volcanii/latest_assembly_versions/GCF_000025685.1_ASM2568v1/*
```

#### quality control - trim reads
```shell
# run in raw_seqs directory
for fq in */raw_seqs/*.gz; do gzip -d $fq & done

module load scythe/0.991
module load sickle/1.33

mkdir -p ../trimming

output="../trimming"
for fwdseq in *_R1*; do
  revseq="${fwdseq/_R1/_R2}"
  scythe -q sanger -a ../../adapters.fasta \
      -o ${output}/${fwdseq}.adapt ${fwdseq};
  scythe -q sanger -a ../../adapters.fasta \
      -o ${output}/${revseq}.adapt ${revseq};
  sickle pe -f ${output}/${fwdseq}.adapt -r ${output}/${revseq}.adapt \
      -t sanger -o ${output}/${fwdseq} -p ${output}/${revseq} \
      -s /dev/null;
  mv ${output}/${fwdseq} ${output}/${fwdseq/fastq/fq};
  mv ${output}/${revseq} ${output}/${revseq/fastq/fq};
done
```

#### Bowtie2 alignment to reference genome
```shell
#run in directory containing downloaded genome directory
for genome in GCF*; do gzip -d $genome/*.gz & done

module load bowtie2/2.3.1
module load samtools/1.3.1

mkdir -p indexes
for genome in GCF*; do \
  bowtie2-build --threads 20 \
    ${genome}/${genome}_genomic.fna \
    indexes/${genome}

#run in directory containing transcript directories and genome subdirectory
#experiment="/home/CAM/alouyakis/transcriptomes/exp_*"
#mkdir -p $experiment"/coverage"

index="hfx_volc_ref/indexes/GCF_000025685.1_ASM2568v1"

for exp in exp*; do \
  for seq in ${exp}/trimming/*R1*fq; do \
    seqbase=${seq/_R1_001.fq/_}
    bowtie2 -q --threads 20 -x ${index} \
      -1 ${seqbase}R1_001.fq -2 ${seqbase}R2_001.fq | samtools view -b - | \
      samtools sort -@ 20 -m 6G -o ${exp}/coverage/$(basename ${seqbase/_L001_/.bam});
  done;
done
```

#### housekeeping
```shell 
# moved files from xanadu to imac pro

cd /Users/gogartenlab/Desktop/transcriptomes/exp_mating/trimming
for f in *.fq; do mv "$f" "$f.tmp"; mv "$f.tmp" "`echo $f | tr "[:upper:]" "[:lower:]"`"; done
for f in *.fq; do mv "$f" "${f/_001.fq/.fq}"; done
for f in *.fq; do mv "$f" "${f/_l001_r/_r}"; done
for f in *.fq; do mv "$f" "${f/_s[0-9]_/_}"; done
for f in *.fq; do mv "$f" "${f/_s[0-9][0-9]_/_}"; done
for f in *.fq; do mv "$f" "${f/h_/_}"; done
```

#### Salmon alignment and quantification to reference genome
```shell
salmon index \
  -t Haloferax_volcanii/GCF_000025685.1_ASM2568v1/GCF_000025685.1_ASM2568v1_cds_from_genomic.fna.gz \
  -i hfx_volc_ref/salmon_index/hfx_02_genomic_fna

salmon_index="hfx_volc_ref/salmon_index/hfx_02_genomic_fna"
data="exp_mating/trimming/"
for fn in exp_mating/trimming/*r1.fq; do
  base=$(basename -s _r1.fq "$fn")
  samp=${base/_r1.fq/_}
  echo "Processing sample $samp"
  salmon quant -i ${salmon_index} -l A \
    -1 ${data}/${samp}_r1.fq \
    -2 ${data}/${samp}_r2.fq \
    -p 8 -o quants/mating/${samp}_quant
done

## run in quants/mating/
salmon quantmerge --quants ap0_quant a0_quant a2_quant a4_quant a8_quant a24_quant bp0_quant b0_quant b2_quant b4_quant b8_quant b24_quant -o mating.mat
```

#### differential expression
Run differential expression analysis in R - see hvo_mating.Rmd

```shell
###### make maps for annotating output ######
#move into directory containing genome directories

## annotation map prepared using GCF_000025685.1_ASM2568v1_feature_table.txt [find script]
grep ">" GCF_000025685.1_ASM2568v1_cds_from_genomic.fna > \
  ../../annot_hvo/extra/GCF_000025685.1_ASM2568v1_cds_from_genomic.fna.headers
sed 's/ \[/\t\[/g' GCF_000025685.1_ASM2568v1_cds_from_genomic.fna.headers > \
  GCF_000025685.1_ASM2568v1_cds_from_genomic.fna.headers.tabs
awk -F "\t" '{OFS="\t"} {print $1,$2,$4}' GCF_000025685.1_ASM2568v1_cds_from_genomic.fna.headers.tabs > \
  GCF_000025685.1_ASM2568v1_cds_from_genomic.fna.headers.tabs.col1-2-4
sed 's/\[locus_tag=//g; s/\[protein=//g; s/\]//g; s/^>//g' GCF_000025685.1_ASM2568v1_cds_from_genomic.fna.headers.tabs.col1-2-4 > \
  GCF_000025685.1_ASM2568v1_cds_from_genomic.fna.headers.tabs.col1-2-4.clean
## not all rows print out correctly, but few enough errors to correct by hand. 
grep "gene=" GCF_000025685.1_ASM2568v1_cds_from_genomic.fna.headers
## and correct those that were out of column alignment; will find an automated method later to avoid this
head annomap.txt
## how did i make hfx_volc_kegg.txt?
sed 's/^hvo\://' hfx_volc_kegg.txt > hfx_volc_kegg.tmp

### filter the gff file for genome attributes separated by tabs and prepare map of old and new locus tags
awk -F "\t" '{print $9}' hfx_volc_ref/GCF_000025685.1_ASM2568v1/GCF_000025685.1_ASM2568v1_genomic.gff |\
   sed 's/;/\t/g' | grep "ID=gene" | grep "old_locus_tag=" |   awk -F "\t" '{OFS="\t"} {print $6,$7}' |\
    sed 's/old_locus_tag=//g; s/locus_tag=//g'   > locus_tag_map.txt.tmp
    ###hand correct and remove .tmp

### use kegg website to get hvo ids and corresponding kegg ids; save in a text file
### http://www.genome.jp/kegg-bin/get_htext?hvo00001 [click fourth arrow to open all and save page]
### click download htext [hvo00001.keg]
## make map
sed 's/<b>//g; s/<\/b>//g; s/^D      /\t\t\t/; s/^C    /\t\t/; s/^B  /\t/; s/^B//; s/^A//; s/^\#//; /^$/d; /^\s*$/d; s/"//g' hvo00001.keg > hvo_kegg.txt
awk -F "\t" '{OFS="\t"} {
for (i=1;i<=NF;++i) if ($i != "") a[i] = $i;
if (na < NF) na = NF;
for (i=1;i<na;++i) printf "%s\t", a[i]
printf "%s\n", a[na];
}' hvo_kegg.txt > hvo_kegg_fill.txt

## calculate gene counts
awk -F "\t" '{OFS="\t"}{print $4,$5,$9}' hfx_volc_ref/GCF_000025685.1_ASM2568v1/GCF_000025685.1_ASM2568v1_genomic.gff |\
  sed 's/;/\t/g' | grep "ID=gene" |\
    sed 's/ID=//g; s/Dbxref=GeneID://g; s/Name=//g; s/gbkey=//g; s/gene_biotype=//g; s/locus_tag=//g; s/old_locus_tag=//g' > hvo_gff.txt
## in excel - tidy up or just remove all pseudogenes with grep
  #sed  gene_lengths.txt > genes.tmp
  #awk -F "\t" '{OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' genes.tmp > gene_lengths.txt

## make clean files for complete mapping
grep "ID=gene" GCF_000025685.1_ASM2568v1_genomic.gff | \
  sed 's/end_range=.*\.;//; s/partial=true;//; s/;pseudo=true//; s/;start_range=.*//; s/Name=.*\;gbkey/gbkey/; s/gene=.*\;gene_biotype/gene_biotype/; s/;/\t/g' | \
  sed 's/ID=gene//g; s/Dbxref=GeneID://g; s/Name=//g; s/gbkey=//g; s/gene_biotype=//g; s/old_locus_tag=//g; s/locus_tag=//g' > gff_clean.txt

grep ">" GCF_000025685.1_ASM2568v1_cds_from_genomic.fna | \
  sed 's/\[gene.*] \[locus/ \[locus/' | \
  awk '{OFS="\t"}{print $1,$2,$3}' > cds_ids.txt
  # used textwrangler to find replace >, [, and ] because I'm lazy - need to write into shell script for future easy manipulating

#grep ">" hfx_volc_ref/GCF_000025685.1_ASM2568v1/GCF_000025685.1_ASM2568v1_cds_from_genomic.fna | \
#  grep "_WP_" | \
#  sed 's/ \[/\t\[/g' | \
#  awk -F"\t" '{OFS="\t"} {print $1,$6}' | \
#  sed 's/\[location=complement(//g; s/\[location=//g; s/)\]//g; s/\.\./\t/g; s/\]//g' > newlengths.txt


# from GCF_000025685.1_ASM2568v1_assembly_report.txt
# Sequence-Name Sequence-Role   Assigned-Molecule       Assigned-Molecule-Location/Type GenBank-Accn    Relationship    RefSeq-Accn     Assembly-Unit   Sequence-Length UCSC-style-name
#pHV1    assembled-molecule      pHV1    Plasmid CP001957.1      =       NC_013968.1     Primary Assembly        85092   na
#pHV2    assembled-molecule      pHV2    Plasmid CP001954.1      =       NC_013965.1     Primary Assembly        6359    na
#pHV3    assembled-molecule      pHV3    Plasmid CP001953.1      =       NC_013964.1     Primary Assembly        437906  na
#pHV4    assembled-molecule      pHV4    Plasmid CP001955.1      =       NC_013966.1     Primary Assembly        635786  na
#ANONYMOUS       assembled-molecule      na      Chromosome      CP001956.1      =       NC_013967.1     Primary Assembly        2847757 na

# from GCF_000025685.1_ASM2568v1_cds_from_genomic.fna
#>lcl|NC_013968.1_cds_WP_013035720.1_5 [locus_tag=HVO_RS19310] [db_xref=GeneID:8926954] [protein=hypothetical protein] [protein_id=WP_013035720.1] [location=4640..5341] [gbkey=CDS]

## using the entrez gene id, acquired GO terms and other data from david
## https://david.ncifcrf.gov/conversion.jsp
```

thiberio http://projetos.lbi.iq.usp.br/hgt/interactive_plots.html











### end ###
