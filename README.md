# Qiime_steps

# QIIME 2 tutorial for 16S rRNA or ITS data analysis
The only two things that need to be changed:
1. Primer sequences in "Step2"
2. Database along with the primers in "Step6".

Rest everything will be the same in the analysis of both kind of datasets.

## This tutorial can be used on Demultiplexed Paried end sequences
Deatiled instructions can be followed from https://docs.qiime2.org/2021.4/tutorials/moving-pictures/

### Load Qiime
``` r
module load miniconda/qiime/3/2019.7
```
### Step1: Import sequences
Before running ** qiime tools import ** user need to create a Manifest file: a simple csv file wich have three columns. Check the formatting of Manifest file
``` r
head Manifest.csv
sample-id,absolute-filepath,direction
IS013,/home/gomeza/sharm646/Raw/IS013_S198_R1_001_filtered.fastq.gz,forward
IS013,/home/gomeza/sharm646/Raw/IS013_S198_R2_001_filtered.fastq.gz,reverse
IS014,/home/gomeza/sharm646/Raw/IS014_S205_R1_001_filtered.fastq.gz,forward
IS014,/home/gomeza/sharm646/Raw/IS014_S205_R2_001_filtered.fastq.gz,reverse
IS015,/home/gomeza/sharm646/Raw/IS015_S212_R1_001_filtered.fastq.gz,forward
IS015,/home/gomeza/sharm646/Raw/IS015_S212_R2_001_filtered.fastq.gz,reverse
```
This manifest file will be used to import these sequence data files into a QIIME 2 artifact. The final output **.qzv** artifcat can be visualized at https://view.qiime2.org

``` r
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Manifest.csv --output-path demux.qza --input-format PairedEndFastqManifestPhred33
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
```

### Step2: Remove primers using cutadapt
It is important to know which primers were used for the amplicaiton. Here, we have used forward and reverse primer sequences of 16S V4 region. These primer sequences need to be changed depending upon the study.
``` r
qiime cutadapt trim-paired --i-demultiplexed-sequences demux.qza --p-cores 20 --p-front-f GTGYCAGCMGCCGCGGTAA --p-front-r GGACTACNVGGGTWTCTAAT --p-overlap 10 --p-discard-untrimmed False --o-trimmed-sequences trimmed-seqs.qza --verbose
qiime demux summarize --i-data trimmed-seqs.qza --o-visualization trimmed-seqs.qzv
```
### Step3: Feature table construction using Dada2
Upload trimmed-seqs.qzv file at https://view.qiime2.org to check the quality of basepairs in foward and reverse reads and decide the length of the read that should be trimmed based on **q30 scores**. 
``` r
qiime dada2 denoise-paired --verbose --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 240 --p-trunc-len-r 180 --o-table table.qza --o-representative-sequences rep-seqs.qza --p-n-threads 12 --o-denoising-stats dada-stats.qza
```

### Step4: Generate a tree for phylogenetic diversity analyses
Align and mask (elimination of phylogenetically uninformative or misleading sites) the alignments.
``` r
qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza
qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
```
Construct both unrooted and rooted phylogenetic trees
``` r
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza
```

### Step5: Alpha and beta diversity analysis - Detiled anlysis will be done in R on final output files.
To quickly access the patterns in the data Alpha and beta diversity analysis can be performed at a rarefy depth of min 1000 (User can select this as per the sequecning depth in their own dataset)
``` r
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth 1000 --m-metadata-file Metadata.tsv --output-dir core-metrics-results
```

### Step6: Taxonomic Assignment
First step is to download the database along with the taxonomy file. Here is an example of GreenGenes database which was downloaded from https://docs.qiime2.org/2021.4/data-resources/
**Or** User can download the pretrained classifiers also but they are available only for specific regions such v4 16S.
For ITS of Fungal analysis Unite database can be downloaded and processed using the same steps or follow steps here to train the classifier
http://john-quensen.com/tutorials/training-the-qiime2-classifier-with-unite-its-reference-sequences/

1. Next we import these data into QIIME 2 Artifacts
``` r
qiime tools import --type 'FeatureData[Sequence]' --input-path /home/sharmaa4/Databases/gg_13_8_otus/rep_set/99_otus.fasta --output-path 99_otus.qza 
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path /home/sharmaa4/Databases/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt --output-path ref-taxonomy.qza 
```
2. Extract reference reads using primer sequences
``` r
qiime feature-classifier extract-reads --i-sequences 99_otus.qza --p-f-primer GTGYCAGCMGCCGCGGTAA --p-r-primer GGACTACNVGGGTWTCTAAT --o-reads ref-seqs.qza 
```
3. Train the classifier
``` r
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs.qza --i-reference-taxonomy ref-taxonomy.qza --o-classifier classifier.qza 
```
4. Predict taxonomy using Trained Classifier
``` r
qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
```
5. Make taxa table at different taxonomic levels
``` r
Speices: qiime taxa collapse --i-table table.qza  --i-taxonomy taxonomy.qza --p-level 7 --o-collapsed-table table-l7.qza
Genus: qiime taxa collapse --i-table table.qza  --i-taxonomy taxonomy.qza --p-level 6 --o-collapsed-table table-l6.qza
Family: qiime taxa collapse --i-table table.qza  --i-taxonomy taxonomy.qza --p-level 5 --o-collapsed-table table-l5.qza
Phylum: qiime taxa collapse --i-table table.qza  --i-taxonomy taxonomy.qza --p-level 2 --o-collapsed-table table-l2.qza
```
6. Unzip qza files and then convert them in txt format using follwoing commands
``` r
unzip table.qza
biom convert -i feature-table.biom -o feature_table.txt --to-tsv
```

### Step7: Final files- that will be used to contruct the Phyloseq object
table.qza, rooted-tree.qza, taxonomy.qza, Metadata.tsv can be used to contruct the phyloseq object and for downstream statistical analysis using Phyloseq and DESeq2

Use the following steps
``` r
library(qiime2R)
library(phyloseq)
physeq<-qza_to_phyloseq(features="data/table.qza", tree="data/rooted-tree.qza", taxonomy="data/taxonomy.qza", metdata="data/sample-metadata.qza")
```
**Or** another way to Constuct the phyloseq object:
1. Load the ASV table
``` r
asv_table_in <- read.csv("feature_table.txt", sep = "\t", row.names = 1)
asv_table_in_t <- data.frame(t(asv_table_in))
asv_table_in_t <- asv_table_in_t[order(rownames(asv_table_in_t)),]
```
Additional steps to filter ASVs- ASV sum and in how many samples that ASV should be present can be changed as per the data
``` r
library (labdsv)
asv.sums <- colSums(asv_table_in_t)
asv_table_filtered <- asv_table_in_t[ , which(asv.sums > 100)] #---Check ASV sum should be more than 100
asv_table_filtered_filtered <- dropspc(asv_table_filtered, 10) #--- ASV should be present in more than 10 samples
asv_table_filtered_t <- data.frame(t(asv_table_filtered_filtered))
asv_table_in <- as.matrix(asv_table_filtered_t)
```
2. Read in taxonomy
``` r
taxonomy <- read.csv("taxonomy.tsv", sep = "\t", row.names = 1)
taxonomy <- as.matrix(taxonomy)
```
3. Read in metadata
``` r
metadata <- read.table("Metadata.txt", sep="\t", row.names = 1, header=T)
metadata_t <- data.frame(t(metadata))
metadata <- data.frame(t(metadata_t))
metadata <- metadata[order(rownames(metadata)),]
```
4. Read in tree
``` r
phy_tree <- read_tree("tree.nwk")
```
5. Check everything should be in an order
``` r
ASV <- otu_table(asv_table_in, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
META <- sample_data(metadata)
```
Check for consistent OTU names
``` r
taxa_names(TAX)
taxa_names(ASV)
taxa_names(phy_tree)
```
Check for consistent sample names
``` r
sample_names(ASV)
sample_names(META)
```
Finally merge to create Phyloseq object!
``` r
ps <- phyloseq(ASV, TAX, META, phy_tree)
```

### Step8: Subsetting, rarefaction, removal of samples less than particular depth without rarefaction and normalization in the Phyloseq Object
1. Subsetting phyloseq object- based on any metadata column here For example we are extracting based on Sample_Type
``` r
ps_filtered <- subset_samples(ps, Group%in%c("Controls", "Cases"))
```
2. Either: Rarefy data at a particular depth- can be decided by ASV counts in each sample
``` r
ps_filtered.rare = rarefy_even_depth(ps_filtered, sample.size = 1000)
```
3. Or: Don't rarefy only remove samples having less depth here example is 1000
``` r
ps_filtered.prune1k = prune_samples(sample_sums(ps_filtered) > 1000, ps_filtered)
```
4. Normalize to calculate realtive abundances- Either the output of Step2 and Step3 can be used for normalization
``` r
ps_normalized  = transform_sample_counts(ps_filtered.rare, function(x) x / sum(x) )
```

### Step8: Alpha, beta and taxonomic analysis
1. Alpha diversity analyis
``` r
erich <- estimate_richness(ps_filtered, measures = c("Observed", "Chao1", "Shannon"))
```
2. Beta diversity analyis
``` r
library (vegan)
library (ape)
beta_div_dist <- distance(physeq = ps_normalized, method = "bray")
beta_div_dist <- as(beta_div_dist, "matrix")
bray_pcoa <-pcoa(beta_div_dist)
bray_pcoa$values[1:2,]
```
2. Taxonomic analyis: Convert phyloseq object in to DESeq
``` r
diagdds = phyloseq_to_deseq2(ps_filtered, ~ Group)
diagdds <- estimateSizeFactors(diagdds, type = 'poscounts')

diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_filtered)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
```
For taxonomic - indval analysis can be done using labdsv package or important taxa can be identified using randomForest package in R.
