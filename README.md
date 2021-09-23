# Sarek Installation and Documentation
Installation and setup documents for sarek pipeline

## Introduction
Sarek is a pipeline built in Nextflow designed to detect variants on whole genome or targeted sequencing data. Sarek can also handle tumor/normal pairs…
More information can be found at: 
https://nf-co.re/sarek

## Installing & Setup

1. [Install Nextflow](#install-nextflow)
2. [Install Singularity](#install-singularity)
3. [Pull or Clone Sarek](#pull-or-clone-sarek)
4. [Get & configure resources](#get-and-configure-resources)
5. Configure Nextflow
6. Configure processes
7. Modify code
8. [Available tools](#available-tools)

## Install Nextflow

Install the latest version of Nextfow (v20.10+ is needed for the latest v2.7.1 sarek release) using the [these instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation)

## Install Singularity

Singularity is most likely already available on the HPC, if not follow [these instructions](https://sylabs.io/guides/3.5/admin-guide/)

## Pull or Clone Sarek

Use the following command to clone the Sarek pipeline in a sarek subfolder in the current working directory:

```
git clone https://github.com/nf-core/sarek.git
```

## Get and Configure resources

The default reference genome source for Sarek is [Illumina's iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html) which can alternatively be installed using [Amazon Web Services (AWS)](https://ewels.github.io/AWS-iGenomes)

It is also possible to use specific local genome resources.

### Downloading reference genome directly (GRCh37 example)

Download the reference genome of choice from [https://support.illumina.com/sequencing/sequencing_software/igenome.html](https://support.illumina.com/sequencing/sequencing_software/igenome.html) and extract the downloaded archive

```
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz
tar -xzvf Homo_sapiens_Ensembl_GRCh37.tar.gz
```

### Downloading reference genome using AWS (GRCh37 example)

Download the reference genome of choice from https://ewels.github.io/AWS-iGenomes/<br><br>
For Genome -> Source -> Build -> Type choose:<br>
Homo_Sapiens -> GATK -> GRCh37

Set up AWS for ad-hoc usage using guix

```
guixr environment --ad-hoc awscli
```

Change to target folder
```
cd /hpc/ubec/projects/sarek_v2.7.1_rollout/resources
```

Download selected genome (using the specified Sync command
```
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/ ./Homo_sapiens/GATK/GRCh37
```

### Customize config files

In order to use the downloaded reference genome, the config files in the sarek pipeline must be adjusted

For example, to use and specify a custom config with the name 'umcu_grch37'<br>

#### Create custom config file

First we create a conf_custom subfolder in the sarek folder where the custom configs will be stored<br>
Next in that folder we create a new file the name 'umcu_grch37.config' with the following contents<br>

```
//Profile config names for nf-core/configs
params {
  config_profile_description = 'UMCU HPC Cluster GRCh37 config.'
  config_profile_contact = 'Firstname Lastname'
  config_profile_url = 'https://www.umcutrecht.nl/'

  save_reference = true

  // illumina iGenomes reference file paths on UMCU HPC
  igenomes_base = '/hpc/ubec/resources/igenomes/'
  genome = 'GRCh37'

  //keep mapped bams by default (preferable in case pipeline needs to be rerun)
  save_bam_mapped = true
  
  //increased mem usage due to step crashing otherwise
  markdup_java_options = '"-Xms16g -Xmx32g"' // Established values for markDuplicates memory consumption, see https://github.com/SciLifeLab/Sarek/pull/689 for details

  email = 'name@mail.nl'
}

//UMCU HPC Singularity settings
singularity {
  enabled = true
  autoMounts = true
  runOptions = '-B /hpc -B $TMPDIR:$TMPDIR'
  cacheDir = '/hpc/local/CentOS7/cog_bioinf/singularity_cache'
}

//resources customized for processes
process {
  executor = 'slurm'

  //customized resource settings, fine tune as required
   withName : 'get_software_versions' {
      time = '1h'
      memory = '2G'
   }
   
   withName : 'FastQCFQ' {
      time = '8h'
      memory = '8G'
   }
   
   withName : 'MarkDuplicates' {
      time = '48h'
   }
   
   withName : 'VEP' {
      time = '48h'
      memory = '16G'
   }  
   withName : 'VEP' {
      time = '48h'
      memory = '16G'
   }  
}

//Tells nextflow to cleanup work directory after succesful run
cleanup = true
```

#### Modify nextflow.config

The next step is to modify the nextflow.config in order for sarek to be able to properly find and use the custom config file. The nextflow.config file can be found in the sarek folder and in order for custom config files to properly work a line must be added for each specific custom config file in the section: 

```
profiles {
...
  test                { includeConfig 'conf/test.config' }
  test_annotation     { includeConfig 'conf/test_annotation.config' }
  test_use_gatk_spark { includeConfig 'conf/test_use_gatk_spark.config' }
...
}
```

Add the following line to that section in order to tell sarek where to find the config we just created:

```
umcu_grch37         { includeConfig 'conf_custom/umcu_grch37.config' }
```
This tells sarek that the config can be found in the 'conf_custom' subfolder and has the filename 'umcu_grch37.config'

## Available tools
* ASCAT<br>
https://github.com/VanLoo-lab/ascat<br><br>
Infer tumor purity, ploidy and allele specific copy number profiles

* CNVkit*<br>
https://github.com/etal/cnvkit<br><br>
Detecting copy number variants and alterations genome wide from high throughput sequencing

* ControlFREEC*<br>
https://github.com/BoevaLab/FREEC<br><br>
Copy number and genotype annotation for whole genome and whole exome sequencing data. Detecting copy number changes and imbalances

* FreeBayes<br>
https://github.com/freebayes/freebayes<br><br>
A haplotype-based variant detector for small polymorphisms (snps, indels mnp and complex events smaller then the length of a short-read sequencing alignment

* HaplotypeCaller<br>
https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller<br><br>
Call germline SNPs and indels via local re-assembly of haplotypes

* Manta<br>
https://github.com/Illumina/manta<br><br>
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads and is optimized for analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.

* Mpileup<br>
http://www.htslib.org/doc/samtools-mpileup.html<br><br>
produces "pileup" textual format from an alignment

* MSIsensor<br>
https://github.com/ding-lab/msisensor<br><br>
Microsatellite instability detection user tumor only or paired tumor-normal data

*	Mutect2<br>
https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2<br><br>
Call somatic SNVs and indels via local assembly of haplotypes

*	Strelka2<br>
https://github.com/Illumina/strelka<br><br>
Small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal pair samples

*	TIDDIT<br>
https://github.com/SciLifeLab/TIDDIT<br><br>
A tool used to identify chromosomal rearrangements using mate pair or paired end sequencing data. TIDDIT identifies intra and inter-chromosomal translocations, deletions, tandem-duplications and inversions, using supplementary alignments as well as discordant pairs.

*	snpEff<br>
http://pcingola.github.io/SnpEff<br><br>
Genetic variant annotation and functional effect prediction toolbox. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes).

*	VEP<br>
https://github.com/Ensembl/ensembl-vep<br><br>
(Variant Effect Predictor) predicts the functional effects of genomic variants.

*	Merge ???<br>
