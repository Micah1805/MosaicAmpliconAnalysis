# Script: 01_generate_asv_tables_prokaryotes.R
# Purpose: Generates ASV (Amplicon Sequence Variant) tables for prokaryotic samples.
# Input: Fastqs
# Main output: taxonomy_prokaryotes.csv, asv_table_prokaryotes.csv
# Author: [Insert your name]
# Date: [Insert date or leave blank]
# -----------------------------------------------------

#This script prepares amplicon sequence (ASV) tables 
#from Fastq reads for each sample 
#Input: 16S Fastq reads (forward and reverse redas)
#Output: (Amplicon) sequence table seqtab2

####Housekeeping####

#Load packages
library(dada2)
packageVersion("dada2")
library(tidyverse)
library(ape)
library(readxl)
library(stringr)
library(ShortRead)
library(Biostrings)

#Prepare environment
#Define project and database directories
#Define working directory as project directory

project_dir <- "./data"
path="./data/16S"

#Define runs as input
#Define all relevant analysed runs (each directory contains forward and reverse 16S fastq reads)
#If all fastqs are in one directory:
runID.1 <- file.path(path,"/Fastqs")

#If fastqs are in multiple subdirectories (replace dir_x by subdirectory name):
#runID.1 <- file.path(path,"/Fastqs/dir_x")
#runID.2 <- file.path(path,"/Fastqs/dir_x")
#runID.3 <- file.path(path,"/Fastqs/dir_x")
#runID.4 <- file.path(path,"/Fastqs/dir_x")

#If all fastqs are in one directory:
runs <- list(runID.1)

#If there are multiple subdirectories
#runs <- list(runID.1, runID.2, runID.3, runID.4)

#Set sequencing run to denoise (name of sequencing run in raw data directory)
seq_run <- runs

#Raw data to work with
seq_dir <- file.path(seq_run)

#Get timestamp to distinguish between analysis runs
time.stamp <- gsub(" ","_",gsub(":","-",Sys.time()))

#Pathes and directories for intermediate results 
#Construct needed pathes
temp_dir <- file.path(path,paste(time.stamp,"seq_runs","out",sep="_"))
preFilt_dir <- file.path(temp_dir,"preFilt")
filtered <- file.path(temp_dir,"filterted")
path.cut <- file.path(temp_dir, "cutadapt")

#Create directories
if(!dir.exists(temp_dir)) dir.create(temp_dir)
if(!dir.exists(preFilt_dir)) dir.create(preFilt_dir)
if(!dir.exists(filtered)) dir.create(filtered)
if(!dir.exists(path.cut)) dir.create(path.cut)

####Check quality of raw data####
#Fast check if all samples are included
#Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(seq_dir, pattern="_R1_001.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(seq_dir, pattern="_R2_001.fastq.gz$", full.names = TRUE))  

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnRs), "_S"), `[`, 1)
sample.names

#Construct needed file lists
fnFs.preFilt <- file.path(preFilt_dir,basename(fnFs))
fnRs.preFilt <- file.path(preFilt_dir,basename(fnRs))
fnFs.filtered <- file.path(filtered,basename(fnFs))
fnRs.filtered <- file.path(filtered,basename(fnRs))
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

#Plot Phred quality score profiles of subset of raw files; data are aggregated; always subset to reduce computational costs
#Forward reads
Qual_raw_Fs <- plotQualityProfile(fnFs[1:2],aggregate=TRUE)
#Reverse reads
Qual_raw_Rs <- plotQualityProfile(fnRs[1:2],aggregate=TRUE)

ggsave("Quality_raw_Fs.pdf",Qual_raw_Fs, device="pdf",path=temp_dir)
ggsave("Quality_raw_Rs.pdf",Qual_raw_Rs, device="pdf",path=temp_dir)


####Remove primers using cutadapt and prefilter data####
#Define primer
FWD_PRIMER="GTGYCAGCMGCCGCGGTAA"
REV_PRIMER="CCGYCAATTYMTTTRAGTTT"

#Create reverse complements of primer sequences
FWD_PRIMER.RC <- dada2:::rc(FWD_PRIMER)
REV_PRIMER.RC <- dada2:::rc(REV_PRIMER)

#Prepare primer checks; create all possible primer versions (fwd, fwd_complement, rev, rev_complement)
allOrients <- function(primer) {                                                # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)                                                      # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), 
               Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  
}

FWD_PRIMER.orients <- allOrients(FWD_PRIMER)
REV_PRIMER.orients <- allOrients(REV_PRIMER)

FWD_PRIMER.orients 
REV_PRIMER.orients 



##Prefiltering (remove sequences with Ns to increase accuracy of mapping short primer sequences)

prefilterOut <- filterAndTrim(fnFs,fnFs.preFilt,fnRs,fnRs.preFilt,
                              minLen=100,multithread = FALSE)
prefilterOut


#Count number of times  the primers appear in forward and reverse read
primerHits <- function(primer, fn) {                                            # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, 
      fn = fnFs.preFilt[[1]]), 
      FWD.ReverseReads = sapply(FWD_PRIMER.orients,
      primerHits, fn = fnRs.preFilt[[1]]), 
      REV.ForwardReads = sapply(REV_PRIMER.orients, primerHits,
      fn = fnFs.preFilt[[1]]), REV.ReverseReads = sapply(REV_PRIMER.orients, 
      primerHits, fn = fnRs.preFilt[[1]]))


#Remove all primers with cutadapt
#Define location of cutadapt (needs to be installed locally if not installed on a server)
cutadapt <- "cutadapt.exe" 
#wWich cutadapt version? (check whether cutadapt is found)
system2(cutadapt, args = "--version")


#Remove primers at 5'-end with cutadapt
#Primer sequences
FWD_PRIMER
REV_PRIMER

#Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_PRIMER, "-a", REV_PRIMER.RC) 
#Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_PRIMER, "-A", FWD_PRIMER.RC) 


for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,                       # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i],              # output files
                             fnFs.preFilt[i], fnRs.preFilt[i]))                 # input files
}



rbind(FWD.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, 
      fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD_PRIMER.orients,
      primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV_PRIMER.orients, primerHits,
      fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV_PRIMER.orients, 
                                primerHits, fn = fnRs.cut[[1]]))

#Forward and reverse Fastq filenames have the format
cutFs <- sort(list.files(path.cut, pattern ="_L001_R1_001.fastq.gz", 
                         full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern ="_L001_R2_001.fastq.gz", 
                         full.names = TRUE))

#Extract sample names, assuming filenames have format
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#Read quality
Qual_CutFS<-plotQualityProfile(cutFs[1:2],aggregate=TRUE)
Qual_CutRS<-plotQualityProfile(cutRs[1:2],aggregate=TRUE)

ggsave("Qual_CutFS.pdf",Qual_CutFS, device="pdf",path=temp_dir)
ggsave("Qual_CutRS.pdf",Qual_CutRS, device="pdf",path=temp_dir)

####Filter and trimm and dereplicate####

#Place filtered files in filtered/ subdirectory

filtFs <- file.path(filtered,basename(cutFs))
filtRs <- file.path(filtered,basename(cutRs))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the 
#maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.

#maxEE define quality score. lower = better quality
#To be conservative use 1.5. this means reads with average quality score of 25 and sd 4.5 will be discarded
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,#truncLen = c(275,200),
                     maxN=0,maxEE=c(2.60,2.10),
                     truncLen=c(260,210),
                     verbose = TRUE, rm.phix = TRUE, 
                     compress = TRUE, multithread = FALSE)


#Sequences lost by quality filtering
out   


#Dereplicate identical sequences to reduce computational costs
fnFs.deRep <- derepFastq(fnFs.filtered)
fnRs.deRep <- derepFastq(fnRs.filtered)
names(fnFs.deRep) <- sample.names
names(fnRs.deRep) <- sample.names



#### Learn the Error Rates####
errF <- learnErrors(fnFs.deRep, multithread=TRUE)
errR <- learnErrors(fnRs.deRep, multithread=TRUE)
plotErrors_F <- plotErrors(errF, nominalQ=TRUE)
plotErrors_R <- plotErrors(errR, nominalQ=TRUE)

ggsave("errorplot_F.pdf",plotErrors_F, device="pdf",path=temp_dir)
ggsave("errorplot_R.pdf",plotErrors_R, device="pdf",path=temp_dir)


####Sample Interference####
#Apply the core sample inference algorithm to the filtered and trimmed sequence data.
#Check how many reads per unique sequence

dadaFs <- dada(fnFs.deRep, err=errF, multithread=TRUE)
dadaRs <- dada(fnRs.deRep, err=errR, multithread=TRUE)


####Merge paired reads (merge forward and reverse reads)####

mergers <- mergePairs(dadaFs, fnFs.deRep, dadaRs, fnRs.deRep, 
                      verbose=TRUE,minOverlap = 15)

#Inspect the merger data.frame from the first sample
head(mergers[[1]])


####Construct ASV table####
#construct an amplicon sequence variant table (ASV) table, a higher-resolution 
#version of the OTU table produced by traditional methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
ht=hist(nchar(colnames(seqtab)),100)

#Expected sequence length after merging considering the overlap von 12 ~20 bp
#Consider sequence length as reliable if > than 500 samples 
#Get length with max abundance
sldist=names(which.max(table(nchar(getSequences(seqtab)))))
(sldist)

#Take a deviation of 1% left and right from max value
lft=as.numeric(sldist)-0.025*as.numeric(sldist)
lft=round(lft,digits = 0)
rgt=as.numeric(sldist)+0.025*as.numeric(sldist)
rgt=round(rgt,digits = 0)
show(c(lft,rgt))
abline(v=c(lft,rgt),col=2,lwd=2)

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% lft:rgt]
table(nchar(getSequences(seqtab2)))
save(seqtab2,file = paste0(temp_dir,"/seqtab_Prok.Rdata"))

#### Remove chimeras #####
seqtab.nochim <- removeBimeraDenovo(seqtab2,method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#Frequency of non-chimeric sequences
sum(seqtab.nochim)/sum(seqtab2)


#####Track reads through the pipeline####
#As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: 
#e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- 
  c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

save(seqtab.nochim,file = paste0(temp_dir,"/seqtab_nochim_Prok.Rdata"))


#### TAXA assignment####
#Data for taxonomy and species assignment need to be downloaded from Silva and
#stored locally (whwn not working on a server)
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", 
multithread=TRUE)

taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")

#Inspect taxonomic assignments
taxa.print <- taxa                                                              # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

save(taxa,file=paste0(temp_dir,"/TMPtaxa_Prok.Rdata"))

prok_taxa=taxa

####Write OTU sequences to fasta file####
seq=colnames(seqtab.nochim)
names(seq)=paste0("prok_asv",seq(length(seq)))
fa=sapply(seq,function(x)strsplit(x,''))
fa=as.DNAbin.list(fa)
write.FASTA(fa,file=paste0(temp_dir,"/otu_seq_Prok.fa"))

#Create R object with taxa and abund mat file
mi=match(rownames(taxa),seq)
rownames(prok_taxa)=names(seq)[mi]
save(prok_taxa,file=paste0(temp_dir,"/taxa_Prok.Rdata"))
write.csv(prok_taxa,file="./output/taxonomy_prokaryotes.csv", row.names=TRUE)

abundMatRaw=seqtab.nochim
mi=match(colnames(seqtab.nochim),seq)
colnames(abundMatRaw)=names(seq)[mi]
abundMatRaw=t(abundMatRaw)
save(abundMatRaw,file=paste0(temp_dir,"/asv_table_prokaryotes.Rdata"))
write.csv(abundMatRaw,file="./output/asv_table_prokaryotes.csv", row.names=TRUE)

