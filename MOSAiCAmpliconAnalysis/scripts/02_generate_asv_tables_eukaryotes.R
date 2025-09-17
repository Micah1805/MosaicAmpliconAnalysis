# Script: 02_generate_asv_tables_eukaryotes.R
# Purpose: Generates ASV (Amplicon Sequence Variant) tables for eukaryotic samples.
# Input: (Specify expected input file(s) in data/)
# Output: (Specify generated output file(s) in output/)
# Author: [Insert your name]
# Date: [Insert date or leave blank]
# -----------------------------------------------------

################################################################################
######Generate amplicon sequence (ASV) tables for eukaryotes using DADA2#######
################################################################################

#This script prepares amplicon sequence (ASV) tables 
#from Fastq reads for each sample 
#Input: 18S Fastq reads (forward and reverse redas)
#Output: (Amplicon) sequence table seqtab2

####Preparation####
##Prepare environment
#define project and database directories
#define working directory as project directory
project_dir <- file.path(getwd())
path="C:/Users/jahuem001/OneDrive/Dokumente/AWI/MOSAiC/Sequencing_Prok_Euk/Fastqs/Euks/Water_only"
#set software pathes
#overwrite personal R library pathes: definition of system library pathes and 
#common project library path
.libPaths(c('/usr/local/lib/R/site-library',
            'C:/Users/jahuem001/AppData/Local/R/win-library/4.3'
            ,'C:/Program Files/R/R-4.3.2/library'))

#Load packages
library(dada2)
packageVersion("dada2")
library(tidyverse)
library(ape)
library(readxl)
library(stringr)
library(ShortRead)
library(Biostrings)

#Define runs as input
#Define all relevant analysed runs
#<time_stamp>_<seq_run>
#example: "2020-10-27_13-31-04_M03457_0049_000000000-CRH8D"
runID.1 <- file.path(path,"/Lauf_1")
#runID.2 <- file.path(path,"/Lauf_2")
#runID.3 <- file.path(path,"/Lauf_3")
#runID.4 <- file.path(path,"/M03457_0091_000000000-JV6HN")
#runID.5 <- file.path(path,"/M03690_0336_000000000-JLGJW")
#runID.6 <- file.path(path,"/M03690_0339_000000000-JN8C7")

runs<-list(runID.1)
           #,runID.2,runID.3,runID.4,runID.5,runID.6)

#Set sequencing run to denoise (name of sequencing run in raw data directory)
seq_run <- runs

#Raw data to work with
seq_dir <- file.path(seq_run)

#Get timestamp to distinguish between analysis runs
time.stamp <- gsub(" ","_",gsub(":","-",Sys.time()))


##Pathes and directories for intermediate results
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
# Forward and reverse fastq filenames have format: 
#SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <-sort(list.files(seq_dir, pattern="_R1_001.fastq.gz$", full.names = TRUE))
fnRs <-sort(list.files(seq_dir, pattern="_R2_001.fastq.gz$", full.names = TRUE))  

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

#Plot Phred quality score profiles of subset of raw files; 
#data are aggregated; always subset to reduce computational costs
#Forward reads
Qual_raw_Fs<-plotQualityProfile(fnFs[1:2],aggregate=TRUE)
#Reverse reads
Qual_raw_Rs<-plotQualityProfile(fnRs[1:2],aggregate=TRUE)

ggsave("Quality_raw_Fs.pdf",Qual_raw_Fs, device="pdf",path=temp_dir)
ggsave("Quality_raw_Rs.pdf",Qual_raw_Rs, device="pdf",path=temp_dir)

####Remove primers using cutadapt and prefilter data####
#Define primer sequences
FWD_PRIMER="GCGGTAATTCCAGCTCCAA"
REV_PRIMER="ACTTTCGTTCTTGATYRR"

#Create reverse complements of primer sequences
FWD_PRIMER.RC <- dada2:::rc(FWD_PRIMER)
REV_PRIMER.RC <- dada2:::rc(REV_PRIMER)

#PREPARE primer checks; create all possible primer versions 
#(fwd, fwd_complement, rev, rev_complement)
allOrients <- function(primer) {                                                #Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)                                                      #The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), 
               Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  
}

FWD_PRIMER.orients <- allOrients(FWD_PRIMER)
REV_PRIMER.orients <- allOrients(REV_PRIMER)

FWD_PRIMER.orients 
REV_PRIMER.orients 


#Prefiltering 
#(remove sequences with Ns to increase accuracy of mapping short primer sequences)

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


#Removeall PRIMERS with cutadapt
#Define location of cutadapt 
#(needs to be installed locally if not instlled on a server)
cutadapt <- "C:/Users/jahuem001/OneDrive/Dokumente/AWI/MOSAiC/Sequencing_Prok_Euk/Fastqs/cutadapt.exe" 
#Wich cutadapt version? (check whether cutadapt is found)
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
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,                       #-n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i],              #Output files
                             fnFs.preFilt[i], fnRs.preFilt[i]))                 #Input files
}



rbind(FWD.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, 
      fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD_PRIMER.orients,
      primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV_PRIMER.orients, primerHits,
      fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV_PRIMER.orients, 
      primerHits, fn = fnRs.cut[[1]]))

#Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern ="_L001_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern ="_L001_R2_001.fastq.gz", full.names = TRUE))

#Extract sample names, assuming filenames have format:
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

# standard filtering parameters: maxN=0 (DADA2 requires no Ns), 
#truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the 
#maximum number of “expected errors” allowed in a read, which is 
#a better filter than simply averaging quality scores.

## maxEE define quality score. lower = better quality
## to be conservative use 1.5. this means reads with average quality score of 
#25 and sd 4.5 will be discarded
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
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

####Learn the Error Rates####
errF <- learnErrors(fnFs.deRep, multithread=TRUE)
errR <- learnErrors(fnRs.deRep, multithread=TRUE)
plotErrors_F<-plotErrors(errF, nominalQ=TRUE)
plotErrors_R<-plotErrors(errR, nominalQ=TRUE)

ggsave("errorplot_F.pdf",plotErrors_F, device="pdf",path=temp_dir)
ggsave("errorplot_R.pdf",plotErrors_R, device="pdf",path=temp_dir)

####Sample Interference####
#Apply the core sample inference algorithm to the 
#filtered and trimmed sequence data.
#Check how many reads per unique sequence

dadaFs <- dada(fnFs.deRep, err=errF, multithread=TRUE)
dadaRs <- dada(fnRs.deRep, err=errR, multithread=TRUE)

####Merge paired reads (merge forward and reverse reads)####
mergers <- mergePairs(dadaFs, fnFs.deRep, dadaRs, fnRs.deRep, 
                      verbose=TRUE,minOverlap = 15)

#Inspect the merger data.frame from the first sample
head(mergers[[1]])

####Construct ASV table####
#Construct an amplicon sequence variant table (ASV) table, a higher-resolution 
#version of the OTU table produced by traditional methods
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
ht=hist(nchar(colnames(seqtab)),100
        
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
save(seqtab2,file = paste0(temp_dir,"/seqtab_Euk.Rdata"))

#### Remove chimeras #####
seqtab.nochim <- removeBimeraDenovo(seqtab2,method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#frequency of non-chimeric sequences
sum(seqtab.nochim)/sum(seqtab2)

#####Track reads through the pipeline####
#As a final check of our progress, we’ll look at the number of reads
#that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

save(seqtab.nochim,file = paste0(temp_dir,"/seqtab_nochim_Euk.Rdata"))

####Taxa assignmennt####
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/jahuem001/OneDrive/Dokumente/AWI/MOSAiC/Sequencing_Prok_Euk/Fastqs/Euks/Pr2/pr2_version_5.0.0_SSU_dada2.fasta.gz",multithread = T,
                       taxLevels = c("Kingdom","Supergroup","Division","Subdivision","Class","Order","Family","Genus","Species"))

#Inspect taxonomic assignments
taxa.print <- taxa                                                              #Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


save(taxa,file=paste0(temp_dir,"/TMPtaxa_Euk.Rdata"))

euk_taxa=taxa

####Write OTU sequences to fasta file####
seq=colnames(seqtab.nochim)
names(seq)=paste0("prok_asv",seq(length(seq)))
fa=sapply(seq,function(x)strsplit(x,''))
fa=as.DNAbin.list(fa)
write.FASTA(fa,file=paste0(temp_dir,"/otu_seq_Euk.fa"))

#Create R object with taxa and abund mat file
mi=match(rownames(taxa),seq)
rownames(prok_taxa)=names(seq)[mi]
save(prok_taxa,file=paste0(temp_dir,"/taxa_Euk.Rdata"))

abundMatRaw=seqtab.nochim
mi=match(colnames(seqtab.nochim),seq)
colnames(abundMatRaw)=names(seq)[mi]
abundMatRaw=t(abundMatRaw)
save(abundMatRaw,file=paste0(temp_dir,"/RawAbundanceMat_Euk.Rdata"))

