# Name: 08_CrisprVariants.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 10/10/2018
# Desc: import the target/off target regions, reference sequences, bam files and call variants


source('header.R')

### load the sample information and bam file names from database
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
dbListFields(db, 'Sample')
# get the query
g_did
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.group3, Sample.title, Sample.description, File.* from Sample, File
           where (Sample.idData = 35) AND (File.idSample = Sample.id AND File.type like "%sort%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)

#### set working directory to appropriate location with bam files
setwd('dataExternal/raw/Aligned/')
csFiles = list.files('.', pattern = '*.bam$')
# check if these files match the file names in database
table(dfSample$name %in% csFiles)

dfSample = dfSample[dfSample$name %in% csFiles, ]
dim(dfSample)
csFiles = dfSample$name

##### load the locations for the target and off target sequences
setwd(gcswd)
df = read.csv('dataExternal/targets.csv', header=T, stringsAsFactors = F, strip.white = T)
str(df)
df$Sequence = gsub(' ', '', df$Sequence)

library(GenomicRanges)
## create genomic ranges objects
f1 = function(x){
  s = strsplit(x$Location, ':')
  s = unlist(strsplit(x$Location, ':'))
  s2 = as.numeric(unlist(strsplit(s[2], '-')))
  s = paste0('chr', s[1])
  gr = GRanges(s, IRanges(s2[1], s2[2]), strand = x$Strand)
  d = data.frame(WGE.ID=x$WGE.ID, seq=as.character(x$Sequence), mismatches=x$Mismatches, type=x$Type)
  elementMetadata(gr) = d
  return(gr)
}

oGRtargets = lapply(1:nrow(df), function(i) f1(df[i,]))
oGRtargets = do.call(c, oGRtargets)

oGRtargets = resize(oGRtargets, width(oGRtargets)+10, fix = 'center')

## load the sequences for these positions from the reference genome
library(BSgenome.Hsapiens.UCSC.hg38)
oDNASSreference = getSeq(BSgenome.Hsapiens.UCSC.hg38, oGRtargets)

## all required data loaded - call variants
library(CrispRVariants)
setwd('dataExternal/raw/Aligned/')

oCrispr = readsToTarget(csFiles, target = oGRtargets[1], reference= oDNASSreference[1], target.loc=22, names=dfSample$title)

lCrispr = lapply(seq_along(oGRtargets), function(x) { readsToTarget(csFiles, target = oGRtargets[x], 
                                                                    reference= oDNASSreference[x], target.loc=22, names=dfSample$title)})

setwd(gcswd)
save(lCrispr, file='temp/lcrispr.rds')

names(lCrispr) = as.character(oGRtargets$WGE.ID)

## drop regions with no alignments
f = sapply(lCrispr, is.null)
which(f)
oGRtargets$bNoAlignments = f

lCrispr = lCrispr[!f]

## extract regions with some samples having no alignments
## this could be due to transfection efficiency
lEfficiency = lapply(lCrispr, mutationEfficiency)

f = sapply(lEfficiency, length)
f = f < 7
table(f)


mEfficiency.Full = do.call(rbind, lEfficiency[!f])
mEfficiency.Partial = do.call(rbind, lEfficiency[f])

# > mEfficiency.Partial
# Uni4 Average Median Overall StDev ReadCount
# 997978255     0       0      0       0    NA         2
# 1061385025    0       0      0       0    NA         6
# 1079010249    0       0      0       0    NA         7
# 1038040742    0       0      0       0    NA         2
# 1159292531    0       0      0       0    NA         2
# 1150575423    0       0      0       0    NA         1
## mark these as regions with no alignments
f = as.character(oGRtargets$WGE.ID) %in% rownames(mEfficiency.Partial)
oGRtargets$bNoAlignments[f] = TRUE

## drop these from the crispr results list
f = sapply(lEfficiency, length)
f = f < 7
table(f)

lCrispr = lCrispr[!f]

f = rep(NA, times=length(oGRtargets))
metadata(oGRtargets) = list(efficiency=mEfficiency.Full, crisprObject=oCrispr)

setwd(gcswd)
n = make.names(paste('oGRanges with Crispr Results for S198 run for Alex data id 35 rds'))
n2 = paste0('~/Data/MetaData/', n)
save(oGRtargets, file=n2)

# comment out as this has been done once
library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
                comment='oGRanges with Crispr Results for S198 run for Alex data id 35 rds')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)
