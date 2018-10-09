# File: 04_fastQC_trimmed.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the fastq files after trimming
# Date: 9/10/2018


## set variables and source libraries
source('header.R')
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CFastqQuality/master/CFastqQuality.R'
download(url, 'CFastqQuality.R')

# load the required packages
source('CFastqQuality.R')
# delete the file after source
unlink('CFastqQuality.R')

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'Sample')
# another way to get the query, preferred
g_did
dfSample = dbGetQuery(db, "select * from Sample where idData=35;")
# remove any whitespace from the names
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
head(dfSample)
## get the file names from the files table
dbListFields(db, 'File')

## write query to get file names
q = paste0('select * from File where idSample=', dfSample$id, ' AND type="fastq"', ';')

df = lapply(q, function(x) dbGetQuery(db, x))
dfFiles = do.call(rbind, df)
# close connection after getting data
dbDisconnect(db)

#### get the names of the fastq files for first sequencing run
setwd('dataExternal/raw/Trimmed/')
csFiles = list.files('.', pattern = '*.gz')

# sanity check if all files present in directory
dfFiles$name = paste('trim_', dfFiles$name, sep ='')
table(dfFiles$name %in% csFiles)

# # split the file names into batches by list, to run analysis on each sample
# lFilesIndex = split(dfFiles$name, dfFiles$idSample)
# names(lFilesIndex) = dfSample$title

## perform the analysis one sample at a time
## function to write the qa files
write.qa = function(fls, indir, title){
  wd = getwd()
  setwd(indir)
  ob = CFastqQuality(fls, title)
  setwd(wd)
  cat(paste('done', title, '\n'))
  return(ob)
}

ivFilesIndex = seq_along(csFiles)

lOb = lapply(ivFilesIndex, function(x){
  tryCatch(write.qa(dfFiles$name[x], getwd(), as.character(dfFiles$id[x])), error=function(e) NULL)
})

names(lOb) = as.character(dfFiles$id)
setwd(gcswd)
n = make.names(paste('CFastqQuality trimmed S198 run for Alex CRISPR Cas9 sequencing for off target analysis data id 35'))
lOb$meta.1 = dfSample
lOb$meta.2 = dfFiles
lOb$desc = paste('CFastqQuality trimmed S198 run for Alex CRISPR Cas9 sequencing for off target analysis data id 35', date())
n2 = paste0('~/Data/MetaData/', n)
save(lOb, file=n2)

## note: comment out as this entry has been made in db
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/', comment='CFastqQuality S198 run after trimming for Alex CRISPR Cas9 sequencing for off target analysis data id 35')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

### create the plots of interest
getwd()
lOb$desc = NULL
lOb$meta.1 = NULL
lOb$meta.2 = NULL

cvSampleNames = gsub('(Uni\\d+).+', '\\1', dfFiles$name)
cvSampleNames = paste(cvSampleNames, gsub('.+_(R[1|2])_.+', '\\1', dfFiles$name), sep='_')
names(lOb) = cvSampleNames

pdf(file='results/qa.trimmed.fastq.pdf')

iReadCount = sapply(lOb, CFastqQuality.getReadCount)
iReadCount = iReadCount/1e+6

barplot(iReadCount, las=2, main='Trimmed Read Count', ylab = 'No. of Reads in Millions', cex.names =0.8, col=grey.colors(2))

mQuality = sapply(lOb, function(x){
  m = mGetReadQualityByCycle(x)
  m = colMeans(m, na.rm = T)
  return(m)
})

matplot(mQuality, type='l', main='Trimmed base quality', ylab = 'Mean Score', xlab='Position in Read')

lReadWidth = lapply(lOb, iGetReadWidth)
boxplot(lReadWidth, las=2, main='Trimmed Read Width', ylab = 'Read Width', col=grey.colors(2), outline=F, xaxt='n')
axis(1, at=1:length(lReadWidth), labels = names(lReadWidth), cex.axis=0.8, las=2)

## plot all the alphabets by cycle for each forward and reverse reads
i = grep('_R1_', dfFiles$name)

lAlphabets = lapply(i, function(x){
  m = t(mGetAlphabetByCycle(lOb[[x]]))
  m = m[,c('A', 'T', 'G', 'C')]
  r = rowSums(m)
  m = sweep(m, 1, r, '/')
  return(m)
})

matplot(lAlphabets[[1]], type='l', main='Trimmed Sequence Content - Forward Strands', ylab = 'Proportion of Base count', xlab='Position in Read')
temp = lapply(lAlphabets[-1], function(x)
  matlines(x, type='l'))
legend('topleft', legend = colnames(lAlphabets[[1]]), lty=1:4, col=1:4, ncol=2, lwd=2)

# reverse strands
i = grep('_R2_', dfFiles$name)

lAlphabets = lapply(i, function(x){
  m = t(mGetAlphabetByCycle(lOb[[x]]))
  m = m[,c('A', 'T', 'G', 'C')]
  r = rowSums(m)
  m = sweep(m, 1, r, '/')
  return(m)
})

matplot(lAlphabets[[1]], type='l', main='Trimmed Sequence Content - Reverse Strands', ylab = 'Proportion of Base count', xlab='Position in Read')
temp = lapply(lAlphabets[-1], function(x)
  matlines(x, type='l'))
legend('topleft', legend = colnames(lAlphabets[[1]]), lty=1:4, col=1:4, ncol=2, lwd=2)


### some diagnostic plots
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
g_did
q = paste0('select File.id as fid, File.name, File.type, Sample.* from File, Sample where Sample.idData = 35 AND (File.idSample = Sample.id AND File.type like "%fastq%")')
dfBatches = dbGetQuery(db, q)
dbDisconnect(db)

fReadDirection = rep(NA, times=nrow(dfBatches))
i = grepl('_R1_', dfBatches$name)
fReadDirection[i] = 'R1'
fReadDirection[!i] = 'R2'
fReadDirection = factor(fReadDirection)

fBatch = dfBatches$group1
mBatch = mQuality
colnames(mBatch) = paste0(dfBatches$fid, '-', dfBatches$title, '-', as.character(fReadDirection))

oDiag = CDiagnosticPlots(mBatch, 'trimmed Base Quality')
l = CDiagnosticPlotsGetParameters(oDiag)
l$PCA.jitter = F; l$HC.jitter = F;
oDiag = CDiagnosticPlotsSetParameters(oDiag, l)

plot.mean.summary(oDiag, fReadDirection)
plot.sigma.summary(oDiag, fReadDirection)
boxplot.median.summary(oDiag, fReadDirection)
plot.PCA(oDiag, fReadDirection)
plot.dendogram(oDiag, fReadDirection, labels_cex = 0.8)

dev.off(dev.cur())



