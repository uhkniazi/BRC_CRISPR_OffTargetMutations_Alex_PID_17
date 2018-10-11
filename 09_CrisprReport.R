# Name: 09_CrisprReport.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 11/10/2018
# Desc: produce a report from the saved results from previous analysis steps

source('header.R')

### load the sample information and bam file names from database
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
# get the query
g_did
q = paste0('select MetaFile.* from MetaFile where idData=', g_did, ';')
df = dbGetQuery(db, q)
dbDisconnect(db)

n = paste0(df$location[4], df$name[4])
load(n)

library(CrispRVariants)
library(GenomicRanges)
library(lattice)
oGRtargets
l = metadata(oGRtargets)
mEff = l$efficiency
oCrispr = l$crisprObject

dfData = data.frame(mEff)
df = stack(dfData[,1:2])
df$names= as.character(rownames(mEff))
str(df)
dotplot(values ~ names, groups=ind,  data=df, auto.key=list(space='top', columns=2, cex=0.8),
        scales=list(x=list(rot=90, cex=0.7)), ylab='Mutation Efficiency')

## extract the location i.e. add another covariate to the data
m = elementMetadata(oGRtargets)
rownames(m) = as.character(m$WGE.ID)
m = m[rownames(dfData), ]

dfData = cbind(dfData, m)
dfData = data.frame(dfData)
df = stack(dfData[, 1:2])
df$names = rownames(dfData)
df$genomic = dfData$type

str(df)
dotplot(values ~ names | genomic, groups=ind,  data=df, auto.key=list(space='top', columns=2, cex=0.8),
        scales=list(x=list(rot=90, cex=0.35)), ylab='Mutation Efficiency')

plotAlignments(oCrispr)
plotVariants(oCrispr)

m = variantCounts(oCrispr)
write.csv(m, file='results/variants.csv')

