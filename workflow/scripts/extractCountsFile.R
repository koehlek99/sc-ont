library(seurat)

count_transcript_list <- readRDS(snakemake@input[['countList']])

samples <- read.table('config/samples.tsv', sep =',', header = TRUE)$sample

##create transcript-level expression objects
seuratList_transcripts = list()
for(sample in samples){
    message(sample)
    tmp_transcript_count = count_transcript_list[[sample]]
    rownames(tmp_transcript_count) = count_transcript_list[[sample]]$gffID
    seuratList_transcripts[[sample]] = CreateSeuratObject(tmp_transcript_count[-c(1,2,ncol(tmp_transcript_count))])
    }
##merge transcript-level expression objects
merged_transcripts_NS = Merge_Seurat_List(seuratList_transcripts)

##add sample origin to seurat object
sampleOrigin <- sapply(strsplit(rownames(merged_transcripts_NS[[]]),"_"), `[`, 2)
samplesOrigin <- samples[as.integer(sampleOrigin)]
merged_transcripts_NS <- AddMetaData(merged_transcripts_NS, samplesOrigin, col.name = 'orig.ident')

##aggregate transcript expression by orig.ident
merged_transcripts_NS_perSample <- AggregateExpression(merged_transcripts_NS, assays = 'RNA', slot = 'counts', group.by = 'orig.ident')
merged_transcripts_NS_perSample <- merged_transcripts_NS_perSample$RNA
rownames(merged_transcripts_NS_perSample) <- lapply(rownames(merged_transcripts_NS_perSample), str_replace, '-', '_')


write.table(merged_transcripts_NS_perSample, snakemake@output[['countsFile']], quote = FALSE, sep = ',')