library(dplyr)

#############################helper functions####################################

prepareTranscriptDict <- function(path, allSampleNames){
    transcript_dict = read.table(path)

    for(col in c(5:ncol(transcript_dict))){
        transcript_dict[,col] = gsub('.*transcript:(.*?)\\|.*', '\\1', transcript_dict[,col])
    }
    
    
    names(transcript_dict)[5:(ncol(transcript_dict))] = allSampleNames

    transcript_dict$gffID = sapply(strsplit(transcript_dict$V1,"\\|"), `[`, 1)
    transcript_dict$geneID = sapply(strsplit(transcript_dict$V3,"\\|"), `[`, 1)
        
    transcript_dict = transcript_dict[,c('gffID', 'geneID', allSampleNames)]
        
    return(transcript_dict)
}

findDuplicateMapping <- function(id, df){
    mapping = df[df$dup1==id,]$dup2
    if(mapping %in% df$dup1) {findDuplicateMapping(mapping, df)}
    else return(mapping)
}

buildDuplicateMappingDf <- function(df){
    mappedDF = df
    mappedDF$dup2 = lapply(df$dup1, findDuplicateMapping, df = df)
    return(mappedDF)
}

renameDuplicates <- function(transcript_counts, duplicate_dict_path){
    
    ##read duplicate dateframe
    duplicate_df = read.table(duplicate_dict_path, header = TRUE)

    ##adapt dataframe to directly map to final ID 
    duplicate_df = buildDuplicateMappingDf(duplicate_df)

    ##get indices in count dataframe for duplicated IDs
    indices = match(duplicate_df$dup1,transcript_counts$transcript_id)
    
    ##replace duplicated IDs with mappings
    transcript_counts[indices,]$transcript_id = duplicate_df$dup2
    
    return(transcript_counts)
}

addGffID <- function(transcript_counts, transcript_dict_path, sample, allSampleNames){
    
    ##read transcript mapping dataframe
    transcript_dict = prepareTranscriptDict(transcript_dict_path, allSampleNames)
    
    ##filter gffIDs
    filtered_transcript_counts = transcript_counts[transcript_counts$transcript_id %in% transcript_dict[[sample]],]

    ##match transcriptIDs between count and transcript dict matrix 
    gffIDs_indices = match(filtered_transcript_counts$transcript_id, transcript_dict[!transcript_dict[sample]=='-',][[sample]])
    message('#NAs: ', sum(is.na(gffIDs_indices)))
    gffIDs_indices = na.omit(gffIDs_indices)
    
    ##filter and add gffIDs
    filtered_transcript_counts$gffID = transcript_dict[!transcript_dict[sample]=='-',][gffIDs_indices,]$gffID

    return(filtered_transcript_counts)
}

aggregateDuplicateCounts <- function(transcript_counts){
    
    ##separate into duplciates and non-duplicates
    nonDuplicated_ids = transcript_counts[!(duplicated(transcript_counts$gffID)| duplicated(transcript_counts$gffID, fromLast = TRUE)),]
    Duplicated_ids = transcript_counts[duplicated(transcript_counts$gffID) | duplicated(transcript_counts$gffID, fromLast = TRUE),]
    
    ##just aggregate the duplicates and add transcript/gene ids
    aggregated_counts = aggregate(.~gffID, Duplicated_ids[,-c(1,2)], FUN = sum)
    aggregated_counts[,c("transcript_id", "gene_id")] = Duplicated_ids[match(aggregated_counts$gffID, Duplicated_ids$gffID),][,c("transcript_id", "gene_id")]
    
    ##combine non-duplicates and aggregates duplicates
    transcript_counts = rbind(nonDuplicated_ids,aggregated_counts[,names(nonDuplicated_ids)])
    
    return(transcript_counts)
}

####################################################################################


samples <- read.table('config/samples.tsv', sep =',', header = TRUE)$sample
count_transcript_list = list()


##collaps counts per sample based on gffIDs
for(sample in samples){
    message(sample)
    
    message('reading count matrix')
    expCounts <- paste0("results/FLAMES/sc_long_pipeline/ed", as.character(snakemake@params['edist']),"/", sample ,"/transcript_count.csv.gz")
    tmp_raw_counts<-read.table(file=expCounts,sep=",", header = TRUE)
    
    ##rename duplicates within sample
    message('renaming duplicates')
    dupDf <- paste0("results/gffcompare/ed",as.character(snakemake@params['edist']),"/duplicates/duplicates_", sample,".tsv")
    tmp_raw_counts = renameDuplicates(tmp_raw_counts, dupDf)
    
    ##add gffID for between sample collapse
    message('adding gffID')
    print(snakemake@input['trackingDf'])
    tmp_raw_counts = addGffID(tmp_raw_counts, snakemake@input[['trackingDf']], sample, allSampleNames = samples)
    
    ##collaps counts of duplicated gffIDS
    message('collapsing counts')
    tmp_raw_counts = aggregateDuplicateCounts(tmp_raw_counts)
    
    count_transcript_list[[sample]] = tmp_raw_counts
}


##save collapsed count matrix list as rds 
saveRDS(count_transcript_list, snakemake@output[['collapsedExpressionMatrices']])