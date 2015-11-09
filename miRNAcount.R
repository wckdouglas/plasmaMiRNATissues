#!/usr/bin/env Rscript

library(biomaRt)
library(RColorBrewer)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(grid)
library(stringi)

# setup functions
rename_subtype <- function(x){
    paste((unlist(strsplit(x,'-')))[1:3],collapse='-')
}

split_tissues <- function(x){
    stri_replace(x,replacement=' ',fixed='_')
}

remove_hsa <- function(x) {
    paste(unlist(strsplit(x,'-'))[2:3],collapse='-')
}

removeAsterisks <- function(x){
    gsub('[*]','',x)
}

cellType <- function(x){
    x <- as.character(x)
    words <- unlist(strsplit(x,''))
    paste(words[1:(length(words)-1)],collapse='')
}

assignColor <- function(x){
	ifelse(grepl(',',x), 'Orange',
		ifelse(grepl('RBC',x),'red',
		   ifelse(grepl('platelets',x),'green','black')))
}

# read miRNA  database 
filename <- 'D.counts'
filename <- 'MCBZD.counts'
dbpath <- '../tissuesDB'
datapath <- '..//counts'
figurepath <- datapath 
samplename <- unlist(strsplit(filename,'\\.'))[1]
expProf <- dbpath %>%
	stri_c('miRNATissue.csv',sep='/') %>%
    read_csv() %>%
    setNames(make.names(names(.))) %>%
    mutate(Brain = hsa_Cerebellum.adult + hsa_Frontal.cortex.adult + hsa_Midbrain.adult + hsa_Hippocamp.adult) %>%
    mutate(Skin = hsa_Fibrobl.CMV ) %>%
    mutate(Podocytes = hsa_Podocytes.Moins.undiff + hsa_Podocytes.Moins.diff) %>% 
    mutate(Liver = hsa_Liver) %>%
    mutate(Heart = hsa_Heart )%>%
    mutate(Spleen = hsa_Spleen) %>%
    mutate(Peripheral_leukocytes =  hsa_T.cell.CD4 + hsa_T.cell.CD4.2 + hsa_T.cell.CD4.A301 + hsa_T.cell.CD4.ACH2 ) %>%
	mutate(Peripheral_leukocytes=Peripheral_leukocytes + hsa_T.cell.CD4.ACH2.stim + hsa_T.cell.CD4.naive + hsa_T.cell.CD4.effector) %>% 
	mutate(Peripheral_leukocytes = Peripheral_leukocytes + hsa_T.cell.CD4.memory + hsa_T.cell.CD8.2 + hsa_T.cell.CD8 ) %>%
	mutate(Peripheral_leukocytes = Peripheral_leukocytes + hsa_T.cell.CD8.naive + hsa_B.cell.CD19 + hsa_B.cell.CD19.pool ) %>%
	mutate( Peripheral_leukocytes = Peripheral_leukocytes + hsa_B.cell.CD19.2 + hsa_NK.CD56 + hsa_Monocytes.CD14 ) %>%
	mutate( Peripheral_leukocytes = Peripheral_leukocytes +	hsa_Granulocytes.CD15 ) %>%
    mutate(Pituitary_gland = hsa_Pituitary)%>%
    mutate(Thyroid = hsa_Thyroid ) %>%
    mutate(Pancreatic_islets = hsa_Pancreatic.islets)%>%
    mutate(Ovary = hsa_Ovary)%>%
    mutate(Testis = hsa_Testis)%>%
    mutate( Uterus = hsa_Uterus)%>%
    mutate(Placenta = hsa_Placenta)%>%
    mutate(Epididymis = hsa_Epididymis)%>%
    mutate(Prostate = hsa_Prostata ) %>%
    select(matureform,Brain,Skin,Podocytes,Liver,Heart,
            Spleen,Peripheral_leukocytes,Pituitary_gland,
            Thyroid,Pancreatic_islets,Testis,Epididymis)  %>%
    mutate(matureform = sapply(matureform,rename_subtype)) %>%
    group_by(matureform) %>%
    summarise_each(funs(sum)) %>%
	mutate(matureform = stri_replace(matureform,replacement='',fixed='hsa-')) %>%
	mutate(matureform = ifelse(grepl('let',matureform),stri_replace(matureform,replacement='',fixed='-'),matureform)) %>%
	mutate(matureform = ifelse(!grepl('^miR',matureform),paste0('miR',matureform),matureform)) %>%
	mutate(matureform = stri_replace(matureform,replacement='MIR',fixed='miR-')) %>%
	mutate(matureform = toupper(matureform)) %>%
	tbl_df

#read count data
countData <- datapath %>%
	stri_c(filename,sep='/')  %>%
	read_delim(delim='\t',col_names=F,col_types='cnncncccn') %>%
    select(4,9) %>% 
    setNames(c('gene_name','counts')) %>%
    filter(counts>0) %>%
	mutate(gene_name = stri_replace(gene_name,replacement='MIR',fixed='hsa-mir-')) %>%
	mutate(gene_name = ifelse(grepl('[A-Z][0-9]$',gene_name), stri_sub(gene_name,1,nchar(gene_name)-1) ,gene_name)) %>%
	group_by(gene_name) %>%
	summarize(counts = sum(counts)) %>%
	rename(matureform = gene_name) %>%
	inner_join(expProf) %>%
    gather(tissues,expVal,-counts,-matureform) %>%
    arrange(-counts) %>%  
    mutate(rank = 1:n()) %>%
    mutate(matureform = factor(matureform,levels=rev(unique(matureform)))) %>%
    arrange(-expVal) %>%
    mutate(tissues = split_tissues(tissues)) %>%
    group_by(tissues) %>%
    do(data.frame(matureform = .$matureform, 
                expVal = rescale(.$expVal,to=c(1,10)),
                counts = .$counts,
                rank = .$rank))

tissueArrange <- countData %>% 
        group_by(tissues) %>% 
        summarize(summ = sum(expVal)) %>% 
        arrange(-summ) 

countData <- countData %>%  
	ungroup() %>%
    mutate(tissues = factor(tissues,levels=unique(tissueArrange$tissues)))

#incorporate blood cells info
refData <- dbpath %>%
    #stri_c('/bloodComponentsMiRNA.csv') %>%
    stri_c('/journal.pone.0041561.s005.csv') %>%
    read_csv() %>%
	rename(sample=Sample) %>%
    mutate(sample = sapply(sample,remove_hsa)) %>%
	filter(!grepl('-NA',sample)) %>%
    mutate(sample = sapply(sample,removeAsterisks)) %>% 
    gather(cell,count,-sample) %>% 
	replace_na(list(count=0)) %>%
    filter(count > 0 ) %>% 
    mutate(cells = sapply(cell,cellType)) %>%
    group_by(sample,cells) %>%
    summarize(counts = sum(count)) %>%
    ungroup() %>%
    mutate(matureform = stri_replace(sample,replacement = 'MIR',fixed='miR-')) %>%
    group_by(matureform) %>%
    summarize(cell = paste(unique(cells),collapse=','),
			  counts = sum(counts)) %>%
    mutate(matureform = ifelse(!grepl('^MIR',matureform),paste0('MIR',matureform),matureform)) %>%
	mutate(matureform = stri_replace(matureform, fixed='-',replacement='')) %>%
	mutate(matureform = toupper(matureform)) %>%
	arrange(-counts) %>%
	group_by(cell) %>%
	do(data.frame(top_n(.,length(.$counts)/10,counts))) %>%
	select(-counts) %>%
    right_join(countData) %>%
    group_by(matureform) %>%
    do(data.frame(rank = min(.$rank),
                cell = unique(.$cell),
                counts = unique(.$counts))) %>%
	ungroup() %>%
    arrange(rank,-counts) %>%
    mutate(matureform = factor(matureform,levels=unique(countData$matureform))) %>%
    mutate(cell = ifelse(is.na(cell),'none',cell)) %>%
    mutate(color = assignColor(cell)) 
    
colors <- colorRampPalette(c('yellow','navy','blue'))(1000)
p <- ggplot(data=countData,aes(x=as.character(tissues),y=matureform,fill=log10(expVal))) +
    geom_tile() +
    theme(axis.text.x = element_text(angle=75,hjust=1,size=17,color='black'),
          legend.key.width = unit(1.5,'cm'),
          legend.key.height = unit(0.4,'cm'),
          legend.position="bottom",
          legend.text = element_text(size=17),
          legend.title = element_text(size=17),
          axis.text.y = element_text(size = 12 , color = rev(as.character(refData$color))),
          axis.ticks.y = element_blank(),
          panel.background=element_blank())+
    labs(y=' ',x=' ',fill='Normalized Expression')+
    scale_fill_gradientn(colours=colors,
                breaks=c(0,max(log10(countData$expVal))),
                labels=c(0,1),
                guide = guide_colorbar(title.position='top',
                                title.hjust=0,
                                title.vjust=1))
figurename <- stri_c(figurepath ,'/', samplename,'_heatmap_miRNA.pdf')
ggsave(p,file=figurename,width=8,height=16)
print (paste0('saved ',figurename))
        
