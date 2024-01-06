library(gggenes)
library(ggplot2)
library(export)
library(dplyr)
library(openxlsx)
library(gggenomes)
library(readxl)
library(ggrepel)

# Import manual data
Table_S4_Desulfitobacterium_nos_cluster_annotation <- read_excel("Table S4_Desulfitobacterium nos cluster annotation.xlsx")
output_first_ids <- read_excel("output_first_ids.xlsx")
NosZ_AAI <- read_excel("NosZ_AAI.xlsx")

Table_S4_Desulfitobacterium_nos_cluster_annotation$color<-ifelse(Table_S4_Desulfitobacterium_nos_cluster_annotation$Gene=="nosZ","red","blue")
# Specify the directory containing GFF files
directory_path <- "~/Documents/OneDrive - University of Tennessee/Other microcosm/CladeIII_Nature_template/Fig_1/gff_extract"

# Get a list of all GFF files in the directory
file_names <- list.files(path = directory_path, pattern = "\\.gff$", full.names = TRUE)

# Initialize an empty dataframe
combined_df <- data.frame()

# Loop through each file
for (file_name in file_names) {
  # Read GFF file
  current_df <- read_gff3(file_name)
  # Concatenate dataframes
  combined_df <- bind_rows(combined_df, current_df)
}

#add genome and gene name
combined_df<-as.data.frame(inner_join(combined_df,output_first_ids,by="seq_id"))
combined_df<-as.data.frame(inner_join(combined_df,
Table_S4_Desulfitobacterium_nos_cluster_annotation[,c(1,2,3,4,6)],
by=c("seq_id","start","end","strand")))

combined_df$seq_id<-combined_df$Genome

#add sequence similarity data
#NosZ_AAI<-read_links("NosZAAI.blast",format = "blast")

# write.xlsx(NosZ_AAI,"NosZ_AAI.xlsx") # export to manually adjust

# The aligment auto produced start position is not correct
NosZ_AAI$start<-NosZ_AAI$start1+1 
NosZ_AAI$start2<-NosZ_AAI$start2_1+1 

# add gene name to link data
NosZ_AAI<-as.data.frame(inner_join(NosZ_AAI,
                                 Table_S4_Desulfitobacterium_nos_cluster_annotation[,c(1,2,3,4,6)],
                                 by=c("seq_id","start","end")))

NosZ_AAI<-as.data.frame(inner_join(NosZ_AAI,
                                   Table_S4_Desulfitobacterium_nos_cluster_annotation[,c(6:9)],
                                   by=c("seq_id2","start2","end2")))

link_df<-subset(NosZ_AAI, Gene.x==Gene.y) # subset link data only retain the same gene alignment


# add bin_id1 and bin_id2 to output_first_ids

output_first_ids$seq_id2<-output_first_ids$seq_id

# replace genome as bin_id in link_df
link_df<-as.data.frame(inner_join(link_df,output_first_ids[,c(1,2)],by="seq_id"))
link_df<-as.data.frame(inner_join(link_df,output_first_ids[,c(2,3)],by="seq_id2"))

link_df$seq_id<-link_df$Genome.x
link_df$seq_id2<-link_df$Genome.y



# Gene cluster with orientation as in the gff file
Orientation_cluster<-gggenomes(combined_df,links = link_df)+
  geom_gene(aes(fill=Gene))  +    # draw contig/chromosome lines
  geom_bin_label()+
  geom_link(aes(fill=Gene.x),show.legend = F,offset = 0.3)+
  theme(axis.line.x = element_blank(),axis.ticks = element_blank(),
        axis.text.x = element_blank())

Orientation_cluster%>% flip(3,9,10,13)


# extract the graph position information
combine_df_no_orientation<-Orientation_cluster%>% flip(3,9,10,13)%>% pull_genes()




# change the locus start and end position, make all orientation in forward
combine_df_no_orientation$start<-combine_df_no_orientation$x
combine_df_no_orientation$end<-combine_df_no_orientation$xend
combine_df_no_orientation$strand<-"+"

# make link data no orientation

link_df_no_orientation<-Orientation_cluster%>% flip(3,9,10,13)%>% pull_links()


link_df_no_orientation$start<-if_else(link_df_no_orientation$xend-link_df_no_orientation$x>0,
                                      link_df_no_orientation$x,link_df_no_orientation$xend)
link_df_no_orientation$end<-if_else(link_df_no_orientation$xend-link_df_no_orientation$x>0,
                                      link_df_no_orientation$xend,link_df_no_orientation$x)
link_df_no_orientation$start2<-if_else(link_df_no_orientation$xmax-link_df_no_orientation$xmin>0,
                                      link_df_no_orientation$xmin,link_df_no_orientation$xmax)
link_df_no_orientation$end2<-if_else(link_df_no_orientation$xmax-link_df_no_orientation$xmin>0,
                                    link_df_no_orientation$xmax,link_df_no_orientation$xmin)

link_df_no_orientation$strand<-"+" # remove the cross links demonstrating the strands






# Now, combined_df contains the concatenated data from all GFF files in the specified directory
data = subset(link_df_no_orientation,
              bin_id=="Desulfitobacterium hafniense DCB-2" &
                yend=="10")
data$middle<-(data$x+(data$xend-data$x)/2)

data_Sab = subset(combine_df_no_orientation,
              bin_id=="Candidatus Desulfitobacterium nitrosoreducens Sab" )
data_Sab$Gene.x<-c("hemE","S","sasA","atoC","nosZ","Cyt-b","Cyt-b","HP","MacB","ABC")
data_Sab$middle<-(data_Sab$x+(data_Sab$xend-data_Sab$x)/2)


Desulfito_NosZ<-gggenomes(combine_df_no_orientation,links = link_df_no_orientation)+
  geom_gene(aes(fill=Gene),show.legend = F,
            size=4,shape = 6)+  # draw contig/chromosome line
  geom_bin_label(size = 5)+
  geom_link(aes(alpha=cut(pident,c(100,80,60,40,0))),fill="gray",color="black",offset = 0.3)+
  scale_alpha_manual(
    values = c(0.01, 0.15, 0.3, 0.7),
    name = 'Amino acid identity (%)',
    breaks = c('(0,40]', '(40,60]', '(60,80]', '(80,100]'),
    labels = c('pident < 40', '40 <= pident < 60', '60 <= pident < 80', 'pident >= 80')
  )+
  theme(axis.line.x = element_blank(),axis.ticks = element_blank(),
        axis.text.x = element_blank(),legend.position = c(0.9,0.5))+
  geom_text(data,
             mapping=aes(x=middle,y=0,label=Gene.x),angle=90)+
  geom_text(data_Sab,
            mapping=aes(x=middle,y=15,label=Gene.x),angle=90)


graph2ppt(x=  Desulfito_NosZ,file="Figure",margins=c(0,0,0,0),upscale=T,
          append=T,width=20,height=12)

write.xlsx(output_first_ids,"output_first_ids.xlsx")


save.image(file = "Desulfito_Nos_comparison.Rdata")


