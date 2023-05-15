#!/usr/bin/env Rscript

#loading libraries
library('ggplot2')
library('tidyr')
library("optparse")

#Script to plot the bubble plot for the top 10 genefamilies/pathways/taxa in an experiment.
#The samples can be grouped into single metadata if provided.
#Command line arguments parsing
#Test if there is at least one argument: if not, return an error
#bubble_plot.R --version=0.1

option_list = list(
   make_option(c("-i", "--input"), default = NULL, type ="character", help = "Genefamilies/Pathways/taxa Relative abundance file"),
   make_option(c("-d", "--datatype"), type="character", default = "GF", help = "Input datatype:For GeneFamiles -> GF, Pathways -> PW, Taxa -> T, [default=%default]"),
   make_option(c("-m", "--metadata"), action = "store_false", default = FALSE, help = "If INPUT has metadata in it.[default= %default]"),
   make_option(c("-q", "--qcfile"), type = "character", default = NULL, help = "QC file with SampleName, Total_Reads, Unmapped, Mapped"),
   make_option(c("-o", "--output"), type = "character", default = "bubble_plot", help = "Prefix for output files [defualt=%default]")
 )
opt = parse_args(OptionParser(option_list=option_list))

# Set the datatype

if (toupper(opt$datatype) == "GF"){
   datatype = "Genefamilies" 
} else if (toupper(opt$datatype) == "PW") {
  datatype = "Pathways"
} else if (toupper(opt$datatype) == "T") {
    datatype = "Taxa"
} else {
    datatype = "Genefamilies/Pathways/Taxa"
}


# Processing the QC file
if (!is.null(opt$qcfile)) {
  # Importing the QC file
  qctable = read.csv(opt$qcfile, header = TRUE) # Import QC file
  colnames(qctable)[1]<-"Sample"
  qctable=subset(qctable, select=-c(Total_reads))
  N=nrow(qctable)
  qctable_long = qctable %>%
    gather(key = "labels", value = "n_reads",colnames(qctable)[2]:colnames(qctable)[ncol(qctable)],factor_key = TRUE)
  
  # QC barplot
  qc_plot=ggplot(qctable_long, aes(fill=labels, y=n_reads, x=Sample,width=0.5)) + 
    geom_bar(position="stack", stat="identity")+
    coord_flip()+
    labs(x="Samples", y="# reads",fill="Reads")+
    ggtitle("# mapped vs unmapped reads")+
    theme(
      axis.text.x =element_text(colour = "black", size = 7, face = "bold"),
      axis.text.y = element_text(colour = "black", face = "bold", size = 7),
      legend.title = element_text(size = 7, face = "bold"),
    )
  #Saving the plots
  ggsave(file.path(basename(c("qc_plots.svg"))), plot=qc_plot,width=10,height=0.7*N,device="svg")
}


# Processing abundance file
reldata=read.csv(opt$input, header = TRUE) # Reading abundance file
colnames(reldata)[1]<-"samples" # Renaming the first column
rel_sub=subset(reldata, select=-c(Others)) # Dropping column 'Others'
rel_sub=rel_sub %>% replace(is.na(.), 0)

#Check if the metadata is provided in the input file
if (opt$metadata == "TRUE") {
    
    metadata_col=tail(colnames(rel_sub),n=1) #Metadata data column name
    colnames(rel_sub)[ncol(rel_sub)] <- "metadata"
    
    N=ncol(reldata)-3 #Number of genefamilies/pathways/taxa
    
    #rel_sub$samples=factor(rel_sub$samples, levels = rel_sub$samples)
    #Converting data frame from wide to long.
    reldata_long=rel_sub %>% 
        gather(key="metadata", value="metadata_value",metadata) %>%
        gather(key = "gfpata", value = "abundance",colnames(rel_sub)[2]:colnames(rel_sub)[ncol(rel_sub)-1],factor_key = TRUE)
    reldata_long=reldata_long[order(reldata_long$metadata_value),] #Sort based on metadata

    reldata_long$abundance=round(reldata_long$abundance*100,2) #Rounding-off the abundance to 2 decimal
    
    #write.csv(gf,file="out.csv",quote = FALSE,row.names = FALSE)
    
    #Make plot abundance
    bubble_plot=ggplot(reldata_long, aes(y=samples, x=gfpata)) +
        geom_point(aes(size=abundance,fill=metadata_value),alpha=0.8,shape=21) +
        labs(x="", y="",size="Relative abundance (%)",color="",fill=metadata_col) +
        guides(fill=guide_legend(override.aes = list(size=4)))+ #change the legend size in fill
        ggtitle(paste0("Top ",N," ",datatype,c(" in the experiment"))) +
        scale_size(limits = NULL, range = c(.1,6), breaks = waiver()) +
        theme(legend.key = element_blank(),
              legend.key.size = unit(0.008, 'cm'),
              axis.text.x =element_text(colour = "black", size = 8, angle = 90,hjust=1),
              axis.text.y = element_text(colour = "black", size = 8),
              plot.title = element_text(size = 10, face = "bold"),
              legend.title = element_text(size = 8, face = "bold"),
              panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
              panel.grid.major.y = element_line(colour = "grey95"),
        ) +
        scale_x_discrete(position = "bottom")
} else {
    N=ncol(reldata)-2 #Number of gene families/pathways/taxa
    #Converting data frame from wide to long.
    reldata_long=rel_sub %>% 
        gather(key = "gfpata", value = "abundance",colnames(rel_sub)[2]:colnames(rel_sub)[ncol(rel_sub)],factor_key = TRUE)
    reldata_long$abundance=round(reldata_long$abundance*100,2) #Rounding-off the abundance to 2 decimal
    #write.csv(gf,file="out.csv",quote = FALSE,row.names = FALSE)
    
    #Make plot abundance
    bubble_plot=ggplot(reldata_long, aes(y=samples, x=gfpata)) +
        geom_point(aes(size=abundance,fill=gfpata),alpha=0.8,shape=21) +
        labs(x="", y="",size="Relative abundance (%)",color="",fill=datatype) +
        guides(fill=guide_legend(override.aes = list(size=4)))+ #change the legend size in fill
        ggtitle(paste0("Top ",N," ",datatype, c(" in the experiment"))) +
        scale_size(limits = NULL, range = c(.1,6), breaks = waiver()) +
        theme(legend.key = element_blank(),
              legend.key.size = unit(0.008, 'cm'),
              axis.text.x =element_text(colour = "black", size = 8, angle = 90, hjust=1),
              axis.text.y = element_text(colour = "black", size = 8),
              legend.title = element_text(size = 8, face = "bold"),
              plot.title = element_text(size = 10, face = "bold"),
              panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
              panel.grid.major.y = element_line(colour = "grey95"),
        ) +
        scale_x_discrete(position = "bottom")
}


#ggsave(file.path(paste0(basename(opt$output),c(".svg"))),plot=bubble_plot,width=10,height=0.7*N,device="svg")
ggsave(file.path(paste0(basename(opt$output),c(".svg"))),plot=bubble_plot,device="svg",limitsize=FALSE,units="px",width=750*N,height=200*N)