# plot species phenologies 

# INPUTS
# - list of phenologies

# OUTPUTS
# - list of phenology plots

library(tidyverse)

# -------------------------------------------------------------------------
# save plots?

save.to.disk <- TRUE

# -------------------------------------------------------------------------

load("data/phenology_list_CBMS_UBMS.Rdata")

n.sp <- length(phenology.list)
pheno.plot.list <- list()
sindex.plot.list <- list() # not used
collated.index.plot.list <- list() # not used

# Phenology plot for TRANSECTS_WALKS:

for(i.sp in 1:n.sp){
  # if(!is.na(phenology.list[[i.sp]])){
  
  # -------------------------------------------------------------------------
  # 1 - species phenologies
  
  pheno <- phenology.list[[i.sp]]
  
  if("trimWEEKNO" %in% names(pheno)){
    pheno.plot <- ggplot(pheno,aes(x = trimWEEKNO, y = NM)) + 
      labs(x = "Monitoring Week", y = "Relative abundance")
  }else{
    pheno.plot <- ggplot(pheno,aes(x = trimDAYNO, y = NM)) + 
      labs(x = "Monitoring Day", y = "Relative abundance")
  }
  
  pheno.plot <- pheno.plot + 
    geom_line(aes(color = M_YEAR)) +
    scale_color_discrete(name = "Year") +
    theme_bw() +
    ggtitle(pheno$SPECIES[1]) +
    NULL
  
  pheno.plot.list[[i.sp]] <- pheno.plot
  
  if(save.to.disk){
    ggsave(paste("results/images/",pheno$SPECIES[1],".tiff",sep=""),pheno.plot,width = 6,height = 4,dpi = 300)
  }
  
}# for every species

names(pheno.plot.list) <- names(as.character(phenology.list))

# -------------------------------------------------------------------------
# store in disk

save(pheno.plot.list,file = "results/images/phenology_plot_list_CBMS_UBMS.Rdata")
