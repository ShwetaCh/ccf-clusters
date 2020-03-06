library(cluster)
library(plyr)
library(dplyr)
library(stringr)
library(data.table)

#mydist <- function(x) as.dist((1-cor(t(x)))/2)
#mycluster <- function(x, k) list(cluster=cutree(hclust(mydist(x), method = "ward.D2"),k=k))
#pam1 = function(x, k){list(cluster = pam(x,k, cluster.only=TRUE))}

############### DECLARE ALL FILE NAMES ###############
wesmaf = "~/wes_mafanno_ccf.maf"
somePDFPath = "~/wes_ccf_clusters_030520.pdf"
sometextoutput = "~/wes_ccf_clusters_030520.txt"
######################################################

######################################################
maf = wesmaf
sampleid = sort(unique(wesmaf$Tumor_Sample_Barcode))
head(sampleid)
N = length(unique(sampleid))
pdf(file=somePDFPath)  
allplots = list() #*
######################################################

for(i in 1:N)
  local({
  par(mfrow=c(3,1))
  set.seed(100)
  p1 = NULL; p2 = NULL; p3 = NULL;
  #wes
  data = NULL
  data = filter(maf, Tumor_Sample_Barcode == sampleid[i]) %>% 
    #filter(Variant_Classification %in% exonic_variant_class) %>% 
    select(ccf_expected_copies_em) %>% filter(ccf_expected_copies_em >=0) %>%
    mutate_all( ~replace(., is.na(.), 0))
  dim(data)
  nmut = dim(data)[1]
  
  if(nmut >=10){
    
            p1 = plot(data, main = sampleid[i], xlab = "CCF")
            
            density_data = NULL
            density_data = density(data$ccf_expected_copies_em)
            p2 = plot(density_data, main = sampleid[i])
            
            formatted_data = NULL; gap = NULL; k = NULL; cl =NULL; 
            formatted_data = as.data.frame(cbind(density_data$x, density_data$y)); names(formatted_data) = c('CCF','Density')
            head(formatted_data)
            
            # if(nmut<65){
            #   kmax = floor(nmut/15); 
            #   if(kmax <= 1){ kmax = 2} #Errors out when K.max is 1, can't do clustering anymore
            #   } else{
            #     kmax = ceiling(nmut/10)  
            #   }
            
            gap = clusGap(formatted_data, kmeans, d.power = 2, K.max = kmax, B=100)
            k = maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method="Tibs2001SEmax") #Tibs2001SEmax #globalSEmax #firstSEmax
            #kmeans(formatted_data, centers=5, iter.max = 100, trace=FALSE)
            
            # if(nmut >= 65){
            #   cl = k+1} else{
            #     cl = k}
            
            cl = k+1
            
            p3 = plot(gap, main = paste(sampleid[i],"PREDICTED CLUSTERS=",cl))
            p1
            p2
            p3
            
            p = cowplot::plot_grid(plotlist = list(p1,p2,p3), nrow = 3) #*any()
            allplots[[i]] <<- p #*
            print(p) #*
            print(paste(i,"-->",sampleid[i],"nmut=",nmut,"clusters=",cl,"...Done!"))
            df = paste(sampleid[i],nmut,cl,sep="\t")
            write.table(df,sometextoutput, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
    } else {
            df = paste(sampleid[i],nmut,"skip-NA",sep="\t")
            write.table(df,sometextoutput, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE);
           }
})
dev.off()
######################################################
######################################################
