library(cluster)
#mydist <- function(x) as.dist((1-cor(t(x)))/2)
#mycluster <- function(x, k) list(cluster=cutree(hclust(mydist(x), method = "ward.D2"),k=k))
#pam1 = function(x, k){list(cluster = pam(x,k, cluster.only=TRUE))}

sampleid = c(
  's_C_000012_T001_d','s_C_000167_T001_d','s_C_000170_T001_d','s_C_001395_P001_d','s_C_002747_T001_d','s_C_L3KR5D_M001_d','s_C_001747_P001_d','s_C_001387_P001_d','s_C_NT4RPH_M001_d','s_C_NT4RPH_M002_d','s_C_AE20L6_M001_d'
             ,'s_C_006799_P001_d','s_C_006757_P002_d','s_C_000796_P001_d','s_C_001566_M001_d','s_C_006846_M001_d','s_C_001699_P001_d')
length(sampleid)
for(i in 1:length(sampleid)){
    par(mfrow=c(3,1))
  
    #wes
    set.seed(100)
    data = NULL
    data = filter(wesmaf, Tumor_Sample_Barcode == sampleid[i]) %>% 
              select(ccf_expected_copies_em) %>% filter(ccf_expected_copies_em >=0) %>%
              mutate_all( ~replace(., is.na(.), 0))
    dim(data)
    nmut = dim(data)[1]
    p1 = plot(data, main = sampleid[i], xlab = "CCF")
    
    density_data = NULL
    density_data = density(data$ccf_expected_copies_em)
    p2 = plot(density_data, main = sampleid[i])
    formatted_data = as.data.frame(cbind(density_data$x, density_data$y)); names(formatted_data) = c('CCF','Density')
    head(formatted_data)
    
    gap = NULL #
    gap <- clusGap(formatted_data, kmeans, K.max=ceiling(nmut/15), B=100)
    p3 = plot(gap)
    k = NULL
    k <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method="Tibs2001SEmax") #Tibs2001SEmax #globalSEmax #firstSEmax #firstmax #globalmax
    print(paste(sampleid[i],"nmut=",nmut,"clusters=",k,"...Done!\n")) #2
    p1
    p2
    p3
}