# function CK.enrich.v1
# author: Xiaoping Li
# CK.enrich.v2(gff3, cancer, cytoband, ...)

# modifications:
# 4/8/2017: change "chrLen = maxEnd - minStart" to "chrLen = maxEnd - 0"

# 4/8/2017: change chrLen using cytoband chromosome coordinates

CK.enrich.v1 <- function(gff3, cancer, cytoband, ...){
        
        # load libraries
        library(Ckmeans.1d.dp)
        library(biomaRt)
        
        # build data set from gff3 (gencode chromosome dataset), type == "gene"
        # input format modification
        
        genes.gff3 <- gff3[gff3$type == "gene", ]
        genes.gff3$chr <- as.factor(genes.gff3$chr)
        
        
        ###========================= using CKmean to construct zones =========== 
        
        
        
        # to use Ckmeans.1d.dp, we need to get two arguments 
        
        # 1. x: a vector of start and end position for each chromosome
        # 2. k: min k and max k
        
        # max k is calculated as max(20, min(Ni/5, Li/delta))
        
        # to get length (Li) for each chromosome, we need to know min start, and max end in each
        # chromosome and then substract
        
        # Delta (resolution) = 1000000
        
        # (Ni) would be the total number of genes, tapply(gene.gff3$type, gene.gff3$chr, sum)
        
        
        minStart <- tapply(genes.gff3$start, genes.gff3$chr, function(e){
                min = min(e)
               min
        })
        
        maxEnd <- tapply(genes.gff3$end, genes.gff3$chr, function(e){
               max(e)
        })
        
        # length would be maxEnd for each chromosome - inital start position
        # length use information from cytoband
        
        cyto1 <- cytoband[which(cytoband$chr %in% as.character(genes.gff3$chr)), ]
        
        chrLen <- c()
        for(i in 1:length(levels(genes.gff3$chr))){
                
                index <- levels(genes.gff3$chr)[i]
            
                        
                chrLen[i] <- max(c(t(cyto1[cyto1$chr == index, c("start", "end")])))
                
                
        }
        
        
        #total number of genes for each chromosome
        Ni <- tapply(genes.gff3$type, genes.gff3$chr, length)
        
        # array k construct
        
        kmax <- array()
        
        for( i in 1: length(chrLen)){
                delta <- 1000000
                kmax[i] <- max(20, min(Ni[i]/5, chrLen[i]/delta))
                
        }
        
        
        kmin <- rep(2, 25)
        
        # find a way to reform k to c(kmin, kmax), later use
        df <- data.frame(kmin = kmin, kmax = kmax)
        
        
        # use Ckmean.1d.dp function, get zones for each chromosome
        
        #================================   ZONES ==========================================!!!!
        
        zones <- list()
        
        for(i in 1:length(levels(genes.gff3$chr))){
                
                ks <- c(t(df[i, c("kmin", "kmax")])) # this step change the format of k into a vector
                
                
                # iterate through chromosomes
                gene.cord <- c(t(genes.gff3[genes.gff3$chr == (levels(genes.gff3$chr)[i]), c("start", "end")]))
                zones[[i]] <- Ckmeans.1d.dp(gene.cord, ks)
                names(zones)[i] <- levels(genes.gff3$chr)[i]
        }
        
        #============================= ZONE BOUNDARIES ==========================================!!!!
        
        
        # the boundaries are defined as the mid points between each center
        
        boundaries <- vector( mode = "list", length(zones))
        for(i in 1:length(levels(genes.gff3$chr))){
                
                centers <- zones[[i]]$centers
                
                
                for( k in 1:length(centers)){
                        
                        if(k > 1){
                                
                                boundaries[[i]][k] <- mean(c(centers[k], centers[k-1]))
                                
                        }else{
                                
                                next
                                
                        }
                        
                }
                
                boundaries[[i]][1] <- minStart[i]    # do not include these two in the loop
                boundaries[[i]][length(centers) + 1] <- maxEnd[i]
                
                #boundaries[[i]][1] <- 0
                #boundaries[[i]][length(centers + 1)] <- chrLen[i]
                
                
                names(boundaries)[i] <- names(zones)[i]
                
                
                
        }
        
        
        
        ##================================ cancer dataset modification ================
        
        # cancer$chr contains chr "-", we need to dig more information about chromosome
        
        # use entrez identifier to retrieve the chromosome number
        
        # build mart
        mart <- useMart("ensembl", dataset =  "hsapiens_gene_ensembl")
        
        simpleCancer <- cancer[, c("chromosome", "tstart", "tend", "entrez", "symbol")]
        
        for( i in 1: nrow(simpleCancer)){
                
                chr <- simpleCancer$chromosome
                
                if( chr[i] == "-" | is.na(chr[i] )){
                        
                        values <- simpleCancer$entrez[i]
                        
                        g <- getBM(c("chromosome_name", "start_position", "end_position"), filters = "entrezgene", values = values, mart) # retrieve from the data base
                        
                        if(is.na(g[1,1])){
                                
                                
                                simpleCancer$chromosome[i] <- g[1, 1]
                                simpleCancer$tstart[i] <- g[1, 2]
                                simpleCancer$tend[i] <- g[1, 3]
                                
                        } else if(is.numeric(g[1,1])){
                                
                                simpleCancer$chromosome[i] <- paste0("chr", g[1, 1])
                                simpleCancer$tstart[i] <- g[1, 2]
                                simpleCancer$tend[i] <- g[1, 3]
                                
                        } else if(g[1,1] == "X"){
                                
                                simpleCancer$chromosome[i] <- paste0("chr", g[1, 1])
                                simpleCancer$tstart[i] <- g[1,2]
                                simpleCancer$tend[i] <- g[1,3]
                                
                        } else if(g[1,1] == "MT"){
                                
                                simpleCancer$chromosome[i] <- "chrM"
                                simpleCancer$tstart[i] <- g[1, 2]
                                simpleCancer$tend[i] <- g[1, 3]
                                
                        } else if(nchar(g[1,1] > 4)){
                                
                                pattern <- "(CHR_HS)(CHR[0-9]{1,2})(_.*)"
                                
                                
                                m <- regexec(pattern, g[1,1])
                                
                                v <- regmatches(g[1,1], m)
                                
                                simpleCancer$chromosome[i] <- tolower(v[[1]][3])
                                simpleCancer$tstart[i] <- g[1, 2]
                                simpleCancer$tend[i] <- g[1, 3]
                                
                        } 
                        
                        
                } else {
                        
                        next
                        
                }
                
                
        }
        
        
        # reduce the dataset
        cancer_genes <- simpleCancer[, c("chromosome", "tstart", "tend", "symbol")]
        
        names(cancer_genes) <- c("chr", "start", "end", "symbol")
        
        cancer_genes$chr <- as.factor(cancer_genes$chr)
        
        # remove NA [[[[[[[[[[[[[[[[[[[[[[[[]]]]]]]]]]]]]]]]]]]]]]]]
        cancer_genes <- cancer_genes[complete.cases(cancer_genes$chr),]
        
        #### ============================ Cytoband datasset ==================================
        #cytoband
        cyto <- cytoband[cytoband$chr %in% unique(cancer_genes$chr), ]
        
        
        ## create a data.frame based on cancer_genes dataset, adding all other information
        
        
        
        zone_number <- c()  # the zone in chromosome (gff) for each cancer gene 
        total_zones <- c()  # the total number of zones on certain chromosome
        ngenes_zone <- c()  # number of genes in that perticular zone of a particular chr
        cyto_bands_1 <- c() # cytobands determined by start position of a zone that a cancer gene falls into
        cyto_bands_2 <- c() # cytobands determined by the end position of a zone that above cancer gene falls inot
        cyto_bands <- c() # integrate cyto_bands_1 and cyto_bands_2
        zone_start_chr <- c() # a particular zone start position on a chromosome
        zone_end_chr <- c() # a particular zone end position on a chromosome
        
        for( i in 1: nrow(cancer_genes)){
                
                start_p <- as.numeric(cancer_genes$start[i])
                
                end_p <- as.numeric(cancer_genes$end[i])
                
                chr_index <- cancer_genes$chr[i]
                
                chr_boundary <- boundaries[[chr_index]]
                
                if(start_p >= min(chr_boundary) & end_p <= max(chr_boundary)){
                        
                        mid <- (start_p + end_p)/2 #use mid point to determing which zone this gene is belonging to if it crosses two zones
                        zone_number[i] <- max(which(mid > chr_boundary))
                        ngenes_zone[i] <- round((0.5*zones[[chr_index]]$size)[zone_number[i]], 0)
                        zone_start_chr[i] <- chr_boundary[zone_number[i]]
                        zone_end_chr[i] <- chr_boundary[zone_number[i] + 1]
                        
                } else {
                        
                        zone_number[i] <- NA
                        ngenes_zone[i] <- NA
                        zone_start_chr[i] <- NA
                        zone_end_chr[i] <- NA
                        
                }
                
                
                cyto_start <- cyto$start[cyto$chr == chr_index]
                cyto_end <- cyto$end[cyto$chr == chr_index]
                cyto_name <- cyto$name[cyto$chr == chr_index]
                
                if(is.na(zone_start_chr[i]) & is.na(zone_end_chr[i])){
                        
                        cyto_bands_1[i] <- NA
                        cyto_bands_2[i] <- NA
                        
                } else if(zone_start_chr[i] >= min(cyto_start) & zone_end_chr[i] <= max(cyto_end)){
                        
                        cyto_bands_1[i] <- cyto$name[max(which(zone_start_chr[i] >= cyto_start))]
                        cyto_bands_2[i] <- cyto$name[min(which(zone_end_chr[i] <= cyto_end))]
                        
                        
                } else {
                        
                        cyto_bands_1[i] <- NA
                        cyto_bands_2[i] <- NA
                }
                
                
                # need to integrate cytobands
                if(identical(cyto_bands_1[i], cyto_bands_2[i])){
                        
                        cyto_bands[i] <- cyto_bands_1[i]
                        
                } else {
                        
                        cyto_bands[i] <- paste0(cyto_bands_1[i], "-", cyto_bands_2[i])
                }
                
                
                
                total_zones[i] <- length(chr_boundary) - 1
                
        } 
        
        
    
        
        # 3 genes from chr10 out of boundary !!!!!!!!!!!!
        
        # integrate all the elements
        
        mod_cancer <- cbind(cancer_genes, zone_number, zone_start_chr, zone_end_chr, total_zones, ngenes_zone, cyto_bands)
        
        
        #remove NAs [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
        #mod_cancer <- mod_cancer[which(complete.cases(mod_cancer)),]
        
        
        #============= to get total number of zones for p.adjust ( x number of hypothesis)
        
        p_length <- lapply(zones, function(x)length(x$centers)) # bonferroni alpha/p_length for each chromosome
        
        
        
        
        # to get number of cancer genes in the zone
        
        ###====================== geting probability for each zone ==================
        
        #hypergeometric distribution
        
        
        # N is the total number of genes in the human
        # M is the total number of cancer genes
        # ngenes is the number of genes fall inot a given genomic zone
        # m is the number of the cancer gene fall into a given genomic zone
        # probability p(m|n, N, M) = Pr(X = m | n, N, M) = (M choose m) * ( N - M choose n - m) / ( N choose n)
        
        
        
        mgenes_zone <- c()
        probs <- c()
        p.adjusted <- c()
        
        
        
        for( i in 1 : nrow(mod_cancer)){
                
                M <- nrow(cancer)
                N <- nrow(genes.gff3)
                
                
                mgenes_zone[i] <- length(which(mod_cancer$chr == mod_cancer$chr[i] & mod_cancer$zone_number == mod_cancer$zone_number[i]))
                
                
                probs[i] <- phyper(mgenes_zone[i], M, N - M, mod_cancer$ngenes_zone[i], lower.tail = FALSE)
                
                p.adjusted[i] <- p.adjust(probs[i], method = "bonferroni", n = p_length[[mod_cancer$chr[i]]])
                
        }
        
        
        final_cancer <- cbind(mod_cancer, mgenes_zone, probs, p.adjusted)
        
        
        significant <- final_cancer[which(final_cancer$p.adjusted <= 0.05), ]
        
        return(significant)
      
}


