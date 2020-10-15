# An S4 Helper class to allow NULL values    
setClassUnion("characterOrNULL", c("character","NULL"))    
# An S4 Helper class to allow NULL values    
setClassUnion("numericOrNULL", c("numeric","NULL"))    
# An S4 Helper class to allow NULL values    
setClassUnion("listOrNULL", c("list", "NULL"))    
# An S4 Helper class to allow NULL values    
setClassUnion("dataframeOrNULL", c("data.frame", "NULL"))    
# An S4 Helper class to allow NULL values    
setClassUnion("matrixOrNULL", c("matrix", "NULL"))    

#' A S4 class representing a single image
#'
#' @slot cellprofiler 
hypeRimg <- setClass("hypeRimg", 
		     slots=c(id="character",
			     wd="character",
			     dataFull="dataframeOrNULL",
			     expr="dataframeOrNULL",
			     subExpr="dataframeOrNULL",
			     dimRedData="dataframeOrNULL",
			     dimRedMethod="characterOrNULL",
			     cluster="numericOrNULL",
			     tiffs="listOrNULL"))


#' @title Constructor hypeRimg
#'
#' @description
#' 
#' @return hypeRimg instance
setMethod("initialize", "hypeRimg",
	  function(.Object,
		   id=character) {
		.Object@id <- id
		.Object@wd <- getwd()
		.Object
	  })


#' @title Load cellprofiler data
#' 
#' @description Load cell profiler output for cell segmentation and 
#' intensity features per cell, usually found in cells.txt
#' Only features starting with Intensity_ are retained, for visualization
#' also the features (columns) starting with Location_ are needed.
#'
#' @param file filename 
#'
#' @import data.table
#' @export
loadCellP <- function(obj, file) {
    print("Reading data ... ")
    tmp <- data.frame(fread(file))
    print("Filter data ... ")

    ## TODO: check if necessary information is available:
    # Location and Intensity cols

    obj@dataFull <- tmp
    tmp <- tmp[, which(grepl("Intensity_", colnames(tmp)) & !grepl("Location_", colnames(tmp)))]
    tmp <- data.frame(log(data.matrix(tmp)))
    obj@expr <- tmp

    obj
}


#' @title Load tiffs
#' 
#' @description Load tiff file of any channel, used to visualize clusters.
#'
#' @import tiff
#'
#' @export
loadTiffs <- function(obj, name) {
    if (file.exists(name)) {
	if (!dir.exists(name)) {
	    fl <- name    
	} else {
	    fl <- list.files(name, full.names=T)
	}
	l <- list()
	for (f in fl) {
	    l[[length(l)+1]] <- readTIFF(f)
	}
	names(l) <- fl
	#TODO: append if not null
	obj@tiffs <- l
    } else {
	warning("Folder/file does not exist!")
    }
    obj
}

#' @title Feature selection for final characterization
#'
#' @description Utilizes normalmixE; implemented in the 
#' mixtools package to test for the presence of > 1 population 
#' and identifies minima. Cells are then assigned to high/low 
#' groups. 
#'
#' @param obj hypeRimg object
#' @param k parameter k from mixtools::normalmixEM()
#'
#' @import mixtools
#'
#' @export 
featSel <- function(obj, k=2) {
    tmp <- obj@expr[,]

    coll <- list()
    for (i in 1:length(tmp[1,])) {
	mn <- NA
	tryCatch({
	    val <- tmp[,i]
	    val <- val[which(!is.na(val) & !is.infinite(val))]
	    mn <- normalmixEM(val, k=k)$mu
	}, error=function(e) { warning(e) })
	coll[[i]] <- mn
    }

    df <- unlist(lapply(coll, function(x) x[2]-x[1]))
    wDf <- which(df < quantile(df, 0.1, na.rm=T) | df > quantile(df, 0.9, na.rm=T))
    wDFNM <- colnames(tmp)[wDf]
    #### calculate minima
    collSub <- coll[wDf]
    collMin<- list()
    for (i in 1:length(wDf)) {
        tryCatch({
            dns <- density(tmp[,wDf[i]])
            #plot(dns)            
            df1 <- diff(dns$y)
            collSub[i]
            w <- which(((df1[-1] > 0 & df1[-length(df1)] <0) |  df1[-1] < 0 & df1[-length(df1)] > 0 ))
            w2 <- which((dns$x[w]  >= collSub[[i]][1] & dns$x[w] <= collSub[[i]][2]))# | (dns$x[w] <= collSub[[i]][1] & dns$x[w] >= collSub[[i]][2]))
            xCand <-  dns$x[w[w2]]
            xCand <- xCand[order(dns$y[w[w2]])[1]]
            #abline(v=xCand)      
            collMin[[length(collMin)+1]] <- data.frame(xCut=xCand, i=i, wDf=wDf[i], wDFNM=wDFNM[i])
        }, error=function(e) {print(e)})
    }    
    collMin <- do.call(rbind, collMin)
    collMin <- collMin[which(!is.na(collMin$xCut)),]
    print("#######")
    print(collMin$wDFNM)
    print("---")
    print(wDFNM)
    subExpr <- data.frame(tmp[,which(colnames(tmp) %in% collMin$wDFNM)])
    
    for (i in 1:length(subExpr[1,])) {
        subExpr[,i] <- ifelse(subExpr[,i] >= collMin$xCut[i], 1, 0)
    }
    obj@subExpr <- subExpr

    obj
}

#' @title Dimension reduction
#' 
#' @description Performs dimension reduction on preselected
#' features 
#'
#' @param obj hypeRimg object
#' @param type type, can be umap or tSNE (default). The latter 
#' uses FI-tSNE (https://github.com/KlugerLab/FIt-SNE).
#'
#' @import umap
#'
#' @export
dimRed <- function(obj, type="tSNE") {
    
    if (type == "umap") {
	#TODO: native R implementation
    } else if (type == "umap-learn") {
	tmp <- data.matrix(obj@expr)
	tmp[which(is.na(tmp))] <- 0
	tmp[which(is.infinite(tmp))] <- 0
	um <- umap((tmp), method="umap-learn")
	obj@dimRedData <- data.frame(um$layout)
    } else if (type == "tSNE") {
	if (!exists("fftRtsne")) {
	    stop("Could not find FI-tSNE: fftRtsne(). Did you forget to
		 source the data with source('FIt-SNE/fast_tsne.R', chdir=T)?
		 Check https://github.com/KlugerLab/FIt-SNE for details.")
	}
	tmp <- data.matrix(obj@expr)
	tmp[which(is.na(tmp))] <- 0
	tmp[which(is.infinite(tmp))] <- 0
	ts <- fftRtsne((data.matrix(tmp[,])))
	obj@dimRedData <- data.frame(ts)
    }
    obj
}

#' @title Plot clusters
#'
#' @description Plot dimensionality reduced data with 
#' identified clusters.
#'
#' @param obj hypeRimg instance
#' 
#' @import ggplot2
#' @export
plotClusters <- function(obj) {
    df <- data.frame(obj@dimRedData, CLUSTER=factor(obj@cluster))
    g <- ggplot(df, aes(y=X1, x=X2, color=CLUSTER)) + geom_point() + ggtitle(obj@id)
    print(g)
    return(g)
}

#' @title Identify clusters
#'
#' @description Perform clustering analysis 
#' 
#' @type hcl for hierarchical cluster analysis or km for 
#' kmeans clustering
#'
#' @n number of clusters
#' 
#' @export
identCluster <- function(obj, type="hcl", method="ward.D2", n=5) {
    ### TODO: implement additional methods, automatic #cluster identification
    if (type == "hcl") {
	hc <- hclust(dist(obj@dimRedData), method=method)
	obj@cluster <- cutree(hc, n)
    } else if (type == "km") {
	km <- kmeans(obj@dimRedData, k=k)
	obj@cluster <- km$cluster
    } else {
	warning("Currently only k-means and hierarchical clustering supported")
	obj
    }

    obj
}

#' @title Visualize clusters
#' 
#' @description Show previously identified clusters in a representative 
#' image 
#'
#' @param obj hypeRimg instance
#' @param split can be T (default) or F 
#' @param off size for visualization of the clusters 
#' @param rgb color for background image 
#'
#' @export 
visualize <- function(obj, ...) {
    vis(obj@cluster, obj@dataFull, obj@tiffs[[1]], ...)     
}

vis <- function(cl, data, dat,split=T, off=3, col=rgb(0.2, 0.1, 0.1, 0.2)) {    
    #### tSNE    
    #plot(res, pch=19,cex=0.2, col=cl, axes=F, xlab="", ylab="")    
    
    #### clusters    
    if (!split) {    
        image(log(dat), alpha=0.2, col=col, add=F, axes=F)    
    }    
    empty <- log(dat)    
    for (j in unique(cl)) {    
        empty[] <- NA    
        #image(empty)    
        emp <- as.matrix(empty, ncol=1)    
        lenY <- length(empty[1,])    
        lenX <- length(empty[,1])    
        x <- round(data$Location_Center_X[which(cl == j)])    
        y <- round(data$Location_Center_Y[which(cl == j)])    
        for (k in -off:off) {    
            pos <-(x+k)*lenX+y    
            #pos[which(pos <= 0)] <- 1    
            pos[which(pos <= off)] <- off
            pos[which(pos >= lenX*lenY)] <- lenX*lenY    
            for (i in -off:off) {    
                emp[pos-i] <- j    
            }    
        }    
    
    
        #emp <- as.matrix(emp, ncol=lenX)    
        emp <- matrix(emp, nrow=lenX)    
        if (!split) {    
            image(emp, axes=F, add=T, col=j)    
        } else {    
            image(log(dat), alpha=0.2, col=col, add=F, axes=F)    
            image(emp, axes=F, add=T, col=j)    
        }    
    }    
}


#' @title Assign clusters
#' 
#' @description Quantifies abundance of single channels per cluster
#'
#' @param obj hypeRimg instance
#'
#' @import gridExtra
#'
#' @export 
assignCluster <-function(obj) {
    f <- split(data.frame(obj@subExpr), f=obj@cluster)
    f <- lapply(f, function(x) { 
		     data.frame(POS=apply(x, 2, sum),
				NCLUSTER=length(x[,1]),   
				FRACT_CLUSTER=apply(x, 2, sum)/length(x[,1]))

	})
    for (i in 1:length(f)) { 
	f[[i]] <- data.frame(f[[i]], CLUSTER=names(f)[i], VAR=rownames(f[[i]]), 
			     VAR2=do.call(rbind, strsplit(rownames(f[[i]]), "_"))[,3] ) 
    }

    f1 <-lapply(f, function(x) { 
		    x <- x[rev(order(x$FRACT_CLUSTER)),]
		    x <- x[which(!duplicated(x$VAR2)),]
			     })
    f2 <-lapply(f, function(x) x[rev(order(x$FRACT_CLUSTER)),])
    f2 <- do.call(rbind, f2)
    
    g1 <- ggplot(f2, aes(y=FRACT_CLUSTER, x=CLUSTER, fill=VAR2)) + geom_bar(stat="identity", position="dodge")       
    g2 <- ggplot(f2, aes(y=FRACT_CLUSTER, x=CLUSTER, color=VAR2)) + geom_bar(stat="identity")       

    grid.arrange(g1,g2, nrow=1)

    return(list(f=f, f1=f1, f2=f2))
}
