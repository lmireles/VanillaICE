.xyplot2 <- function(x, data, range, frame=50e3L, ...){
		  dfList <- vector("list", nrow(range))
		  for(i in seq_len(nrow(range))){
			  rm <- findOverlaps(range[i, ], featureData(data), maxgap=frame) ## RangesMatching
			  mm <- matchMatrix(rm)
			  mm.df <- data.frame(mm)
			  mm.df$featureNames <- featureNames(data)[mm.df$subject]
			  marker.index <- mm.df$subject
			  sample.index <- match(sampleNames(range)[i], sampleNames(data))
			  if(any(is.na(sample.index))) stop("sampleNames in RangedData do not match sampleNames in ", class(data), " object")
			  sample.index <- unique(sample.index)
			  data2 <- data[marker.index, sample.index]
			  mm.df$subject <- match(mm.df$featureNames, featureNames(data2))
			  ##
			  ## coersion to data.frame
			  ##
			  df <- as(data2, "data.frame")
			  df$range <- rep(i, nrow(df))##mm.df$query
			  dfList[[i]] <- df
		  }
		  if(length(dfList) == 1) {
			  df <- dfList[[1]]
		  } else{
			  df <- do.call("rbind", dfList)
		  }
		  df$range <- factor(df$range, ordered=TRUE, levels=unique(df$range))
		  df$id <- factor(df$id, ordered=TRUE, levels=unique(df$id))
		  if("return.data.frame" %in% names(list(...))){
			  return.df <- list(...)[["return.data.frame"]]
			  if(return.df) return(df)
		  }
		  list.x <- as.character(x)
		  i <- grep("|", list.x, fixed=TRUE)
		  if(length(i) > 0){
			  zvar <- list.x[[i]]
			  zvar <- strsplit(zvar, " | ", fixed=T)[[1]][[2]]
			  if(zvar == "range"){
				  tmp <- tryCatch(df$range <- mm.df$query, error=function(e) NULL)
			  }
		  }
		  if("gt" %in% colnames(df)){
			  xyplot(x, df,
				 range=range,
				 id=df$id,
				 gt=df$gt,
				 is.snp=df$is.snp,
				 ...)
		  } else {
			  xyplot(x, df,
				 id=df$id,
				 range=range,
				 is.snp=df$is.snp,
				 ...)
		  }
	  }
setMethod("xyplot2", signature(x="formula",
			       data="gSet",
			       range="RangedDataCNV"),
	  function(x, data, range, frame=50e3L, ...){
		  .xyplot2(x=x, data=data, range=range, frame=frame, ...)
	  })

setMethod("xyplot2", signature(x="formula",
			       data="SnpSet",
			       range="RangedDataCNV"),
	  function(x, data, range, frame=50e3L, ...){
		  .xyplot2(x=x, data=data, range=range, frame=frame, ...)
	  })
setMethod("xyplot", signature(x="formula", data="BeadStudioSet"),
	  function(x, data, ...){
		  if("range" %in% names(list(...))){
			  xyplot2(x, data, ...)
		  } else {
			  callNextMethod()
		  }
})

setMethod("xyplot2", signature(x="formula", data="CNSet", range="RangedDataCNV"),
	  function(x, data, range, frame=50e3L, ...){
		  z <- findOverlaps(range, data, maxgap=frame)
		  mm <- matchMatrix(z)
		  mm.df <- data.frame(mm)
		  mm.df$featureNames <- featureNames(data)[mm.df$subject]
		  marker.index <- unique(mm.df$subject)
		  ##marker.index <- featuresInRange(data, rd, FRAME=frame)
		  sample.index <- match(unique(sampleNames(range)), sampleNames(data))
		  data <- data[marker.index, sample.index]
		  ## now we need to know the indices of
		  ## each range after subsetting
		  ## mm.df$subject <- match(mm.df$featureNames, featureNames(data))
		  ## we assume that each range is from a different sample
		  oligoset <- as(data, "oligoSnpSet")
		  df <- as(oligoset, "data.frame")
		  df$range.index <- mm.df$query
		  xyplot(x, df,
			 range=range,
			 gt=df$gt,
			 is.snp=df$is.snp,
			 ...)
	  })

setMethod("xyplot", signature(x="formula", data="SnpSet"),
	  function(x, data, ...){
		  if("range" %in% names(list(...))){
			  xyplot2(x, data, ...)
		  } else {
			  callNextMethod()
		  }
})

setMethod("xyplot", signature(x="formula", data="BeadStudioSet"),
	  function(x, data, ...){
		  if("range" %in% names(list(...))){
			  xyplot2(x, data, ...)
		  } else {
			  callNextMethod()
		  }
})

my.xypanel <- function(x, y,
		       x0, x1, chr.size,
		       col, border, coverage,
		       chr, show.coverage=TRUE,
		       max.y,
		       chromosomeAnnotation,
		       addCentromere=TRUE,
		       ..., subscripts){
	panel.grid(h=-1, v=10)
	panel.xyplot(x, y, ..., subscripts)
	h <- 0.75
	lrect(xleft=x0[subscripts],
	      xright=x1[subscripts],
	      ybottom=y-h/2,
	      ytop=y+h/2,
	      border=border[subscripts],
	      col=col[subscripts], ...)
	if(show.coverage)
		ltext(x, y,labels=coverage[subscripts], cex=0.6)
	##plot centromere
	if(addCentromere){
		chr <- unique(as.integer(as.character(chr)))
		coords <- chromosomeAnnotation[chr, 1:2]/1e6
		lrect(xleft=coords[1],
		      xright=coords[2],
		      ybottom=0,
		      ytop=max.y+h/2,
		      col="grey",
		      border="grey")
	}
}

prepanel.fxn <- function(x,y, chr.size, ..., subscripts){
	list(xlim=c(0, unique(chr.size[subscripts])), ylim=range(as.integer(as.factor(y[subscripts]))))
}

xypanel <- function(x, y,
		    gt,
		    is.snp,
		    range,
		    col.hom="grey20",
		    fill.hom="lightblue",
		    col.het="grey20" ,
		    fill.het="salmon",
		    col.np="grey20",
		    fill.np="grey60",
		    show.state=TRUE,
		    cex.state=1,
		    col.state="blue",
		    cex.pch=0.3,
		    ..., subscripts){
	panel.grid(v=0, h=4, "grey", lty=2)
	panel.xyplot(x[1], y[1], col="white", cex=cex.pch, ...) ## set it up, but don't plot
	is.snp <- is.snp[subscripts]
	if(!missing(gt)){
		gt <- gt[subscripts]
		hets.index <- which(gt == 2)
		hom.index <- which(gt == 1 | gt == 3)
		if(all(!c("col", "fill") %in% names(list(...)))){
			if(any(!is.snp))
				lpoints(x[!is.snp], y[!is.snp], col=col.np,
					fill=fill.np, cex=cex.pch, ...)
			if(length(hom.index) > 0)
				lpoints(x[hom.index], y[hom.index], col=col.hom,
					fill=fill.hom, cex=cex.pch, ...)

			if(length(hets.index) > 0)
				lpoints(x[hets.index], y[hets.index],
					col=col.het,
					fill=fill.het, cex=cex.pch, ...)
		}
	} else {
		lpoints(x[!is.snp], y[!is.snp], col=col.np,
			fill=fill.np, cex=cex.pch, ...)
		## use whatever col.hom to color SNPs
		lpoints(x[is.snp], y[is.snp], col=col.hom,
			fill=fill.hom, cex=cex.pch, ...)
	}
	j <- panel.number()
	st <- start(range)[j]/1e6
	lrect(xleft=st, xright=end(range)[j]/1e6,
	      ybottom=-10, ytop=10, ...)
	if(show.state){
		## left justify the label to the start of the range
		y.max <- current.panel.limits()$ylim[2]
		ltext(st, y.max, labels=paste("state", state(range)[j]),
		      adj=c(0,1), cex=cex.state, col=col.state)
	}
}

arrangeSideBySide <- function(object1, object2){
	grid.newpage()
	lvp <- viewport(x=0,
			y=0.05,
			width=unit(0.50, "npc"),
			height=unit(0.95, "npc"), just=c("left", "bottom"),
			name="lvp")
	pushViewport(lvp)
	nfigs1 <- length(object1$condlevels[[1]])
	nfigs2 <- length(object2$condlevels[[1]])
	stopifnot(length(nfigs1) == length(nfigs2))
	pushViewport(dataViewport(xscale=c(0,1), yscale=c(0.05,1), clip="on"))
	object1$layout <- c(1, nfigs1)
	print(object1, newpage=FALSE, prefix="plot1", more=TRUE)
	upViewport(0)
	lvp2 <- viewport(x=0.5,
			 y=0.05,
			 width=unit(0.50, "npc"),
			 height=unit(0.95, "npc"), just=c("left", "bottom"),
			 name="lvp2")
	pushViewport(lvp2)
	pushViewport(dataViewport(xscale=c(0,1), yscale=c(0.05,1), clip="on"))
	object2$layout <- c(1, nfigs1)
	object2
	print(object2, newpage=FALSE, prefix="plot2", more=TRUE)
}
