.sss.valid <- function(x){
 el <- list(data="data.frame", coords="data.frame", grid="matrix",
            knots="list", W="matrix", contract="logical", regular="logical")
 requireNamespace("methods",quietly=TRUE)
 exists.slot <- match(names(el),methods::slotNames(x))
 if(!any(is.na(exists.slot))){
    if(!is.null(x@data)){
       if(nrow(x@data)==nrow(x@coords)){
          if(!(x@contract && nrow(x@data)==nrow(x@W) && nrow(x@grid)==ncol(x@W)) || !(!x@contract && nrow(x@data)==ncol(x@W) && nrow(x@grid)==nrow(x@W))){
             paste("contract logical does not match with W matrix.")
          }
       } else {
          paste("Note! nrow(data) != nrow(coords).")
       }
    } else {
       if(!(x@contract && nrow(x@grid)==ncol(x@W)) || !(!x@contract && nrow(x@grid)==nrow(x@W))){
          paste("contract logical does not match with W matrix.")
       }
    }
    if(ncol(x@coords)==2L && ncol(x@grid)==2L && length(x@knots)==2L){
       valid = NULL
       for(i in 1:length(el)){
           requireNamespace("methods",quietly=TRUE)
           valid[i] <- any(grepl(paste("^",el[[i]],"$",sep=""),class(methods::slot(x,names(el)[[i]]))))
       }
     if(all(valid)) TRUE else paste("Slots",paste(names(el)[!valid],sep="",collapse=", "),"not match their class.")
    } else {
     paste("Number of columns in coords/grid, or length in knots is not 2.")
    }
 } else {
  paste("Slots",paste(names(el)[is.na(exists.slot)],sep="",collapse=", "),"are missing.")
 }
}
sss <- setClass("sss",
	slots = list(data="data.frame", coords="data.frame", grid="matrix",
				knots="list", W="matrix", contract="logical", regular="logical"),
	package = "scpm", validity = .sss.valid)
.createGrid <- function(coords, check.regular = TRUE){
	if (!is.matrix(coords))
		coords <- as.matrix(coords)
	if (length(dim(coords)) > 2L || !is.numeric(coords) || is.complex(coords)) 
		stop("cooords must be a numeric, real matrix")
	dnx = dimnames(coords)
	if(is.null(dnx)) dnx <- vector("list", 2)
	o.c = order(coords[,2],coords[,1])
	coords = coords[o.c,]
	x.grid = as.numeric(.getKnots(coords[,1]))
	y.grid = as.numeric(.getKnots(coords[,2]))
	xy.grid = list(abcissa = x.grid, ordinate= y.grid)
	names(xy.grid) = dnx[[2]]
	y.x = matrix(as.numeric(unlist(strsplit(unlist(lapply(y.grid,paste,x.grid,sep=" "))," "))),ncol=2,byrow=TRUE)
	grid.s = cbind(y.x[,2],y.x[,1])
	colnames(grid.s) = dnx[[2]]
	attr.gs = attributes(grid.s)
	attr.gs$grid.list = xy.grid
	attr.gs$W = .incidenceSpatial(coords,grid.s)
	attr.gs$ordering = o.c
	attr.gs$regular = nrow(coords)==nrow(grid.s) && all(apply(coords,1,paste,sep="",collapse=":")==apply(grid.s,1,paste,sep="",collapse=":"))
	attributes(grid.s) = attr.gs
	return(grid.s)
}
.fixLoc <- function(data,smallest=0.00000010,largest=0.00000125,digits=8){
    if(class(data)=="sss"){
        coords <- data@coords
    } else {
        coords <- data
    }
    nc <- colnames(coords)
    locsId <- paste(coords[,1],coords[,2],sep=",")
    requireNamespace("stats",quietly=TRUE)
    locs <- stats::xtabs(~locsId)
    locx <- coords[,1]
    locy <- coords[,2]
    for(i in 1:length(locs)){
        isrep <- NULL
        for(j in 1:nrow(coords)){
            if(locsId[j]==names(locs)[i]){
                isrep <- c(isrep,i)
                if(length(isrep) < locs[i]){
                    usign <- runif(2,0,1)
                    locx[j] <- coords[j,1] + ifelse(usign[1]>0.5,1,-1)*round(runif(1,smallest,largest),digits)
                    locy[j] <- coords[j,2] + ifelse(usign[2]>0.5,1,-1)*round(runif(1,smallest,largest),digits)
                }
            }
        }
    }
    locsNewId <- paste(locx,locy,sep=",")
    coordsNew <- data.frame(locx,locy)
    names(coordsNew) <- nc
    if(class(data)=="sss"){
        dataOld <- data.frame(data@data,coords)
    	data <- .asNewsss(dataOld, coords = coordsNew, coords.col = NULL, data.col = NULL, contract = TRUE, ordering = "grid")
    } else {
        dataOld <- data.frame(coords)
        data <- data.frame(dataOld,coordsNew)
    }
    return(data)
}
"create.sss" <- function(coords, data, ...){
	if(missing(coords)) stop("A two-columns matrix or data-frame of coordinates must be specified")
	if(missing(data)) data <- data.frame() else data <- as.data.frame(data)
	if(is.data.frame(coords))
		coords <- as.matrix(coords)
	if(ncol(coords)!=2)
		stop("coords must be a two-column matrix or data-frame")
	if(nrow(data)>0)
		if(nrow(data)!=nrow(coords)) stop("coords and data must have the same number of rows")
	X <- cbind(coords,data)
	coords.col <- 1:2
	data.col <- NULL
	if(ncol(data)>0)
		data.col <- c(1:ncol(data)) + 2
	contract <- TRUE
	ordering <- "grid"
	if(missing(X) || is.null(X)) stop("X must be specified!")
	if(!is.data.frame(X) && !is.matrix(X) && !is.list(X)) stop("X must be a data.frame, matrix or list of elements.")
	if(is.matrix(X)) X <- as.data.frame(X)
    requireNamespace("methods",quietly=TRUE)
	if(is.data.frame(X)){
		attr.X = attributes(X)
		if(is.null(coords)) coords <- attr.X$coords
		if(is.null(coords) && (is.null(coords.col) || !is.numeric(coords.col))) stop("If X is a data.frame/matrix without a coords attribute then coords.col numeric vector must be specified.")
		if(!is.null(coords.col)) coords <- X[,coords.col]
		if(missing(data.col) || is.null(data.col)) data.col <- 1:ncol(X)
		data.X <- X
	} else {
		if((is.null(coords) && is.null(X$coords) && !methods::.hasSlot(X,"coords")) && dim(coords)[2]!=2) stop("If X is a list without a coords elements/slots then the coords n*2 numeric matrix must be specified")
		if(is.null(coords) && !is.null(X$coords) && !methods::.hasSlot(X,"coords")) coords <- X$coords
		if(is.null(coords) && is.null(X$coords) && !methods::.hasSlot(X,"coords")) coords <- X@coords
		if(is.null(X$data) && !!methods::.hasSlot(X,"data")) data.X <- NULL;dnx <- NULL
		if(!is.null(X$data) && !!methods::.hasSlot(X,"data")) data.X <- as.data.frame(X$data)
		if(is.null(X$data) && !methods::.hasSlot(X,"data")) data.X <- as.data.frame(X@data)
		if(!is.null(coords.col)) coords <- as.matrix(coords)[,coords.col]
		if(!exists("data.X")) data.X <- NULL
		if(!is.null(data.X) && (missing(data.col) || is.null(data.col))) data.col <- 1:ncol(data.X)
		if(!exists("data.col")) data.col <- NULL
		attr.X = attributes(data.X)
	}
	if(!is.null(data.col)){
		data.X <- data.frame(data.X[,data.col],row.names=attr.X$row.names)
		names(data.X) <- attr.X$names[data.col]
	}
	i.o = apply(cbind(ifelse(rep(is.null(data.X),nrow(coords)),1,data.X),coords),1,paste,sep="",collapse=":")
	grid.X = .createGrid(coords, contract)
	attr.new = attributes(grid.X)
	if(!is.null(data.X) && !is.null(data.col)){
		data.X <- data.frame(data.X[attr.new$ordering,],row.names=attr.X$row.names[attr.new$ordering])
		names(data.X) <- attr.X$names[data.col]
	} else {
		data.X <- data.frame(data.X[attr.new$ordering,],row.names=attr.X$row.names[attr.new$ordering])
		names(data.X) <- attr.X$names
	}
	coords.X = coords[attr.new$ordering,]
	grid.data = grid.X
	if(ordering=="original"){
		f.o = apply(cbind(ifelse(rep(is.null(data.X),nrow(coords)),1,data.X),coords.X),1,paste,sep="",collapse=":")
		get.back = match(i.o,f.o)
		if(!is.null(data.X) && !is.null(data.col)){
			data.X <- data.frame(data.X[get.back,],row.names=attr.X$row.names)
			names(data.X) <- attr.X$names[data.col]
		} else {
			data.X <- data.frame(data.X[get.back,],row.names=attr.X$row.names)
			names(data.X) <- attr.X$names
		}
		coords.X = coords.X[get.back,]
		if(attr.new$regular){
			grid.data = grid.data[get.back,]
		}
		if(contract) attr.new$W = attr.new$W[get.back,] else attr.new$W = attr.new$W[,get.back]
	}
	attributes(grid.data) <- attr.new[1:3]
	newdata <- methods::new("sss", data = data.frame(data.X), coords = data.frame(coords.X),
            grid = as.matrix(grid.data), knots = attr.new$grid.list,
            W = attr.new$W, contract = contract, regular = attr.new$regular)
	newdata <- .fixLoc(newdata, ...)
	return(newdata)
}
is.sss <- function(x){
    requireNamespace("methods",quietly=TRUE)
	methods::is(x,"sss")
}
.asNewsss <- function(X, coords = NULL, coords.col = NULL, data.col = NULL, contract = TRUE, ordering = c("original","grid"), ...){
	if(missing(X) || is.null(X)) stop("X must be specified!")
	.defaultArg("ordering","grid")
	if(!is.data.frame(X) && !is.matrix(X) && !is.list(X)) stop("X must be a data.frame, matrix or list of elements.")
	if(is.matrix(X)) X <- as.data.frame(X)
    requireNamespace("methods",quietly=TRUE)
	if(is.data.frame(X)){
		attr.X = attributes(X)
		if(is.null(coords)) coords <- attr.X$coords
		if(is.null(coords) && (is.null(coords.col) || !is.numeric(coords.col))) stop("If X is a data.frame/matrix without a coords attribute then coords.col numeric vector must be specified.")
		if(!is.null(coords.col)) coords <- X[,coords.col]
		if(missing(data.col) || is.null(data.col)) data.col <- 1:ncol(X)
		data.X <- X
	} else {
		if((is.null(coords) && is.null(X$coords) && !methods::.hasSlot(X,"coords")) && dim(coords)[2]!=2) stop("If X is a list without a coords elements/slots then the coords n*2 numeric matrix must be specified")
		if(is.null(coords) && !is.null(X$coords) && !methods::.hasSlot(X,"coords")) coords <- X$coords
		if(is.null(coords) && is.null(X$coords) && methods::.hasSlot(X,"coords")) coords <- X@coords
		if(is.null(X$data) && !methods::.hasSlot(X,"data")) data.X <- NULL;dnx <- NULL
		if(!is.null(X$data) && !methods::.hasSlot(X,"data")) data.X <- as.data.frame(X$data)
		if(is.null(X$data) && methods::.hasSlot(X,"data")) data.X <- as.data.frame(X@data)
		if(!is.null(coords.col)) coords <- as.matrix(coords)[,coords.col]
		if(!exists("data.X")) data.X <- NULL
		if(!is.null(data.X) && (missing(data.col) || is.null(data.col))) data.col <- 1:ncol(data.X)
		if(!exists("data.col")) data.col <- NULL
		attr.X = attributes(data.X)
	}
	if(!is.null(data.col)){
		data.X <- data.frame(data.X[,data.col],row.names=attr.X$row.names)
		names(data.X) <- attr.X$names[data.col]
	}
	i.o = apply(cbind(ifelse(rep(is.null(data.X),nrow(coords)),1,data.X),coords),1,paste,sep="",collapse=":")
	grid.X = .createGrid(coords, contract)
	attr.new = attributes(grid.X)
	if(!is.null(data.X) && !is.null(data.col)){
		data.X <- data.frame(data.X[attr.new$ordering,],row.names=attr.X$row.names[attr.new$ordering])
		names(data.X) <- attr.X$names[data.col]
	} else {
		data.X <- data.frame(data.X[attr.new$ordering,],row.names=attr.X$row.names[attr.new$ordering])
		names(data.X) <- attr.X$names
	}
	coords.X = coords[attr.new$ordering,]
	grid.data = grid.X
	if(ordering=="original"){
		f.o = apply(cbind(ifelse(rep(is.null(data.X),nrow(coords)),1,data.X),coords.X),1,paste,sep="",collapse=":")
		get.back = match(i.o,f.o)
		if(!is.null(data.X) && !is.null(data.col)){
			data.X <- data.frame(data.X[get.back,],row.names=attr.X$row.names)
			names(data.X) <- attr.X$names[data.col]
		} else {
			data.X <- data.frame(data.X[get.back,],row.names=attr.X$row.names)
			names(data.X) <- attr.X$names
		}
		coords.X = coords.X[get.back,]
		if(attr.new$regular){
			grid.data = grid.data[get.back,]
		}
		if(contract) attr.new$W = attr.new$W[get.back,] else attr.new$W = attr.new$W[,get.back]
	}
	attributes(grid.data) <- attr.new[1:3]
	newdata <- methods::new("sss", data = data.frame(data.X), coords = data.frame(coords.X),
		grid = as.matrix(grid.data), knots = attr.new$grid.list,
		W = attr.new$W, contract = contract, regular = attr.new$regular)
	return(newdata)
}
"as.sss" <- function(X, coords, coords.col, data.col, ...){
	if(missing(coords)){coords <- NULL}
	if(missing(coords.col)){coords.col <- NULL}
	if(missing(data.col)){data.col <- NULL}
	newdata <- .asNewsss(X = X, coords = coords, coords.col = coords.col, data.col = data.col, contract = TRUE, ordering = "grid")
	newdata <- .fixLoc(newdata, ...)
	return(newdata)
}
"sss2df" <- function(x){
    if(class(x)=="sss"){
        datum <- data.frame(x@coords,x@data)
        attr2df <- attributes(x)
        attr2df$class <- "data.frame"
        attr2df$contract <- NULL
        attr2df$type <- NULL
        attr(datum,"data") <- attr2df$data
        attr(datum,"coords") <- attr2df$coords
        attr(datum,"grid") <- attr2df$grid
        attr(datum,"knots") <- attr2df$knots
        attr(datum,"W") <- attr2df$W
        attr(datum,"regular") <- attr2df$regular
        return(datum)
    } else {
        stop("Not an element from class sss")
    }
}
linear <- function(formula, data = NULL, contrasts = NULL, intercept = FALSE){
	requireNamespace("stats", quietly = TRUE)
    intercept <- any(grepl("(^1|+1)",formula))
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "contrasts"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    o = match(c("formula"),names(mf))
    names(mf)[o] = "object"
    mf[[1L]] <- quote(stats::model.matrix)
    mf <- eval(mf, parent.frame())
    if(!intercept){
    	if(!is.na(match("(Intercept)",colnames(mf)))){
			nn <- colnames(mf)
    		mf = as.matrix(mf[,-match("(Intercept)",nn)])
			colnames(mf) <- nn[-match("(Intercept)",nn)]
    	}
    }
	return(mf)
}
cp <- function(x, psi, data = NULL, groups = NULL, contrasts = NULL, only.UV = FALSE){
	requireNamespace("stats", quietly = TRUE)
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("x", "data"), names(mf), 0L)
	c0 <- mf[c(1L, m)]
	o = match(c("x"),names(c0))
	names(c0)[o] = "object"
	c0[[1L]] <- quote(stats::model.matrix)
	cv <- eval(c0, parent.frame())
	if(!is.na(match("(Intercept)",colnames(cv)))){
		cv = cv[,-1]
	}
	if(is.null(groups)){
		uv = .matrixUV(cv, psi = psi)
    } else {
		m <- match(c("groups", "data","contrasts"), names(mf), 0L)
		c1 <- mf[c(1L, m)]
        if(class(groups)=="formula"){
			atg <- attributes(stats::terms(groups))
            if(atg$intercept==1){
				c1[[2L]] = stats::formula(paste("~ -1 +",paste(atg$term.labels,sep="",collapse="+")))
            }
			o = match(c("groups"),names(c1))
			names(c1)[o] = "object"
			c1[[1L]] <- quote(stats::model.matrix)
			cc <- eval(c1, parent.frame())
			cc = apply(cc,1,paste,sep="")
        } else {
			cc = groups
        }
        cn = .getKnots(cc)
        cc = .incidenceMatrix(cc)
        colnames(cc) = cn
		if(ncol(cc)!=length(psi)) stop("Number of elements in psi must be equal to the number of groups in factor.")
		cv = apply(cc,2,function(x){return(x*cv)})
		uv = NULL
		for(i in 1:length(psi)){
			uv = cbind(uv,.matrixUV(cv[,i], psi = psi[i]))
		}
    }
	list(X = if(!only.UV){cv}, C = if(!is.null(groups) & !only.UV){cc}, UV = uv, psi = psi, call = mf)
}
s2D = function(data = NULL, penalty = c("none","cs","ps","tps"), is.X = c("none","tensor","tps"), intercept = TRUE,
                ps.order = 2, Ztr = c("RVA","RV","VA","V"), Amethod = c("eigen","svd"),
                Adiag = TRUE, aniso.angle = 0, aniso.ratio = 1, env = .GlobalEnv, ...){
	.defaultArg("penalty","none")
	.defaultArg("Ztr","RV")
	.defaultArg("Amethod","eigen")
    if(missing(is.X)) is.X <- NULL
    if(penalty=="none"){
        if(is.null(is.X))
            stop("If penalty=\"none\" then is.X must be specified and different than \"none\"!")
    } else if(penalty=="cs" | penalty=="ps"){
        if(is.null(is.X)){
            is.X = "tensor"
        } else if(is.X!="tensor"){
            is.X = "tensor"
            message("Changing is.X. If penalty=\"cs\" or \"ps\" then is.X must be \"tensor\"!.")
        }
    } else if(penalty=="tps"){
        if(is.null(is.X)){
            is.X = "tps"
        } else if(is.X!="tps"){
            is.X = "tps"
            message("Changing is.X. If penalty=\"tps\" then is.X must be \"tps\"!.")
        }
    }
    if(missing(data)) data <- NULL
	if(is.null(data)){
		ol = ls(envir = env)
		ol.id = match(c("data"),ol)
		if(all(is.na(ol.id))){
			stop("Input a sss (data) object.")
		} else {
			if(!is.na(match("data",ol))){
				data = get("data", envir = env)
			}
		}
	} else if(!is.sss(data)){
		stop("Input a sss (data) object.")
    }
    if(penalty=="none") smooth <- FALSE else smooth <- TRUE
    if(smooth){
        if(penalty=="ps" | penalty=="cs"){
        	c.Q = .tensorBase(obj = data, grid.list = attributes(data@grid)$grid.list, type = penalty, diff.order = ps.order)
        } else if(penalty=="tps"){
        	c.Q = .tpsBase(data@grid, ...)
        }
    } else {
        c.Q <- NULL
    }
    if(is.X=="none"){
        c.X = NULL
    } else if(is.X=="tensor" | is.X=="tps"){
    	c.xym = data@coords
    	c.X = as.matrix(cbind(if(intercept){1},as.matrix(c.xym)))
    	colnames(c.X) = c(if(intercept){"Intercept"},colnames(as.matrix(c.xym)))
        if(is.X=="tensor"){
        	aux.cX <- colnames(c.X)
        	c.X = as.matrix(cbind(c.X,apply(c.xym,1,prod)))
        	colnames(c.X) = c(aux.cX,paste(colnames(c.xym),sep="",collapse="."))
        }
    }
	c.h = .h2D(as.matrix(data@grid), angle = aniso.angle, ratio = aniso.ratio)$h
	c.use = list(penalty = penalty, is.X = is.X, intercept = intercept, ps.order = ps.order, Ztr = Ztr, Amethod = Amethod, Adiag = ifelse(penalty=="ps", FALSE, Adiag), aniso.angle = aniso.angle, aniso.ratio = aniso.ratio)
	list(X = c.X, h = c.h, Q = c.Q, use = c.use)
}
s2D = function(data = NULL, penalty = c("none","cs","ps","tps"), is.X = c("none","tensor","tps"), intercept = TRUE,
                ps.order = 2, aniso.angle = 0, aniso.ratio = 1, env = .GlobalEnv, ...){
	.defaultArg("penalty","none")
	Ztr <- "V"
	Amethod <- "eigen"
	Adiag = TRUE
    if(missing(is.X)) is.X <- NULL
    if(penalty=="none"){
        if(is.null(is.X))
            stop("If penalty=\"none\" then is.X must be specified and different than \"none\"!")
    } else if(penalty=="cs" | penalty=="ps"){
        if(is.null(is.X)){
            is.X = "tensor"
        } else if(is.X!="tensor"){
            is.X = "tensor"
            message("Changing is.X. If penalty=\"cs\" or \"ps\" then is.X must be \"tensor\"!.")
        }
    } else if(penalty=="tps"){
        if(is.null(is.X)){
            is.X = "tps"
        } else if(is.X!="tps"){
            is.X = "tps"
            message("Changing is.X. If penalty=\"tps\" then is.X must be \"tps\"!.")
        }
    }
    if(missing(data)) data <- NULL
	if(is.null(data)){
		ol = ls(envir = env)
		ol.id = match(c("data"),ol)
		if(all(is.na(ol.id))){
			stop("Input a sss (data) object.")
		} else {
			if(!is.na(match("data",ol))){
				data = get("data", envir = env)
			}
		}
	} else if(!is.sss(data)){
		stop("Input a sss (data) object.")
    }
    if(penalty=="none") smooth <- FALSE else smooth <- TRUE
	l.Z <- NULL
    if(smooth){
        if(penalty=="ps" | penalty=="cs"){
        	c.Q = .tensorBase(obj = data, grid.list = attributes(data@grid)$grid.list, type = penalty, diff.order = ps.order)
			l.Z <- list(c.Q$Qx$V%*%solve(t(c.Q$Qx$V)%*%c.Q$Qx$V), c.Q$Qy$V%*%solve(t(c.Q$Qy$V)%*%c.Q$Qy$V))
			names(l.Z) <- names(data@coords)
        } else if(penalty=="tps"){
        	c.Q = .tpsBase(data@grid, ...)
        }
    } else {
        c.Q <- NULL
    }
	l.X <- NULL
    if(is.X=="none"){
        c.X = NULL
    } else if(is.X=="tensor" | is.X=="tps"){
    	c.xym = data@coords
    	c.X = as.matrix(cbind(if(intercept){1},as.matrix(c.xym)))
    	colnames(c.X) = c(if(intercept){"Intercept"},colnames(as.matrix(c.xym)))
        if(is.X=="tensor"){
        	aux.cX <- colnames(c.X)
        	c.X = as.matrix(cbind(c.X,apply(c.xym,1,prod)))
        	colnames(c.X) = c(aux.cX,paste(colnames(c.xym),sep="",collapse="."))
        }
		l.X <- list(cbind(if(intercept){1},.getKnots(data@coords[,1])), cbind(if(intercept){1},.getKnots(data@coords[,2])))
		names(l.X) <- names(data@coords)
    }
	l.F <- NULL
	if(smooth){
		if(penalty!="tps"){
			l.F <- list(X = l.X[[2]]%x%l.X[[1]], Z = cbind(l.X[[2]]%x%l.Z[[1]], l.Z[[2]]%x%l.X[[1]], l.Z[[2]]%x%l.Z[[1]]))
		}
	}
	c.h = .h2D(as.matrix(data@grid), angle = aniso.angle, ratio = aniso.ratio)$h
	c.use = list(penalty = penalty, is.X = is.X, intercept = intercept, ps.order = ps.order, Ztr = Ztr, Amethod = Amethod, Adiag = ifelse(penalty=="ps", FALSE, Adiag), aniso.angle = aniso.angle, aniso.ratio = aniso.ratio)
	list(F = l.F, X = c.X, h = c.h, Q = c.Q, use = c.use)
}
.scp <- function(formula, data,
					initial = NULL, contrasts = NULL,
					pars = list(fix.nugget = FALSE, fix.kappa = FALSE, is.log.rhos = TRUE,
								model = "exponential", nugget.tol = 1.0e-15),
					method = list(use.reml = FALSE, use.profile = TRUE, linearize = "eigen"),
					control.optim = list(), control.gs = list(), control.ch = list(), print.pars = TRUE){
	if(missing(pars) | is.null(pars)){
		pars = list(fix.nugget = NULL, fix.kappa = NULL, is.log.rhos = TRUE,
					model = "matern", nugget.tol = 1.0e-15)
	}
	if(missing(method) | is.null(method)){
		method = list(use.reml = FALSE, use.profile = TRUE, linearize = "eigen")
	}
	if(!is.null(pars$fix.nugget)){
		if(is.logical(pars$fix.nugget)){
			if(pars$fix.nugget){
				warning("Not value set for fix.nugget!!\n\nTo fix the nugget to a value different than zero set fix.nugget=value.\nTo estimate the nugget set fix.nugget=FALSE or NULL.")
				pars$fix.nugget = 0
			} else {
				pars$fix.nugget = NULL
			}
		}
	}
	if(!is.null(pars$fix.kappa)){
		if(is.logical(pars$fix.kappa)){
			if(all(pars$fix.kappa)){
				warning("Not value set for fix.kappa!!\n\nTo fix kappa to a value different than 0.5 set fix.kappa=value.\nTo estimate kappa set fix.kappa=FALSE or NULL.")
				pars$fix.kappa = rep(0.5,length(pars$fix.kappa))
			} else {
				pars$fix.kappa = NULL
			}
		}
	}
	call.v = match.call()
	requireNamespace("stats", quietly = TRUE)
	mf = stats::terms(formula, specials=c("linear","cp","s1D","s2D"))
	af = attributes(mf)
	cr = rownames(af$factors)
	wl = grepl("linear\\(",cr)
	if(sum(wl)==1){
		cl = as.pairlist(parse(text=cr[af$specials$linear]))[[1]]
		if(is.null(cl$data)) cl$data <- data@data
		if(is.null(cl$contrasts)) cl$contrasts <- contrasts
		el = eval(cl, parent.frame())
	} else el = NULL
	wcp = grepl("cp\\(",cr)
	if(sum(wcp)>0){
		ec = NULL
		for(i in 1:sum(wcp)){
			cc = as.pairlist(parse(text=cr[af$specials$cp][i]))[[1]]
			if(is.null(cc$data)) cc$data <- data@data
			if(is.null(cc$contrasts)) cc$contrasts <- contrasts
			cc$only.UV <- FALSE
			ec[[i]] = eval(cc, parent.frame())
		}
	} else ec = NULL
	ws1D = grepl("s1D\\(",cr)
	if(sum(ws1D)>0){
		eg = NULL
		for(i in 1:sum(ws1D)){
			cg = as.pairlist(parse(text=cr[af$specials$s1D]))[[1]]
			if(is.null(cg$data)) cg$data <- data@data
			eg[[i]] = eval(cg, parent.frame())
		}
	} else eg = NULL
	ws2D = grepl("s2D\\(",cr)
	if(sum(ws2D)==1){
		cs = as.pairlist(parse(text=cr[af$specials$s2D]))[[1]]
		if(!is.null(el)){
			if(any(grepl("\\(Intercept\\)",colnames(el)))){
				if(ncol(el)>1){
					nel <- colnames(el)
					el <- as.matrix(el[,-match("(Intercept)",nel)])
					colnames(el) <- nel[-match("(Intercept)",nel)]
					cs$intercept <- TRUE
				} else {
					el <- NULL
				}
			}
		}
		if(is.null(cs$data)) cs$data <- data
		if(is.null(cs$aniso.angle)){ if(!is.null(pars$aniso.angle)){ cs$aniso.angle <- pars$aniso.angle } else { cs$aniso.angle <- 0 }} else { if(is.null(pars$aniso.angle)){ pars$aniso.angle <- cs$aniso.angle } }
		if(is.null(cs$aniso.ratio)){ if(!is.null(pars$aniso.ratio)){ cs$aniso.ratio <- pars$aniso.ratio } else { cs$aniso.ratio <- 1 }} else { if(is.null(pars$aniso.ratio)){ pars$aniso.ratio <- cs$aniso.ratio } }
		es = eval(cs, parent.frame())
	} else {
        es = list()
        es$use = list(penalty = "none", is.X = "none")
        if(is.null(pars$aniso.angle)) pars$aniso.angle <- 0
        if(is.null(pars$aniso.ratio)) pars$aniso.ratio <- 1
        if(is.null(pars$inverse)) pars$inverse <- "solve"
        if(is.null(pars$tol)) pars$tol <- .Machine$double.neg.eps
        es$use$inverse <- pars$inverse
        es$use$tol <- pars$tol
    	es$h = .h2D(as.matrix(data@grid), angle = pars$aniso.angle, ratio = pars$aniso.ratio)$h
		es$F = list(X = NULL, Z = NULL)
        es$X = NULL
        es$Q = NULL
    }
	zres = stats::formula(paste(cr[1],"~ 1"))
	cz = match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(cz), 0L)
    cz <- cz[c(1L, m)]
	cz$data <- data@data
    cz$formula <- zres
    cz[[1L]] <- quote(stats::model.frame)
    zV <- eval(cz, parent.frame())
	zV <- stats::model.response(zV)
	XM = cbind(el)
	VM <- NULL
	ZM <- NULL
	if(!is.null(es$F)){
		if(!is.null(es$F$Z)){
			ZM = cbind(es$F$Z)
			VM <- .generateV(ZM, itol = pars$tol)
		}
	}
	cpc = NULL
	cppsi = NULL
	UVel = NULL
	if(length(ec)>0){
		for(i in 1:length(ec)){
            ec.aux = ifelse(is.null(XM),0,ncol(XM))
			XM = cbind(XM,ec[[i]]$UV)
			cpc[[i]] = ec[[i]]$call
			cppsi[[i]] = ec[[i]]$psi
			UVel[[i]] = (ec.aux+1):ncol(XM)
		}
		ch.l = list(cp.call = cpc, psi0 = cppsi, UVcols = UVel)
	} else {
		ch.l = NULL
	}
	for(i in 1:length(eg)){
		XM = cbind(XM,eg[[i]]$X)
	}
		if(!is.null(XM)){
			XM = cbind(XM,es$X)
		} else {
			XM = cbind(es$X)
		}
		message("Starting computation")
		fit = .scpInitialize(zV = zV, XM = XM, ZM = ZM, VM = VM, hM = es$h, Qel = es$Q, W = data@W, ch = ch.l,
					use = es$use, pars = pars, initial = initial,
					method = method, control.optim = control.optim,
					control.gs = control.gs, control.ch = control.ch, print.pars = print.pars)
		message("Ending computation")
	message("Organising output")
	fit$penalty$matrices = es$Q
	out <- sssFit()
	out@data <- data
	out@zV <- zV
	out@XL <- ec
	out@XC <- ec
	out@XF <- eg
	out@XS <- es
	out@fit <- fit
	out@call <- call.v
	out
}
scp <- function(formula, data, initial = NULL, contrasts = NULL,
				model = "exponential", fix.nugget = FALSE, fix.kappa = FALSE,
                nugget.tol = 1.0e-15, angle = 0, ratio = 1,
				use.reml = FALSE, use.profile = TRUE, chMaxiter = 20, control = list()){
	mc <- match.call()
	if(missing(model)){ model <- "exponential" } else { if(is.null(model)){ model <- "exponential" } }
	pars <- list(fix.nugget = fix.nugget, fix.kappa = fix.kappa, is.log.rhos = TRUE, model = model, nugget.tol = nugget.tol, aniso.angle = angle, aniso.ratio = ratio)
	if(!is.null(control$tol))
		pars$tol <- control$tol
	method <- list(use.reml = use.reml, use.profile = use.profile)
	control.gs <- control[match(c("is.zero", "is.inf","size.each.grid","size.sample"),names(control))]
	control.ch <- list(maxiter = chMaxiter, maxskip = 10)
	control.optim <- control[match(c("trace", "fnscale","parscale","ndeps","maxit","abstol","reltol","alpha","beta","gamma","REPORT","type","lmm","factr","pgtol","temp","tmax"),names(control))]
	out <- .scp(formula, data,	initial, contrasts, pars = pars, method = method, control.optim = control.optim, control.gs = control.gs, control.ch = control.ch, print.pars = FALSE)
	out@call <- mc
	out
}
.scpLinearApprox <- function(obj,c,proj.method=c("Schur","proj","H-proj","ginv"),proj.to=c("sxz","sz")){
    .defaultArg("proj.method","H-proj")
    .defaultArg("proj.to","sz")
    X0 <- obj$fit$X
    R0 <- obj$fit$R
    Sigma <- obj$fit$sigma2s%x%R0
    H0 <- solve(R0)
    P0 <- X0%*%solve(t(X0)%*%H0%*%X0)%*%t(X0)%*%H0
    In <- diag(1,nrow(P0))
    Vdim <- ifelse(obj$XS$use$is.X=="tensor",4,ifelse(obj$XS$use$is.X=="tps",3,0))
    CM <- X0[,-c((ncol(X0)-(Vdim-1)):ncol(X0))]
    Xs <- X0[,(ncol(X0)-Vdim+1):ncol(X0)]
    Zs <- obj$fit$Zs
    pz <- Zs%*%solve(t(Zs)%*%H0%*%Zs + sum(obj$fit$alpha)%x%obj$fit$penalty$decompose$Ai)%*%t(Zs)%*%H0
    ipz <- In - pz
    px <- Xs%*%solve(t(Xs)%*%H0%*%ipz%*%Xs)%*%t(Xs)%*%H0%*%ipz
    ipx <- In - px
    sxz <- px + pz%*%ipx
    requireNamespace("Matrix",quietly=TRUE)
    if(proj.to=="sz"){
        ss <- Matrix::Schur(pz)
        oo <- order(ss$EValues,decreasing=TRUE)
        Qs <- ss$Q[,oo]
        Qsc <- Qs[,1:c]
        if(proj.method=="Schur")
            SS <- Qsc%*%t(Qsc)
        if(proj.method=="proj")
            SS <- Qsc%*%solve(t(Qsc)%*%Qsc)%*%t(Qsc)
        if(proj.method=="H-proj")
            SS <- Qsc%*%solve(t(Qsc)%*%H0%*%Qsc)%*%t(Qsc)%*%H0
        if(proj.method=="ginv")
            SS <- Qsc%*%.whichInverse("ginv",Qsc)
        SZ <- SS
        ISZ <- In - SZ
        SX <- Xs%*%solve(t(Xs)%*%H0%*%ISZ%*%Xs)%*%t(Xs)%*%H0%*%ISZ
        ISX <- In - SX
        SS <- SX + SZ%*%ISX
    } else {
        SX <- NULL
        SZ <- NULL
        ss <- Matrix::Schur(sxz)
        oo <- order(ss$EValues,decreasing=TRUE)
        Qs <- ss$Q[,oo]
        Qsc <- Qs[,1:c]
        if(proj.method=="Schur")
            SS <- Qsc%*%t(Qsc)
        if(proj.method=="proj")
            SS <- Qsc%*%solve(t(Qsc)%*%Qsc)%*%t(Qsc)
        if(proj.method=="H-proj")
            SS <- Qsc%*%solve(t(Qsc)%*%H0%*%Qsc)%*%t(Qsc)%*%H0
        if(proj.method=="ginv")
            SS <- Qsc%*%.whichInverse("ginv",Qsc)
    }
    if(ncol(CM)>0){
        PS <- CM%*%solve(t(CM)%*%H0%*%(In-SS)%*%CM)%*%t(CM)%*%H0%*%(In-SS)
        MS <- PS + SS%*%(In-PS)
    } else {
        PS <- NULL
        MS <- SS
    }
    B0 <- t(In-P0)%*%solve(Sigma)%*%(In-P0)
    SS0 <- t(obj$zV)%*%B0%*%obj$zV
    df0 <- sum(diag(B0%*%Sigma))
    B1 <- t(In-MS)%*%solve(Sigma)%*%(In-MS)
    SS1 <- t(obj$zV)%*%B1%*%obj$zV
    df1 <- sum(diag(B1%*%Sigma))
    criteria <- (df0 >= df1)
    null <- ifelse(criteria,0,1)
    if(criteria){
        B01 <- B0-B1
    } else {
        B01 <- B1-B0
    }
    SS01 <- t(obj$zV)%*%B01%*%obj$zV
    df01 <- sum(diag(B01%*%Sigma))
    Ft <- ifelse(criteria,(df1/df01)*(SS01/SS1),(df0/df01)*(SS01/SS0))
    pv <- 1 - pf(Ft,df01,ifelse(criteria,df1,df0))
    requireNamespace("mvtnorm",quietly=TRUE)
    loglik0 <- mvtnorm::dmvnorm(as.vector(obj$zV-P0%*%obj$zV),mean=rep(0,length(obj$zV)),sigma=Sigma,log=TRUE)
    aic0 <- -2 * loglik0 - 2 * sum(diag(P0))
    bic0 <- -2 * loglik0 - sum(diag(P0)) * log(length(obj$zV))
    logliks <- mvtnorm::dmvnorm(as.vector(obj$zV-MS%*%obj$zV),mean=rep(0,length(obj$zV)),sigma=Sigma,log=TRUE)
    aics <- -2*logliks - 2*sum(diag(SS))
    bics <- -2*logliks - sum(diag(SS))*log(length(obj$zV))
    o <- list(c=c,null=null,"tr(M0)"=sum(diag(P0)),"tr(S(phi))"=sum(diag(obj$fit$S)),"tr(S(c))"=sum(diag(SS)),"tr(M(phi))"=sum(diag(obj$fit$M)),"tr(M(c))"=sum(diag(MS)),F=round(Ft,5),p=round(pv,6),loglik0=loglik0,aic0=aic0,bic0=bic0,loglik=logliks,aic=aics,bic=bics,df0=df0,df1=round(df1,5),df01=round(df01,5),SS0=round(SS0,5),SS1=round(SS1,5),SS01=round(SS01,5))
    attr(o,"SX") <- SX
    attr(o,"SZ") <- SZ
    attr(o,"PX") <- px
    attr(o,"PZ") <- pz
    attr(o,"PS") <- PS
    attr(o,"A") <- Qsc
    attr(o,"Sigma") <- Sigma
    attr(o,"M0") <- P0
    attr(o,"P1") <- SS
    attr(o,"M1") <- MS
    attr(o,"B0") <- B0
    attr(o,"B1") <- B1
    attr(o,"B01") <- B01
    attr(o,"H0") <- ifelse(null==0,"Linear","Spline")
    attr(o,"H1") <- ifelse(null==0,"Spline","Linear")
    return(o)
}
.scpLinearApprox <- function (object, c, proj.method = c("Schur", "proj", "H-proj", "ginv"), proj.to = c("sxz", "sz"), sigma = NULL, R = NULL, tol = .Machine$double.neg.eps*1.0e-10){
    .defaultArg("proj.method", "H-proj")
    .defaultArg("proj.to", "sz")
    if(is.null(object@zV)){
        aux <- object
        object <- sssFit()
        object@zV <- aux$y
        object@fit <- list()
        object@fit$X <- aux$Xss
        object@fit$R <- aux$R
        object@fit$sigma2s <- aux$sigma2
        object@fit$W <- aux$W
        object@fit$Zs <- aux$W%*%aux$Z
        object@fit$alpha <- aux$alpha
        object@fit$psi <- aux$psi
        object@fit$M <- aux$M
        object@fit$S <- aux$Sa
        object@fit$penalty <- list()
        object@fit$penalty$decompose <- list()
        object@fit$penalty$decompose$Ai <- solve(.csBase(.getKnots(object@fit$X[,ncol(object@fit$X)]))$A,tol = tol)
        object@XS <- list()
        object@XS$use <- list()
        object@XS$use$is.X <- "uni"
    }
    if(!is.null(sigma) & !is.null(R)){
        object@fit$sigma2s <- sigma**2
        object@fit$R <- R
    }
    X0 <- object@fit$X
    R0 <- object@fit$R
    Sigma <- object@fit$sigma2s %x% R0
    H0 <- solve(R0,tol = tol)
    P0 <- X0 %*% solve(t(X0) %*% H0 %*% X0,tol = tol) %*% t(X0) %*% H0
    In <- diag(1, nrow(P0))
    Vdim <- ifelse(object@XS$use$is.X == "tensor", 4, ifelse(object@XS$use$is.X == "tps", 3, ifelse(object@XS$use$is.X == "uni", 2, 0)))
    CM <- cbind(X0[, -c((ncol(X0) - (Vdim - 1)):ncol(X0))])
    if(object@XS$use$is.X != "uni"){
        Xs <- X0[, (ncol(X0) - Vdim + 1):ncol(X0)]
    } else {
        Xs <- X0[, (ncol(X0) - Vdim + 1):ncol(X0)]
    }
    Zs <- object@fit$Zs
    pz <- Zs %*% solve(t(Zs) %*% H0 %*% Zs + sum(object@fit$alpha) %x% object@fit$penalty$decompose$Ai,tol = tol) %*% t(Zs) %*% H0
    ipz <- In - pz
    px <- Xs %*% solve(t(Xs) %*% H0 %*% ipz %*% Xs,tol = tol) %*% t(Xs) %*% H0 %*% ipz
    ipx <- In - px
    sxz <- px + pz %*% ipx
    requireNamespace("Matrix", quietly = TRUE)
    if (proj.to == "sz") {
        ss <- Matrix::Schur(pz)
        oo <- order(ss$EValues, decreasing = TRUE)
        Qs <- ss$Q[, oo]
        Qsc <- Qs[, 1:c]
        if (proj.method == "Schur")
            SS <- Qsc %*% t(Qsc)
        if (proj.method == "proj")
            SS <- Qsc %*% solve(t(Qsc) %*% Qsc,tol = tol) %*% t(Qsc)
        if (proj.method == "H-proj")
            SS <- Qsc %*% solve(t(Qsc) %*% H0 %*% Qsc,tol = tol) %*% t(Qsc) %*% H0
        if (proj.method == "ginv")
            SS <- Qsc %*% .whichInverse("ginv", Qsc)
        SZ <- SS
        ISZ <- In - SZ
        SX <- Xs %*% solve(t(Xs) %*% H0 %*% ISZ %*% Xs,tol = tol) %*% t(Xs) %*% H0 %*% ISZ
        ISX <- In - SX
        SS <- SX + SZ %*% ISX
    }
    else {
        SX <- NULL
        SZ <- NULL
        ss <- Matrix::Schur(sxz)
        oo <- order(ss$EValues, decreasing = TRUE)
        Qs <- ss$Q[, oo]
        Qsc <- Qs[, 1:c]
        if (proj.method == "Schur")
            SS <- Qsc %*% t(Qsc)
        if (proj.method == "proj")
            SS <- Qsc %*% solve(t(Qsc) %*% Qsc,tol = tol) %*% t(Qsc)
        if (proj.method == "H-proj")
            SS <- Qsc %*% solve(t(Qsc) %*% H0 %*% Qsc,tol = tol) %*% t(Qsc) %*% H0
        if (proj.method == "ginv")
            SS <- Qsc %*% .whichInverse("ginv", Qsc)
    }
    if (ncol(CM) > 0) {
        PS <- CM %*% solve(t(CM) %*% H0 %*% (In - SS) %*% CM,tol = tol) %*% t(CM) %*% H0 %*% (In - SS)
        MS <- PS + SS %*% (In - PS)
    }
    else {
        PS <- NULL
        MS <- SS
    }
    B0 <- t(In - P0) %*% solve(Sigma,tol = tol) %*% (In - P0)
    SS0 <- t(object@zV) %*% B0 %*% object@zV
    df0 <- sum(diag(B0 %*% Sigma))
    B1 <- t(In - MS) %*% solve(Sigma,tol = tol) %*% (In - MS)
    SS1 <- t(object@zV) %*% B1 %*% object@zV
    df1 <- sum(diag(B1 %*% Sigma))
    criteria <- (df0 >= df1)
    null <- ifelse(criteria, 0, 1)
    if (criteria) {
        B01 <- B0 - B1
    }
    else {
        B01 <- B1 - B0
    }
    SS01 <- t(object@zV) %*% B01 %*% object@zV
    df01 <- sum(diag(B01 %*% Sigma))
    Ft <- ifelse(criteria, (df1/df01) * (SS01/SS1), (df0/df01) * (SS01/SS0))
    pv <- 1 - pf(Ft, df01, ifelse(criteria, df1, df0))
    requireNamespace("mvtnorm", quietly = TRUE)
    loglik0 <- mvtnorm::dmvnorm(as.vector(object@zV - P0 %*% object@zV),
        mean = rep(0, length(object@zV)), sigma = Sigma, log = TRUE)
    aic0 <- -2 * loglik0 - 2 * sum(diag(P0))
    bic0 <- -2 * loglik0 - sum(diag(P0)) * log(length(object@zV))
    logliks <- mvtnorm::dmvnorm(as.vector(object@zV - MS %*% object@zV),
        mean = rep(0, length(object@zV)), sigma = Sigma, log = TRUE)
    aics <- -2 * logliks - 2 * sum(diag(SS))
    bics <- -2 * logliks - sum(diag(SS)) * log(length(object@zV))
    ndig <- nchar(strsplit(paste(pv),"\\.")[[1]][2])
    o <- list(D = length(object@fit$psi), c = c, null = null, `tr(M0)` = sum(diag(P0)), `tr(S(phi))` = sum(diag(object@fit$S)),
        `tr(S(c))` = sum(diag(SS)), `tr(M(phi))` = sum(diag(object@fit$M)), `tr(M(c))` = sum(diag(MS)),
		F = round(Ft, 5), p = ifelse(ndig>6,signif(pv),round(pv, 6)), loglik0 = loglik0, aic0 = aic0, bic0 = bic0, loglik = logliks,
		aic = aics, bic = bics, df0 = df0, df1 = round(df1, 5), df01 = round(df01, 5), SS0 = round(SS0, 5),
		SS1 = round(SS1, 5), SS01 = round(SS01, 5))
    attr(o, "SX") <- SX
    attr(o, "SZ") <- SZ
    attr(o, "PX") <- px
    attr(o, "PZ") <- pz
    attr(o, "PS") <- PS
    attr(o, "A") <- Qsc
    attr(o, "Sigma") <- Sigma
    attr(o, "sigma") <- sqrt(object@fit$sigma2s)
    attr(o, "M0") <- P0
    attr(o, "P1") <- SS
    attr(o, "M1") <- MS
    attr(o, "B0") <- B0
    attr(o, "B1") <- B1
    attr(o, "B01") <- B01
    attr(o, "H0") <- ifelse(null == 0, "Linear", "Spline")
    attr(o, "H1") <- ifelse(null == 0, "Spline", "Linear")
    return(o)
}
.scpApproxGet <- function(object, tol){
	if(missing(tol)){ tol <- .Machine$double.neg.eps*1.0e-10 } else { if(is.null(tol)){ tol <- .Machine$double.neg.eps*1.0e-10 } }
	cVal <- abs(sum(diag(object@fit$S))-1:nrow(object@fit$M))
	cVal <- match(min(cVal),cVal)
	app <- .scpLinearApprox(object = object, c = cVal, tol = tol)
	att <- attributes(app)
	object@fit$fitted.values <- att$M1%*%object@zV
	object@fit$M <- att$M1
	if(!is.null(att$PS)){
		object@fit$g <- att$P1%*%(diag(1,nrow(att$P1))-att$PS)%*%object@zV
	} else {
		object@fit$g <- att$P1%*%object@zV
	}
	object@fit$approximate <- TRUE
	return(object)
}
.scpTestGet <- function(object, tol){
	if(missing(tol)){ tol <- .Machine$double.neg.eps*1.0e-10 } else { if(is.null(tol)){ tol <- .Machine$double.neg.eps*1.0e-10 } }
	if(!is.null(object@fit$approximate))
		if(object@fit$approximate) stop("Object must be a fit from scp command.")
	if(is.null(object@fit$g))
		stop("Model must include a spline component.")
	cVal <- abs(sum(diag(object@fit$S))-1:nrow(object@fit$M))
	cVal <- match(min(cVal),cVal)
	app <- .scpLinearApprox(object = object, c = cVal, tol = tol)
	att <- attributes(app)
	app$H0 <- att$H0
	app$H1 <- att$H1
	app$G <- app$D
	app <- app[match(c("G","c","H0","H1","df0","df1","df01","SS0","SS1","SS01","F","p"),names(app))]
	tab <- cbind(c(app$df01,app$df1,app$df0),c(app$SS01,app$SS1,app$SS0),c(app$SS01/app$df01,app$SS1/app$df1,app$SS0/app$df0),c(app$F,NA,NA),c(app$p,NA,NA))
	rownames(tab) <- c(ifelse(app$H1=="Spline","Approx. over Linear","Linear over Approx."),ifelse(app$H1=="Spline","Approx.","Linear"),ifelse(app$H1=="Spline","Linear","Approx."))
	colnames(tab) <- c("DF", "SS", "MS", "F(y)", "Prob(F>F(y))")
	print(tab,na.print="")
}
.sssFit.valid <- function(x){
 el <- list(data="sss", zV="numeric", XL="ANY", XC="ANY", XF="ANY", XS="ANY", fit="list", call="call")
 requireNamespace("methods",quietly=TRUE)
 exists.slot <- match(names(el),methods::slotNames(x))
 if(!any(is.na(exists.slot))){
    if(!is.null(x@data)){
       if(!is.sss(x@data)){
          paste("data is not of class \"sss\".")
       }
    } else {
		paste("a data of class \"sss\" must be part of the object.")
	}
    if(length(as.vector(unlist(x@zV)))!=nrow(x@data@coords)){
		paste("Wrong structure in zV and data.")
	}
    if(TRUE){
       valid = NULL
       for(i in 1:length(el)){
			if(!is.na(match(el[[i]],c("sss","numeric","list","call")))){
                valid[i] <- any(grepl(paste("^",el[[i]],"$",sep=""),class(methods::slot(x,names(el)[[i]]))))
			} else {
				valid[i] = TRUE
			}
       }
     if(all(valid)) TRUE else paste("Slots",paste(names(el)[!valid],sep="",collapse=", "),"not match their class.")
    }
 } else {
  paste("Slots",paste(names(el)[is.na(exists.slot)],sep="",collapse=", "),"are missing.")
 }
}
sssFit <- setClass("sssFit",
	slots = list(data="sss", zV="numeric", XL="ANY", XC="ANY", XF="ANY", XS="ANY", fit="list", call="call"),
	package = "scpm", validity = .sssFit.valid)
.sssToSurface = function (data = NULL, grid.list = NULL, values = NULL, names = NULL){
	if(!is.null(data)){
		grid.list = attr(data@grid, "grid.list")
	}
	if(!is.null(grid.list)){
		if(is.matrix(grid.list)){
			grid.list = attr(grid.list, "grid.list")
		}
	}
	if(is.null(names)){
		if(!is.null(values)){
			list(x = grid.list[[1]], y = grid.list[[2]], z = matrix(values, ncol = length(grid.list[[2]]), nrow = length(grid.list[[1]]), byrow = FALSE), xlab = names(grid.list)[1], ylab = names(grid.list)[2])
		} else {stop("values or names must be specified!")}
	} else {
		out = NULL
		for(i in 1:length(names)){
			out[[i]] = list(x = grid.list[[1]], y = grid.list[[2]], z = matrix(as.vector(data@data[[names[i]]]), ncol = length(grid.list[[2]]), nrow = length(grid.list[[1]]), byrow = FALSE), xlab = names(grid.list)[1], ylab = names(grid.list)[2])
		}
		names(out) = names
		return(out)
	}
}
.add3Dpoints <- function(x,y,z, pmat) {
  tr <- cbind(x,y,z,1) %*% pmat
  list(x = tr[,1]/tr[,4], y= tr[,2]/tr[,4])
}
.legendMapSimple <- function(col, lev, val=NULL){
	requireNamespace("graphics", quietly = TRUE)
	opar <- graphics::par
	n <- length(col)
	range.lev = c(min(lev),max(lev))
	bx <- graphics::par("usr")
	box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
	bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
	box.cy <- c(bx[3], bx[3])
	box.sy <- (bx[4] - bx[3]) / n
	unit.lev = (bx[4] - bx[3]) / (range.lev[2] - range.lev[1])
	olev = lev[order(lev)]
	xx <- rep(box.cx, each = 2)
	xx = xx - (xx[4]-xx[1])/2
	graphics::par(xpd = TRUE)
	yb = NULL
	for(i in 1:n){
		yb[i] = box.cy[1] + box.sy/2
		yy <- c(box.cy[1] + (box.sy * (i - 1)),
		box.cy[1] + (box.sy * (i)),
		box.cy[1] + (box.sy * (i)),
		box.cy[1] + (box.sy * (i - 1)))
		yy = yy + box.sy/2
		yb[i] = mean(yy)
		graphics::polygon(xx, yy, col = col[i], border = col[i])
	}
	requireNamespace("stats", quietly = TRUE)
	lev.tick = stats::quantile(lev, probs = seq(0,1,0.1))
	lev.plot = lev[match(paste(lev.tick),paste(lev))]
	yb.plot = yb[match(paste(lev.tick),paste(lev))]
	valb.plot = min(yb) + (val-min(val))*(max(yb)-min(yb))/(max(val)-min(val))
	graphics::axis(side = 4, at = yb.plot, labels = round(lev.plot), las = 3, tick = TRUE, line = 1, tcl = -0.5, cex.axis = 0.85)
	graphics::axis(side = 4, at = valb.plot, labels = FALSE, las = NULL, tick = TRUE, line = 1, lwd.ticks = 2, tcl = -0.25)
	par <- opar
}
.legendMap <- function(col, lev, val=NULL, by.prob = 0.01, out = FALSE){
 requireNamespace("graphics", quietly = TRUE)
 opar <- graphics::par
 n <- length(col)
 range.lev = c(min(lev, na.rm=TRUE),max(lev, na.rm=TRUE))
 bx <- graphics::par("usr")
 box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
 bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
 box.cy <- c(bx[3], bx[3])
 box.sy <- (bx[4] - bx[3]) / n
 unit.lev = (bx[4] - bx[3]) / (range.lev[2] - range.lev[1])
 olev = lev[order(lev)]
 xx <- rep(box.cx, each = 2)
 xx = xx - (xx[4]-xx[1])/2
 graphics::par(xpd = TRUE)
 yb = NULL
 for(i in 1:n){
  yb[i] = box.cy[1] + box.sy/2
  yy <- c(box.cy[1] + (box.sy * (i - 1)),
  box.cy[1] + (box.sy * (i)),
  box.cy[1] + (box.sy * (i)),
  box.cy[1] + (box.sy * (i - 1)))
  yy = yy + box.sy/2
  yb[i] = mean(yy)
  graphics::polygon(xx, yy, col = col[i], border = col[i])
 }
 requireNamespace("stats", quietly = TRUE)
 lev.tick = stats::quantile(lev, probs = seq(0,1,by.prob), na.rm = TRUE)
 lev.plot = lev[match(paste(lev.tick),paste(lev))]
 yb.plot = yb[match(paste(lev.tick),paste(lev))]
 if(!is.null(val)){
  valb.plot = min(yb) + (val-min(val))*(max(yb)-min(yb))/(max(val)-min(val))
 }
 graphics::axis(side = 4, at = yb.plot, labels = round(lev.plot), las = 3, tick = TRUE, line = 1, tcl = -0.5, cex.axis = 0.85)
 if(!is.null(val)){
  graphics::axis(side = 4, at = valb.plot, labels = FALSE, las = NULL, tick = TRUE, line = 1, lwd.ticks = 2, tcl = -0.25)
 }
 par <- opar
 if(out){
  list(usr = yb)
 }
}
.colorVecMat <- function(obj, which = c("rainbow","heat.colors","terrain.colors","topo.colors","cm.colors","colorRamp","colorRampPalette"), col.args = list()){
	.defaultArg("which","rainbow")
	which.fx = c("colorRamp","colorRampPalette")
	if(is.matrix(obj)){
		values = as.vector(obj)
	} else {
		values = obj
	}
	requireNamespace("grDevices", quietly = TRUE)
	cargs = list()
	if(length(col.args)>0){
		for(i in 1:length(col.args)){
			cargs[[i]] = col.args[[i]]
		}
		names(cargs) = names(col.args)
	}
	if(!is.na(match(which,which.fx))){
		af = formals(which)
		fx.args = cargs[!is.na(match(names(cargs),names(af)))]
		fx.cols = do.call(which,fx.args)
		which = "fx.cols"
	}
	if(is.null(cargs$n)){
		cargs$n = length(values)
	}
	af = formals(which)
	cargs = cargs[!is.na(match(names(cargs),names(af)))]
	cols = do.call(which,cargs)
	if(is.matrix(obj)){
		zfacet <- obj[-1, -1] + obj[-1, -ncol(obj)] + obj[-nrow(obj), -1] + obj[-nrow(obj), -ncol(obj)]
		facetcol <- cut(zfacet, length(cols))
		return(cols[facetcol])
	} else {
		oz = values[order(values)]
		icols = match(paste(values),paste(oz))
		return(cols[icols])
	}
}
.scpImage <- function(values, data = NULL, grid.list = NULL,
	which = c("image","image.default","image.plot","image.kriging"),
	which.col = "rainbow", image.args = list(), contour.args = list(),
	legend.args = list(), col.args = list()){
	.defaultArg("which","image.default")
	if(is.list(values)){
		zsurf = values
	} else {
		zsurf = .sssToSurface(data = data, grid.list = grid.list, values = values, names = NULL)
	}
	v.v = as.vector(zsurf$z)
	v.m = zsurf$z
	if(is.null(col.args$start)){ col.args$start = min(.transfTo01(v.v)) }
	if(is.null(col.args$end)){ col.args$end = max(.transfTo01(v.v)) }
	if(is.null(image.args$col)){ image.args$col = .colorVecMat(obj = v.v[order(v.v)], which = which.col, col.args = col.args) }
	if(is.null(image.args$oldstyle)){ image.args$oldstyle = TRUE }
	if(is.null(image.args$useRaster)){ image.args$useRaster = TRUE }
	if(is.null(contour.args$method)){ contour.args$method = "flattest" }
	if(is.null(legend.args$col)){ legend.args$col = image.args$col }
	if(is.null(legend.args$lev)){ legend.args$lev = v.v[order(v.v)] }
	if(is.null(legend.args$val)){ legend.args$val = v.v[order(v.v)] }
	if(is.null(legend.args$by.prob)){ legend.args$by.prob = 0.01 }
	if(is.null(legend.args$out)){ legend.args$out = FALSE }
	af = formals(which)
	image.args = image.args[!is.na(match(names(image.args),names(af)))]
	af = formals(contour.default)
	contour.args = contour.args[!is.na(match(names(contour.args),names(af)))]
	iargs = zsurf
	cargs = zsurf
	for(i in 1:length(image.args)){
		iargs[[i+length(zsurf)]] = image.args[[i]]
	}
	names(iargs) = c(names(zsurf),names(image.args))
	if(grepl("^image.plot$",which))
		if(!requireNamespace("fields", quietly = TRUE)) which <- "image.default"
	if(grepl("^image.kriging$",which))
		which <- "image.default"
	if(grepl("^image(.default){0,1}$",which))
		requireNamespace("graphics", quietly = TRUE)
	do.call(which,iargs)
	for(i in 1:length(contour.args)){
		cargs[[i+length(zsurf)]] = contour.args[[i]]
	}
	names(cargs) = c(names(zsurf),names(contour.args))
	cargs$add = TRUE
	if(!isNamespaceLoaded("graphics"))
		requireNamespace("graphics", quietly = TRUE)
	do.call("contour.default",cargs)
	if(which!="image.plot"){
		do.call(".legendMap",legend.args)
	}
}
.scpPersp = function(values, data = NULL, grid.list = NULL,
	which = c("persp","persp.kriging","persp.grf"),
	which.col = "rainbow",	persp.args = list(), legend.args = list(), col.args = list()){
	.defaultArg("which","persp")
	if(is.list(values)){
		zsurf = values
	} else {
		zsurf = .sssToSurface(data = data, grid.list = grid.list, values = values, names = NULL)
	}
	v.v = as.vector(zsurf$z)
	v.m = zsurf$z
	if(is.null(col.args$start)){ col.args$start = min(.transfTo01(v.v)) }
	if(is.null(col.args$end)){ col.args$end = max(.transfTo01(v.v)) }
	if(is.null(persp.args$col)){ persp.args$col = .colorVecMat(v.m, which = which.col, col.args = col.args) }
	if(is.null(persp.args$theta)){ persp.args$theta = 0 }
	if(is.null(persp.args$phi)){ persp.args$phi = 60 }
	if(is.null(persp.args$ticktype)){ persp.args$ticktype = "detailed" }
	if(is.null(persp.args$border)){ persp.args$border = "transparent" }
	if(is.null(persp.args$shade)){ persp.args$shade = 0.15 }
	if(is.null(persp.args$zlab)){ persp.args$zlab = "z" }
	if(is.null(legend.args$col)){ legend.args$col = .colorVecMat(v.v[order(v.v)], which = which.col,
							col.args = col.args) }
	if(is.null(legend.args$lev)){ legend.args$lev = v.v[order(v.v)] }
	if(is.null(legend.args$val)){ legend.args$val = v.v[order(v.v)] }
	if(is.null(legend.args$by.prob)){ legend.args$by.prob = 0.01 }
	if(is.null(legend.args$out)){ legend.args$out = FALSE }
	pargs = zsurf
	for(i in 1:length(persp.args)){
		pargs[[i+length(zsurf)]] = persp.args[[i]]
	}
	names(pargs) = c(names(zsurf),names(persp.args))
	if(grepl("^persp.(grf|kriging){1}$",which))
		which <- "persp"
	if(grepl("^persp$",which))
		requireNamespace("graphics", quietly = TRUE)
	xypos = do.call(which,pargs)
	do.call(".legendMap",legend.args)
	return(xypos)
}
.scpPoints = function(obj, values = NULL, pos3d = NULL, points.args = list()){
	if(is.sss(obj)){
		coords = obj@coords
	}
	if(is.matrix(obj)){
		coords = obj
	}
	af = formals(points)
	points.args = points.args[!is.na(match(names(points.args),names(af)))]
	pargs = list()
	if(is.null(values) & is.null(pos3d)){
		pargs$x = coords[,1]
		pargs$y = coords[,2]
	} else {
		pargs = .add3Dpoints(coords[,1],coords[,2], values, pos3d)
	}
	if(length(points.args)>0){
		bar = length(pargs)
		nbar = names(pargs)
		for(i in 1:length(points.args)){
			pargs[[i+bar]] = points.args[[i]]
		}
		names(pargs) = c(nbar,names(points.args))
	}
	requireNamespace("graphics", quietly = TRUE)
	do.call("graphics::points",pargs)
}
.scpPlot = function(obj, data = NULL,
					which.col = "rainbow", which.image = "image.default", which.persp = "persp",
					image.args = list(), contour.args = list(), persp.args = list(), col.args = list(),
					add.points = FALSE,
					titles = list()){
	if(is.null(data)){
		data = obj$data
	}
	if(is.null(titles$z)){ titles$z = bquote(Observed~Field) }
	if(is.null(titles$mu)){ titles$mu = bquote(Estimated~Field,~E(Z)==mu[Z]) }
	if(is.null(titles$g)){ titles$g = bquote(Estimated~g(x,y)) }
	if(is.null(titles$cb)){ titles$cb = bquote(Estimated~Cb) }
	if(is.null(titles$xb)){ titles$xb = bquote(Estimated~X*beta) }
	add.one = 0
	c.mu = obj$fit$fitted.values
	if(ncol(obj$fit$X)>4){
		c.cb = obj$fit$X[,1:(ncol(obj$fit$X)-4)]%*%obj$fit$beta[1:(ncol(obj$fit$X)-4)]
		add.one = add.one + 1
	} else {
		c.cb = NULL
	}
	c.gf = obj$fit$g
	c.xb = obj$fit$X[,(ncol(obj$fit$X)-3):ncol(obj$fit$X)]%*%obj$fit$beta[(ncol(obj$fit$X)-3):ncol(obj$fit$X)]
	c.all = c(obj$zV, c.mu, c.cb, c.gf, c.xb)
	out.z = .sssToSurface(data = data, values = obj$zV)
	out.mu = .sssToSurface(data = data, values = c.mu)
	if(ncol(obj$fit$X)>4){
		out.cb = .sssToSurface(data = data, values = c.cb)
	}
	out.gf = .sssToSurface(data = data, values = c.gf)
	out.xb = .sssToSurface(data = data, values = c.xb)
	par(mfrow=c(2,ifelse(add.one>0,5,3)), mar=c(5,4,4,3)+0.1)
	.scpImage(values = out.z, which = which.image, which.col = which.col,
				image.args = image.args, contour.args = contour.args, col.args = col.args)
	title(main = titles$z)
	if(add.one>0){
		.scpImage(values = out.mu, which = which.image, which.col = which.col,
				image.args = image.args, contour.args = contour.args, col.args = col.args)
		title(main = titles$mu)
		.scpImage(values = out.gf, which = which.image, which.col = which.col,
				image.args = image.args, contour.args = contour.args, col.args = col.args)
		title(main = titles$g)
		.scpImage(values = out.cb, which = which.image, which.col = which.col,
				image.args = image.args, contour.args = contour.args, col.args = col.args)
		title(main = titles$cb)
		.scpImage(values = out.xb, which = which.image, which.col = which.col,
				image.args = image.args, contour.args = contour.args, col.args = col.args)
		title(main = titles$xb)
	} else {
		.scpImage(values = out.gf, which = which.image, which.col = which.col,
				image.args = image.args, contour.args = contour.args, col.args = col.args)
		title(main = titles$g)
		.scpImage(values = out.xb, which = which.image, which.col = which.col,
				image.args = image.args, contour.args = contour.args, col.args = col.args)
		title(main = titles$xb)
	}
	if(is.null(persp.args$zlim)){ persp.args$zlim = c(min(c.all),max(c.all)) }
	pp = .scpPersp(values = out.z, which = which.persp, which.col = which.col, persp.args = persp.args, col.args = col.args)
	title(main = titles$z)
	if(add.points){	.scpPoints(obj = data, values = obj$zV, pos3d = pp, points.args = list(pch = 1)) }
	if(add.one>0){
		pp = .scpPersp(values = out.mu, which = which.persp, which.col = which.col, persp.args = persp.args, col.args = col.args)
		title(main = titles$mu)
		if(add.points){	.scpPoints(obj = data, values = obj$zV, pos3d = pp, points.args = list(pch = 1)) }
		pp = .scpPersp(values = out.gf, which = which.persp, which.col = which.col, persp.args = persp.args, col.args = col.args)
		title(main = titles$g)
		if(add.points){	.scpPoints(obj = data, values = obj$zV, pos3d = pp, points.args = list(pch = 1)) }
		pp = .scpPersp(values = out.cb, which = which.persp, which.col = which.col, persp.args = persp.args, col.args = col.args)
		title(main = titles$cb)
		if(add.points){	.scpPoints(obj = data, values = obj$zV, pos3d = pp, points.args = list(pch = 1)) }
		pp = .scpPersp(values = out.xb, which = which.persp, which.col = which.col, persp.args = persp.args, col.args = col.args)
		title(main = titles$xb)
		if(add.points){	.scpPoints(obj = data, values = obj$zV, pos3d = pp, points.args = list(pch = 1)) }
	} else {
		pp = .scpPersp(values = out.gf, which = which.persp, which.col = which.col, persp.args = persp.args, col.args = col.args)
		title(main = titles$g)
		if(add.points){	.scpPoints(obj = data, values = obj$zV, pos3d = pp, points.args = list(pch = 1)) }
		pp = .scpPersp(values = out.xb, which = which.persp, which.col = which.col, persp.args = persp.args, col.args = col.args)
		title(main = titles$xb)
		if(add.points){	.scpPoints(obj = data, values = obj$zV, pos3d = pp, points.args = list(pch = 1)) }
	}
}
.getSemivariogram <- function(x,obj){
    .semiVarF(h=x,phi=obj@fit$phi,kappa=obj@fit$kappa,model=eval(obj@call$model),sill=obj@fit$sigma2,nugget=obj@fit$tau2,tol.nugget=ifelse(is.null(eval(obj@call$nugget.tol)),1.0e-15,eval(obj@call$nugget.tol)),use.cor=TRUE)
}
.plotSemivariogram <- function(obj,h=NULL,hmax=NULL,size=200,plot=TRUE,prange=TRUE,legend=TRUE,...){
    aux.tol <- ifelse(is.null(eval(obj@call$nugget.tol)),1.0e-15,eval(obj@call$nugget.tol))
    if(is.null(hmax) & is.null(h))
        hmax <- max(abs(obj@XS$h)) + (max(abs(obj@XS$h))-min(abs(obj@XS$h)))*0.05
	if(is.null(hmax) & !is.null(h))
		hmax <- max(h)
	if(is.null(h)){
		xv <- seq(0,hmax,l=size)
	} else {
		xv <- h
	}
    sv <- .getSemivariogram(xv,obj)
    xv <- xv[!is.na(sv)]
    sv <- sv[!is.na(sv)]
    xv <- c(aux.tol,xv)
    sv <- c(obj@fit$tau2,sv)
    sv <- sv[order(xv)]
    xv <- xv[order(xv)]
    rr <- xv[match(max(sv),sv)]
    rr <- (xv[sv==as.numeric(obj@fit$sigma2s)])[1]
    rp <- xv[sv <= as.numeric(obj@fit$sigma2s)*0.95]
    rp <- rp[length(rp)]
    if(plot){
        message(paste("Sill=",round(obj@fit$sigma2s,6),sep=""))
        message(paste("Practical range h",ifelse(length(rp)==0,paste(">",hmax,sep=""),paste("=",round(rp,6),sep="")),", 0.95*sill=",round(as.numeric(obj@fit$sigma2s)*0.95,6),sep=""))
        gargs <- list(...)
        if(is.null(gargs$xlim)){xlim <- c(0,max(xv))} else {xlim <- gargs$xlim}
        if(is.null(gargs$ylim)){ylim <- c(0,max(sv))} else {ylim <- gargs$ylim}
        if(is.null(gargs$main)){main <- paste("Fitted semi-variogram model",eval(obj@call$model))} else {main <- gargs$main}
        plot(sv[xv>=aux.tol]~xv[xv>=aux.tol],ylab=expression(hat(gamma(h))),xlab=expression(h),type="l",xlim=xlim,ylim=ylim,main=main)
        points(x=aux.tol,y=obj@fit$tau2,pch=16)
        if(prange){
            rr <- NA
            abline(v=c(if(!is.na(rr)){rr},if(length(rp)!=0){rp}),lty=c(if(!is.na(rr)){4},if(length(rp)!=0){3}),col=c(if(!is.na(rr)){"blue"},if(length(rp)!=0){"blue"}))
            if(legend)
                legend("bottomright",legend=c(if(!is.na(rr)){as.expression(bquote(Range==.(round(rr,4))))},if(length(rp)!=0){as.expression(bquote(Prac.*phantom(0)*Range==.(round(rp,4))))}),lty=c(if(!is.na(rr)){4},if(length(rp)!=0){3}),col=c(if(!is.na("blue")){rr},if(length(rp)!=0){"blue"}),bg="white")
        }
    } else {
        list(x = xv[xv>=aux.tol], y = sv[xv>=aux.tol], sill = obj@fit$sigma2s, prange = rp, psill = obj@fit$sigma2s*0.95)
    }
}
.scpSummary <- function(object,alpha=0.05){
	if(!is.null(object@fit$approximate))
		if(object@fit$approximate) stop("Object must be a fit from scp command.")
	obj <- object
    tv <- obj@fit$beta/sqrt(diag(obj@fit$var.beta))
    ptv <- 2*(1-pt(abs(tv),length(obj@zV)-ncol(obj@fit$X)))
    atv <- qt(1-alpha/2,length(obj@zV)-ncol(obj@fit$X))*sqrt(diag(obj@fit$var.beta))
    ctable <- cbind(obj@fit$beta,sqrt(diag(obj@fit$var.beta)),tv,ptv,obj@fit$beta-atv,obj@fit$beta+atv)
    colnames(ctable) <- c("estimate","std.error","t(y)","p(|t(y)|>t)","LL(95%)","UL(95%)")
	nr <- rownames(ctable)
	nr <- gsub("^Intercept$","(Intercept)",nr)
	rownames(ctable) <- nr
    return(ctable)
}
.scpImagePlot <- function(obj, type = c("obs","fit","g","both"),...){
	if(missing(type)){ type <- "fit"} else { if(is.null(type)){ type <- "fit" } }
    type <- match.arg(type)
    requireNamespace("interp", quietly = TRUE)
	intz <- interp::interp(x=obj@data@coords[,1], y=obj@data@coords[,2], z=obj@zV, method = "linear", linear = TRUE, extrap = FALSE)
    intm <- interp::interp(x=obj@data@coords[,1], y=obj@data@coords[,2], z=obj@fit$fitted.values, method = "linear", linear = TRUE, extrap = FALSE)
    if(!is.null(obj@fit$g))
        intg <- interp::interp(x=obj@data@coords[,1], y=obj@data@coords[,2], z=obj@fit$g, method = "linear", linear = TRUE, extrap = FALSE)
    gargs <- list(...)
    if(type=="both")
        par(mfrow=c(1,2))
    if(type=="obs" | type=="both"){
        image(intz,...)
        if(is.null(gargs$main)){
            title("Observed values")
        }
        contour(intz,add=TRUE)
    }
    if(type=="g" | type=="both"){
        if(!is.null(obj@fit$g)){
            image(intg,...)
            if(is.null(gargs$main)){
                title("Fitted spline")
            }
            contour(intg,add=TRUE)
        }
    }
    if(type=="fit" | type=="both"){
        image(intm,...)
        if(is.null(gargs$main)){
            title("Fitted model")
        }
        contour(intm,add=TRUE)
    }
}
.scpLevelPlot <- function(obj,what=c("obs","fit","g"),level.at="fivenum",colors=c("yellow","red"),...){
	if(missing(what)){ what <- "fit"} else { if(is.null(what)){ what <- "fit" } }
    what <- match.arg(what)
    requireNamespace("interp",quietly=TRUE)
    if(what=="obs")
        values <- obj@zV
    if(what=="fit")
        values <- obj@fit$fitted.values
    if(what=="g")
        values <- obj@fit$g
    int <- interp::interp(x=obj@data@coords[,1], y=obj@data@coords[,2], z=values, method = "linear", linear = TRUE, extrap = FALSE)
    igridxy <- (cbind(1,int$y)%x%cbind(1,int$x))[,2:3]
    igridx <- igridxy[,1]
    igridy <- igridxy[,2]
    igridz <- as.vector(int$z)
    if(!is.null(level.at)){
        if(is.character(level.at)){
            level.at <- as.vector(apply(cbind(igridz),2,level.at))
        }
    } else {
        level.at <- fivenum(igridz)
    }
    lev.color <- .colorVecMat(level.at, which = "colorRampPalette", col.args = list(colors=colors))
    requireNamespace("lattice",quietly=TRUE)
    lattice::levelplot(igridz~igridx*igridy, at = level.at, col.regions = lev.color, panel = lattice::panel.levelplot.raster, xlab = names(obj@data@coords)[1], ylab = names(obj@data@coords)[2], ...)
}
.scpModelPlot <- function(x, what = c("obs","fit","g"), type = c("levelplot", "image", "persp", "persp3d"), which = "colorRampPalette", col.args = list(colors = c("yellow","red")), col.contour = "black", level.at = "fivenum", border = "transparent", theta = 0, phi = 45, shade = 0.1, ...){
	obj <- x
    if(missing(what)){ what <- "fit" } else { if(is.null(what)) { what <- "fit" }}
    what <- match.arg(what)
    if(missing(type)){ type <- "image" } else { if(is.null(type)) { type <- "image" }}
    type <- match.arg(type)
    if(type=="levelplot"){
        .scpLevelPlot(obj, level.at = level.at, what = what)
    } else {
        requireNamespace("interp", quietly = TRUE)
        requireNamespace("lattice", quietly = TRUE)
        if(what=="obs")
            vals <- obj@zV
        if(what=="fit")
            vals <- obj@fit$fitted.values
        if(what=="g"){
            if(!is.null(obj@fit$g)){
            	vals <- obj@fit$g
            } else {
                return(message("No g function in the model!"))
            }
        }
        int <- interp::interp(x=obj@data@coords[,1], y=obj@data@coords[,2], z=vals, method = "linear", linear = TRUE, extrap = FALSE)
        vv <- int$z
        if(type=="image" | type=="levelplot"){
            vv <- as.vector(vv)[order(vv)]
            cols <- .colorVecMat(vv, which = which, col.args = col.args)
        } else if(type=="persp"){
            cols <- .colorVecMat(vv, which = which, col.args = col.args)
        } else if(type=="persp3d"){
            vv <- as.vector(vv)
            cols <- .colorVecMat(vv, which = which, col.args = col.args)
        }
        cf <- match.call(expand.dots=TRUE)
        cto <- cf[-c(na.omit(match(c("obj","what","type","which","col.args","col.contour","level.at",if(type!="persp" & type!="persp3d"){c("border")},if(type!="persp"){c("theta","phi","shade")}),names(cf))))]
        cto$x <- int$x
        cto$y <- int$y
        cto$z <- int$z
        cto$col <- cols
        if(is.null(cto$xlab))
            cto$xlab <- names(obj@data@coords)[1]
        if(is.null(cto$ylab))
            cto$ylab <- names(obj@data@coords)[2]
        if(is.null(cto$zlab) & type!="image" & type!="levelplot")
            cto$zlab <- "z"
        if(type=="image"){
            cto[[1]] <- quote(image)
            eval(cto)
            cto[[1]] <- quote(contour)
            cto$col <- col.contour
            if(is.null(cto$add))
                cto$add <- TRUE
            eval(cto)
        }
        if(type=="persp"){
            if(all(na.omit(cto$z==mean(cto$z,na.rm=TRUE))))
                cto$zlim <- c(min(cto$z,na.rm=TRUE)-sqrt(.Machine$double.neg.eps),max(cto$z,na.rm=TRUE)+sqrt(.Machine$double.neg.eps))
            cto[[1]] <- quote(persp)
            eval(cto)
        }
        if(type=="persp3d"){
            requireNamespace("rgl",quietly=TRUE)
            cto[[1]] <- quote(rgl::persp3d)
            eval(cto)
        }
    }
}
.getInformationCriterion <- function(object,what=c("m-aic","c-aic","j-bic","gic","hq-gic","b-gic","pn-gic","aic","bic","c-bic"),k=2,tol=.Machine$double.neg.eps){
	obj <- object
	if(!is.null(object@fit$approximate))
		if(object@fit$approximate) stop("Object must be a fit from scp command.")
	if(missing(what)){ what <- "aic"} else { if(is.null(what)){ what <- "aic" } }
    what <- match.arg(what)
    np <- ncol(obj@fit$X)
    nep <- sum(diag(obj@fit$M))
    nq <- length(c(obj@fit$sigma2s,if(!is.null(obj@fit$alpha)){obj@fit$alpha},if(!is.numeric(eval(obj@call$fix.nugget))){obj@fit$tau2},obj@fit$phi,if(!is.numeric(eval(obj@call$fix.kappa))){obj@fit$kappa}))
    Sigma <- obj@fit$sigma2s%x%obj@fit$R
    mu <- obj@fit$fitted.values
    Xb <- mu
	V <- Sigma
    if(!is.null(obj@fit$Zs) & !is.null(obj@fit$r)){
        Xb <- mu - obj@fit$Zs%*%obj@fit$r
        V <- Sigma + obj@fit$sigma2s%x%obj@fit$Zs%*%obj@fit$penalty$decompose$A%*%t(obj@fit$Zs)
	}
    if(what=="m-aic"){#Vaida and Blanchard (2005)
        a0 <- k*length(obj@zV)/(length(obj@zV)-np-nq-1)
        npar <- np + nq
        vmean <- Xb
        mvar <- V
    }
    if(what=="c-aic"){#Burnham and Anderson (2002)
        a0 <- k
        npar <- nep + nq
        vmean <- mu
        mvar <- Sigma
    }
    if(what=="aic"){#Akaike's
        a0 <- k
        npar <- np + nq
        vmean <- Xb
        mvar <- V
    }
    if(what=="bic"){#Schwart's
        a0 <- log(length(obj@zV))
        npar <- np + nq
        vmean <- Xb
        mvar <- V
    }
    if(what=="c-bic"){#Mario's
        a0 <- log(length(obj@zV))
        npar <- nep + nq
        vmean <- mu
        mvar <- Sigma
    }
    if(what=="j-bic"){#Jones (2011)
        vone <- rep(1,length(obj@zV))
        U <- diag(diag(V),nrow(V))
        a0 <- log(t(vone)%*%U%*%solve(V,tol=tol)%*%U%*%vone)
        npar <- np + nq
        vmean <- Xb
        mvar <- V
    }
    if(what=="gic"){#Pu and Niu (2006)
        a0 <- k
        npar <- np + nq
        vmean <- Xb
        mvar <- V
    }
    if(what=="hq-gic"){#Hannan-Quinn (1979)
        a0 <- k*log(log(length(obj@zV)))
        npar <- np + nq
        vmean <- Xb
        mvar <- V
    }
    if(what=="b-gic"){#Bozgodan (1987)
        a0 <- (log(length(obj@zV))+1)
        npar <- np + nq
        vmean <- Xb
        mvar <- V
    }
    if(what=="pn-gic"){#Pu and Niu (2006)
        a0 <- sqrt(log(length(obj@zV)))
        npar <- np + nq
        vmean <- Xb
        mvar <- V
    }
    pen <- a0*npar
    requireNamespace("mvtnorm",quietly=TRUE)
    lL <- mvtnorm::dmvnorm(x=obj@zV,mean=vmean,sigma=mvar,log=TRUE)
    cr <- 2*lL - pen
    list(logLik = as.numeric(lL), criterion = as.numeric(-cr), ka0 = as.numeric(a0), numpar = npar, penalty = as.numeric(pen))
}
setGeneric("scpApproximate",function(object, tol){ standardGeneric("scpApproximate") })
setMethod("scpApproximate", signature(object = "sssFit"), function(object, tol){
	.scpApproxGet(object = object, tol = tol)
})
setGeneric("testSurface",function(object, tol){ standardGeneric("testSurface") })
setMethod("testSurface", signature(object = "sssFit"), function(object, tol){
	.scpTestGet(object = object, tol = tol)
})
setGeneric("plot")
setGeneric("summary")
setGeneric("AIC")
setGeneric("BIC")
setGeneric("AICm",function(object, k, only.criterion){ standardGeneric("AICm") })
setGeneric("AICc",function(object, k, only.criterion){ standardGeneric("AICc") })
setGeneric("BICc",function(object, only.criterion){ standardGeneric("BICc") })
setGeneric("BICj",function(object, k, tol, only.criterion){ standardGeneric("BICj") })
setGeneric("GIC",function(object, k, only.criterion){ standardGeneric("GIC") })
setGeneric("GIChq",function(object, k, only.criterion){ standardGeneric("GIChq") })
setGeneric("GICpn",function(object, only.criterion){ standardGeneric("GICpn") })
setGeneric("GICb",function(object, only.criterion){ standardGeneric("GICb") })
setGeneric("Variogram",function(object, distance, plot, ...){ standardGeneric("Variogram") })
setMethod("plot", signature(x = "sssFit", y = "missing"), function(x, what, type, which, col.args, col.contour, level.at, border, theta, phi, shade, ...){
	if(missing(what)){ what <- "fit" } else {if(is.null(what)){ what <- "fit" }}
	if(missing(type)){ type <- "image" } else {if(is.null(type)){ type <- "image" }}
	if(missing(which)){ which <- "colorRampPalette" } else {if(is.null(which)){ which <- "colorRampPalette" }}
	if(missing(col.args)){ col.args <- list(colors=c("yellow","red")) } else {if(is.null(col.args)){ col.args <- list(colors=c("yellow","red")) }}
	if(missing(col.contour)){ col.contour <- "black" } else {if(is.null(col.contour)){ col.contour <- "black" }}
	if(missing(level.at)){ level.at <- "fivenum" } else {if(is.null(level.at)){ level.at <- "fivenum" }}
	if(missing(border)){ border <- "transparent" } else {if(is.null(border)){ border <- "transparent" }}
	if(missing(theta)){ theta <- 0 } else {if(is.null(theta)){ theta <- 0 }}
	if(missing(phi)){ phi <- 45 } else {if(is.null(phi)){ phi <- 45 }}
	if(missing(shade)){ shade <- 0.1 } else {if(is.null(shade)){ shade <- 0.1 }}
	.scpModelPlot(x, what, type, which, col.args, col.contour, level.at, border, theta, phi, shade, ...)
})
setMethod("summary", signature(object = "sssFit"), .scpSummary)
setMethod("AIC", signature(object = "sssFit", k = "ANY"), function(object, k, only.criterion){
	if(missing(k)){ k = 2 } else { if(is.null(k)){ k = 2 } }
	if(missing(only.criterion)){ only.criterion = TRUE } else { if(is.null(only.criterion)){ only.criterion = TRUE } }
	if(only.criterion){
		.getInformationCriterion(object,what="aic",k=k)$criterion
	} else {
		.getInformationCriterion(object,what="aic",k=k)
	}
})
setMethod("BIC", signature(object = "sssFit"), function(object, only.criterion){
	if(missing(only.criterion)){ only.criterion = TRUE } else { if(is.null(only.criterion)){ only.criterion = TRUE } }
	if(only.criterion){
		.getInformationCriterion(object,what="bic")$criterion
	} else {
		.getInformationCriterion(object,what="bic")
	}
})
setMethod("AICm", signature(object = "sssFit", k = "ANY"), function(object, k, only.criterion){
	if(missing(k)){ k = 2 } else { if(is.null(k)){ k = 2 } }
	if(missing(only.criterion)){ only.criterion = TRUE } else { if(is.null(only.criterion)){ only.criterion = TRUE } }
	if(only.criterion){
		.getInformationCriterion(object,what="m-aic",k=k)$criterion
	} else {
		.getInformationCriterion(object,what="m-aic",k=k)
	}
})
setMethod("AICc", signature(object = "sssFit", k = "ANY"), function(object, k, only.criterion){
	if(missing(k)){ k = 2 } else { if(is.null(k)){ k = 2 } }
	if(missing(only.criterion)){ only.criterion = TRUE } else { if(is.null(only.criterion)){ only.criterion = TRUE } }
	if(only.criterion){
		.getInformationCriterion(object,what="c-aic",k=k)$criterion
	} else {
		.getInformationCriterion(object,what="c-aic",k=k)
	}
})
setMethod("BICc", signature(object = "sssFit"), function(object, only.criterion){
	if(missing(only.criterion)){ only.criterion = TRUE } else { if(is.null(only.criterion)){ only.criterion = TRUE } }
	if(only.criterion){
		.getInformationCriterion(object,what="c-bic")$criterion
	} else {
		.getInformationCriterion(object,what="c-bic")
	}
})
setMethod("BICj", signature(object = "sssFit", k = "ANY", tol = "ANY"), function(object, k, tol, only.criterion){
	if(missing(k)){ k = 2 } else { if(is.null(k)){ k = 2 } }
	if(missing(tol)){ tol = .Machine$double.neg.eps } else { if(is.null(tol)){ tol = .Machine$double.neg.eps } }
	if(missing(only.criterion)){ only.criterion = TRUE } else { if(is.null(only.criterion)){ only.criterion = TRUE } }
	if(only.criterion){
		.getInformationCriterion(object,what="j-bic",k=k,tol=tol)$criterion
	} else {
		.getInformationCriterion(object,what="j-bic",k=k,tol=tol)
	}
})
setMethod("GIC", signature(object = "sssFit", k = "ANY"), function(object, k, only.criterion){
	if(missing(k)){ k = 2 } else { if(is.null(k)){ k = 2 } }
	if(missing(only.criterion)){ only.criterion = TRUE } else { if(is.null(only.criterion)){ only.criterion = TRUE } }
	if(only.criterion){
		.getInformationCriterion(object,what="gic",k=k)$criterion
	} else {
		.getInformationCriterion(object,what="gic",k=k)
	}
})
setMethod("GIChq", signature(object = "sssFit", k = "ANY"), function(object, k, only.criterion){
	if(missing(k)){ k = 2 } else { if(is.null(k)){ k = 2 } }
	if(missing(only.criterion)){ only.criterion = TRUE } else { if(is.null(only.criterion)){ only.criterion = TRUE } }
	if(only.criterion){
		.getInformationCriterion(object,what="hq-gic",k=k)$criterion
	} else {
		.getInformationCriterion(object,what="hq-gic",k=k)
	}
})
setMethod("GICpn", signature(object = "sssFit"), function(object, only.criterion){
	if(missing(only.criterion)){ only.criterion = TRUE } else { if(is.null(only.criterion)){ only.criterion = TRUE } }
	if(only.criterion){
		.getInformationCriterion(object,what="pn-gic")$criterion
	} else {
		.getInformationCriterion(object,what="pn-gic")
	}
})
setMethod("GICb", signature(object = "sssFit"), function(object, only.criterion){
	if(missing(only.criterion)){ only.criterion = TRUE } else { if(is.null(only.criterion)){ only.criterion = TRUE } }
	if(only.criterion){
		.getInformationCriterion(object,what="b-gic")$criterion
	} else {
		.getInformationCriterion(object,what="b-gic")
	}
})
setMethod("Variogram", signature(object = "sssFit", distance = "ANY"),
	function(object, distance, plot, ...){
	if(missing(distance)){ distance <- NULL } else { if(!is.null(distance)){ if(!is.numeric(distance)){
	warning("Variogram: distance must be numeric or null. Setting to NULL.")
	distance <- NULL
	} } }
	if(missing(plot)){ plot <- TRUE } else { if(is.null(TRUE)){ plot <- TRUE } }
	.plotSemivariogram(object,h=distance,hmax=NULL,size=10000,plot=plot,prange=FALSE,legend=FALSE, ...)
	})
.tryCatchWE = function (expr){
    W <- NULL
    w.handler <- function(w) {
        W <<- w
        invokeRestart("muffleWarning")
    }
    value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler)
    isE = !is.na(pmatch("Error",paste(value,collapse="")))
    list(value = value, error = isE, warning = W)
}
.defaultArg = function(arg,default){
 env = parent.frame()
 if(!is.character(arg)) stop("arg must be character!")
 txt = NULL
 txt[1] = paste("if(missing(",arg,") || is.null(",arg,")){",arg," = ",shQuote(default),"}",sep="")
 txt[2] = paste("if(!is.character(",arg,")){",arg," = as.character(",arg,")}",sep="")
 txt[3] = paste(arg," = match.arg(",arg,")",sep="")
 txt[4] = paste("if(!is.character(",default,")){",arg," = as.numeric(",arg,")}",sep="")
 txt = txt[c(1:ifelse(is.character(default),3,4))]
 eval(parse(text=txt),envir=env)
}
.invSchur = function(x, Re.values = FALSE, ...){
 requireNamespace("Matrix", quietly = TRUE)
 ex = Matrix::Schur(x, ...)
 if(any(is.complex(ex$EValues))){
  warning("Complex eigenvalues!")
 }
 if(Re.values){
  d.v = Re(ex$EValues)
  ex.Q = Re(ex$Q)
 } else {
  d.v = ex$EValues
  ex.Q = ex$Q
 }
 ex.Q%*%diag(1/d.v)%*%t(ex.Q)
}
.invEigen = function(x, Re.values = FALSE, ...){
 ex = eigen(x, ...)
 if(Re.values){
  d.v = Re(ex$values)
  ex.v = Re(ex$vectors)
 } else {
  d.v = ex$values
  ex.v = ex$vectors
 }
 ex.v%*%diag(1/d.v)%*%t(ex.v)
}
.invGsvd = function(x, tol = sqrt(.Machine$double.eps), ...){
 ex = svd(x, ...)
 r.s = ex$d > tol*ex$d[1]
 if(any(r.s)){
  ex$v[,r.s]%*%(t(ex$u[,r.s])/ex$d[r.s])
 } else {x}
}
.whichInverse = function(name = c("solve","ginv",".invEigen",".invSchur",".invGsvd"), ...){
	.defaultArg("name","solve")
	if(name=="ginv"){
		requireNamespace("MASS", quietly = TRUE)
		name = paste("MASS::",name,sep="")
	}
	l.args = list(...)
	o.args = formals(get("name"))
	c.args = l.args[-1]
	v.args = match(names(c.args),names(o.args))
	if(any(!is.na(v.args))){
		l.args = l.args[c(TRUE,!is.na(v.args))]
	} else {
		l.args = l.args[1]
	}
	do.call(name, l.args)
}
.transfTo01 <- function(y,dom=NULL,to.Inf=.Machine$double.base**(-.Machine$double.min.exp)){
 dom.obs <- c(min(y),max(y))
 if(is.null(dom)){
  dom <- dom.obs
 }
 lim.a <- ifelse(min(dom) <= -to.Inf, -to.Inf, min(dom))
 lim.b <- ifelse(max(dom) >= +to.Inf, +to.Inf, max(dom))
 y.01 <- (y - lim.a)/(lim.b-lim.a)
 return(y.01)
}
.transfFrom01 <- function(y,dom,to.Inf=.Machine$double.base**(-.Machine$double.min.exp)){
	if(missing(dom))
		stop("domain of the original values must be specified (otherwise infinity transformations are available)")
	lim.a <- ifelse(min(dom) <= -to.Inf, -to.Inf, min(dom))
	lim.b <- ifelse(max(dom) >= +to.Inf, +to.Inf, max(dom))
	y.or <- y*(lim.b-lim.a) + lim.a
	return(y.or)
}
.getKnots <- function(x){
	requireNamespace("stats", quietly = TRUE)
    keep.class = class(x)
    aux = stats::xtabs(~x)
    knots.values <- names(aux)
    class(knots.values) <- keep.class
    return(knots.values)
}
.tpsBase = function(x,...){
 s.diff = s.norm(coord=x,...)
 l.s = lower.tri(s.diff,diag=FALSE)
 sl = s.diff[l.s]
 el = (1/(16*pi))*sl*log(sl)
 d.s = !(l.s + t(l.s))
 Em = matrix(0,nrow=nrow(s.diff),ncol=ncol(s.diff))
 Em[l.s] = el
 Em[t(l.s)] = el
 Em
}
.tpsBase = function(x,tol=.Machine$double.neg.eps,...){
    eta.r <- function(r) ifelse(r>0,(1/(16*pi))*r**2*log(r**2),0)
    dist.titj <- .h2D(as.matrix(x),...)$h
    aux.Q <- apply(dist.titj,2,eta.r)
    list(V=diag(1,nrow(aux.Q)), A=solve(aux.Q,tol=tol), Ai=aux.Q, Asr=diag(1,nrow(aux.Q)), K=aux.Q, type="tps")
}
.csBase = function (x){
    knots.values = .getKnots(x)
    nv = length(x)
    kv = length(knots.values)
    hv = knots.values[2:kv] - knots.values[1:(kv - 1)]
    Vm = matrix(0, nrow = kv, ncol = kv - 2)
    Am = matrix(0, nrow = kv - 2, ncol = kv - 2)
    for (k in 1:(kv - 2)) {
        Vm[k, k] = 1/(hv[k])
        Vm[k + 1, k] = -(1/(hv[k]) + 1/(hv[k + 1]))
        Vm[k + 2, k] = 1/(hv[k + 1])
    }
    for (k in 1:(kv - 2)) {
        if (k < kv - 2) {
            Am[k, k + 1] = hv[k + 1]/6
            Am[k + 1, k] = hv[k + 1]/6
            Am[k, k] = (hv[k] + hv[k + 1])/3
        }
        if (k == (kv - 2)) {
            Am[k, k] = (hv[k] + hv[k + 1])/3
        }
    }
    KKm = Vm %*% qr.solve(Am) %*% t(Vm)
    A.im.sr = svd(Am)
    A.im.sr = A.im.sr$v %*% diag(sqrt(A.im.sr$d)) %*% t(A.im.sr$v)
    list(V = Vm, A = Am, K = KKm, Asr = A.im.sr)
}
.psBase = function(x, diff.order = 1:1000000){
 .defaultArg("diff.order",2)
 x.tau = .getKnots(x)
 aux.D = diff(diag(length(x.tau)), differences = diff.order)
 aux.Q = t(aux.D)%*%aux.D
 list(V=t(aux.D), A=diag(1,nrow(aux.D)), Asr=diag(1,nrow(aux.D)), K=aux.Q)
}
.tensorBase = function(obj=NULL, grid.list=NULL, type=c("cs","ps"), diff.order=2){
 .defaultArg("type","cs")
 if(is.null(obj) && is.null(grid.list)) stop("obj or grid.list must be specified!")
 if(!is.null(grid.list) && !is.matrix(grid.list) && is.null(obj)) grid.list = grid.list
 if(!is.null(grid.list) && is.matrix(grid.list) && is.null(obj)) grid.list = attributes(grid.list)$grid.list
 if(!is.null(obj) && is.sss(obj)) grid.list = attributes(obj@grid)$grid.list
 if(is.null(grid.list) && !is.null(obj) && is.sss(obj)) grid.list = obj@knots
 II.x = diag(1,length(grid.list[[1]]))
 JJ.y = diag(1,length(grid.list[[2]]))
 if(type=="cs"){
  QQ.x = .csBase(grid.list[[1]])
  QQ.y = .csBase(grid.list[[2]])
 } else {
  QQ.x = .psBase(grid.list[[1]],diff.order)
  QQ.y = .psBase(grid.list[[2]],diff.order)
 }
 list(Ix=II.x,Qx=QQ.x,Jy=JJ.y,Qy=QQ.y,type=type,difference=if(type=="ps"){diff.order})
}
.tensorPenalty = function(obj, alpha = c(1,1), type=c("Q","P","both")){
 .defaultArg("type","P")
 if(missing(obj) || is.null(obj)) stop("tensor base obj must be specified!")
 if(length(alpha)!=2L || !is.numeric(alpha)) stop("alpha must be a 2*1 vector of smoothing parameters for the coordinates")
 if(type=="Q") return(obj$Jy%x%obj$Qx$K + obj$Qy$K%x%obj$Ix)
 if(type=="P") return(alpha[1]*(obj$Jy%x%obj$Qx$K) + alpha[2]*(obj$Qy$K%x%obj$Ix))
 if(type=="both") list(Qxy = obj$Jy%x%obj$Qx$K + obj$Qy$K%x%obj$Ix, Pxy = alpha[1]*(obj$Jy%x%obj$Qx$K) + alpha[2]*(obj$Qy$K%x%obj$Ix))
}
.decomposePenalty = function(PM, method.eigen = TRUE, diagA = FALSE, thr = .Machine$double.eps, ftr = ifelse(method.eigen,1,10), inverse.method = c("solve","ginv",".invEigen",".invSchur",".invGsvd"), itol = .Machine$double.neg.eps){
 .defaultArg("inverse.method","solve")
 if(missing(PM)) stop("PM must be a penalty matrix!")
 PM = as.matrix(PM)
 if(nrow(PM)!=ncol(PM)) stop("PM must be a squared matrix!")
 d.PM = do.call(ifelse(method.eigen,"eigen","svd"),list(x=PM))
 d.rule = d.PM[[ifelse(method.eigen,"values","d")]] > sqrt(thr)*max(d.PM[[ifelse(method.eigen,"values","d")]])*ftr
 if(!diagA){
  V.el = d.PM[[ifelse(method.eigen,"vectors","u")]][,d.rule]
  Ai.el = diag(d.PM[[ifelse(method.eigen,"values","d")]][d.rule])
  A.el = .whichInverse(inverse.method, Ai.el, tol = itol, Re.values = TRUE)
  A.sr.el = sqrt(A.el)
 } else {
  V.el = (d.PM[[ifelse(method.eigen,"vectors","u")]][,d.rule])%*%sqrt(diag(d.PM[[ifelse(method.eigen,"values","d")]][d.rule]))
  Ai.el = diag(1,sum(d.rule))
  A.el = Ai.el
  A.sr.el = Ai.el
 }
 list(V = V.el, A = A.el, Ai = Ai.el, Asr = A.sr.el)
}
.generateXZ <- function(coords, QList){
	if(missing(QList)) QList <- NULL
	if(is.null(QList)) stop("QList must be specified!")
	X1 <- cbind(1,.getKnots(coords[,1]))
	X2 <- cbind(1,.getKnots(coords[,2]))
	Z1 <- QList$Qx$V%*%solve(t(QList$Qx$V)%*%QList$Qx$V)
	Z2 <- QList$Qy$V%*%solve(t(QList$Qy$V)%*%QList$Qy$V)
	X <- cbind(X2%x%X1)
	Z <- cbind(X2%x%Z1,Z2%x%X1,Z2%x%Z1)
	list(X = X, Z = Z)
}
.generateV <- function(Z, itol = .Machine$double.neg.eps) return(Z%*%solve(t(Z)%*%Z, tol = itol))
.generateAi <- function(Z, Q, itol = .Machine$double.neg.eps) return(t(Z)%*%Q%*%Z)
.generateZs <- function(Rs, V, itol = .Machine$double.neg.eps) return(Rs%*%V%*%solve(t(V)%*%Rs%*%V, tol = itol))
.incidenceMatrix = function(x, by = NULL){
    knots.values = .getKnots(x)
    if(is.null(by)){
     nv = length(x)
    } else {
     knots.i = by(x,by,.getKnots,simplify=TRUE)
     nprod = 1
     for(j in 1:length(by)){
      nprod = nprod*nlevels(factor(by[[j]]))
     }
     nv = nprod*length(knots.values)
     NAknots = rep(NA,length(knots.values))
     NAn = NULL
     for(j in 1:length(knots.i)){
      NAn[[j]] = NAknots
      NAn[[j]][match(knots.i[[j]],knots.values)] = knots.i[[j]]
     }
     NAn = unlist(NAn)
     x = NAn
    }
    kv = length(knots.values)
    W = matrix(0, nrow = nv, ncol = kv)
    xc = paste(round(x, options()$digits))
    kc = paste(round(knots.values, options()$digits))
    for (i in 1:nv) {
        for (r in 1:kv) {
            if (xc[i] == kc[r]) {
                W[i, r] = 1
            }
        }
    }
    return(W)
}
.incidenceMatrix = function(x, by = NULL, add = FALSE){
    knots.values = .getKnots(x)
    if(is.null(by)){
     nv = length(x)
    } else {
     knots.i = by(x,by,.getKnots,simplify=TRUE)
     nprod = 1
     for(j in 1:length(by)){
      nprod = nprod*nlevels(factor(by[[j]]))
     }
     nv = nprod*length(knots.values)
     NAknots = rep(NA,length(knots.values))
     NAn = NULL
     for(j in 1:length(knots.i)){
      NAn[[j]] = NAknots
      NAn[[j]][match(knots.i[[j]],knots.values)] = knots.i[[j]]
     }
     NAn = unlist(NAn)
     x = NAn
    }
    kv = length(knots.values)
    W = matrix(0, nrow = nv, ncol = kv - ifelse(add, 1, 0))
    xc = paste(round(x, options()$digits))
    kc = paste(round(knots.values, options()$digits))
    for (i in 1:nv) {
        for (r in 1:(kv - ifelse(add, 1, 0))) {
            if(add){
                if (xc[i] == kc[r] && xc[i] != kc[kv]) {
                    W[i, r] = 1
                } else if (xc[i] == kc[kv]) {
                    W[i, r] = -1
                }
            } else {
                if (xc[i] == kc[r]) {
                    W[i, r] = 1
                }
            }
        }
    }
    return(W)
}
.incidenceSpatial <- function(coords,grid,obs){
 if(missing(grid)) grid <- NULL
 if(missing(obs)) obs <- NULL
 if(is.null(coords)) stop(".incidenceSpatial: coordinates must be defined!")
 coords.id = apply(coords,1,paste,sep="",collapse=":")
 if(is.null(grid)) grid <- .createGrid(coords)
 grid.id = apply(grid,1,paste,sep="",collapse=":")
 if(is.null(grid.knots <- attributes(grid)$grid.list)){
  n.x = length(.getKnots(coords[,1]))
  n.y = length(.getKnots(coords[,2]))
 } else {
  n.x = length(grid.knots[[1]])
  n.y = length(grid.knots[[2]])
 }
 if(is.null(obs)) obs <- rep(1,length(coords.id))
 row.id = match(coords.id,grid.id)
 W.coords = matrix(0, nrow = n.x*n.y, ncol = length(coords.id))
 for(i in 1:length(row.id)){
  W.coords[row.id[i],i] = ifelse(!is.na(obs[i]),1,0)
 }
 W.coords <- t(W.coords)
 attr.Wxy = attributes(W.coords)
 attributes(W.coords) = attr.Wxy
 return(W.coords)
}
.incidenceSpatial <- function(coords,grid,obs,add=FALSE){
    if(missing(grid)) grid <- NULL
    if(missing(obs)) obs <- NULL
    if(is.null(coords)) stop(".incidenceSpatial: coordinates must be defined!")
    coords.id = apply(coords,1,paste,sep="",collapse=":")
    if(is.null(grid)) grid <- .createGrid(coords)
    grid.id = apply(grid,1,paste,sep="",collapse=":")
    if(is.null(grid.knots <- attributes(grid)$grid.list)){
        n.x = length(.getKnots(coords[,1]))
        n.y = length(.getKnots(coords[,2]))
    } else {
        n.x = length(grid.knots[[1]])
        n.y = length(grid.knots[[2]])
    }
    if(is.null(obs)) obs <- rep(1,length(coords.id))
    row.id = match(coords.id,grid.id)
    W.coords = matrix(0, nrow = n.x*n.y, ncol = length(coords.id))
    for(i in 1:length(row.id)){
        W.coords[row.id[i],i] = ifelse(!is.na(obs[i]),1,0)
    }
    W.coords <- t(W.coords)
    if(add){
        W.i <- cbind(.incidenceMatrix(1:length(coords.id),add=TRUE),0)
        W.coords = W.i%*%W.coords
    }
    attr.Wxy = attributes(W.coords)
    attributes(W.coords) = attr.Wxy
    return(W.coords)
}
.genZ = function(V, R = NULL, Asr = NULL, Ztr = c("V","VA","RV","RVA"), only.z = FALSE,
                 inverse.method = c("solve","ginv",".invEigen",".invSchur",".invGsvd"),
                 itol = .Machine$double.neg.eps){
 .defaultArg("inverse.method","solve")
 if(missing(V)) stop("Matrix V must be specified!")
 .defaultArg("Ztr","RV")
 if(Ztr=="RV" | Ztr=="RVA"){
  VRV = t(V)%*%R%*%V
  VRVi = .whichInverse(inverse.method, VRV, tol=itol, Re.values=TRUE)
  Z = R%*%V%*%VRVi
 }
 if(Ztr=="RVA") Z = Z%*%Asr
 if(Ztr=="V" | Ztr=="VA"){
  VRV = t(V)%*%V
  VRVi = .whichInverse(inverse.method, VRV, tol=itol, Re.values=TRUE)
  Z = V%*%VRVi
 }
 if(Ztr=="VA") Z = Z%*%Asr
 if(only.z) return(Z) else list(Z=Z,VRV=VRV,VRVi=VRVi)
}
.autoCor <- function(model = c("identity","pure.nugget","matern","powered.exponential","spherical","wave","exponential","gaussian","cubic","circular","gencauchy","cauchy","RMmatern","RMwhittle","RMgneiting","RMgengneiting","RMnugget","RMcauchy","RMexp","RMgencauchy","RMgauss","RMspheric","RMstable","RMpoweredexp"), what = c("expr","bounds","all"), ...){
    m.c <- match.call(expand.dots=TRUE)
    .defaultArg("model","RMwhittle")
    .defaultArg("what","expr")
    model.RF = c("RMmatern","RMwhittle","RMgneiting","RMgengneiting","RMnugget","RMcauchy","RMexp","RMgencauchy","RMgauss","RMspheric","RMstable","RMpoweredexp")
    el = list()
    if(what == "expr" | what == "all"){
        add.args = list(...)
        if(is.null(add.args$tol.nugget)) tol.nugget <- 1.0e-15 else tol.nugget <- add.args$tol.nugget
        if(model == "matern") el$rho.h = expression(ifelse(h==0,1,((2**(1-kappa))/gamma(kappa))*((h/phi)**kappa)*besselK(h/phi, nu = kappa, expon.scaled = FALSE)))
        if(model == "powered.exponential") el$rho.h = expression(exp(-(h/phi)**kappa))
        if(model == "spherical") el$rho.h = expression(ifelse(h >= 0 & h <= phi,1-(3/2)*(h/phi)+(1/2)*((h/phi)**3),0))
        if(model == "wave") el$rho.h = expression(ifelse(h>0,((h/phi)**(-1))*sin(h/phi),1))
        if(model == "cauchy") el$rho.h = expression(ifelse(h>0,(1+(h/phi)**2)**(-kappa),1))
        if(model == "gencauchy") el$rho.h = expression(ifelse(h>0,(1+(h/phi)**kappa[2])**(-kappa[1]/kappa[2]),1))
        if(model == "circular") el$rho.h = expression(ifelse(h >= 0 & h < phi,1-2*(min(h/phi,1)*sqrt(1-min(h/phi,1)**2)+asin(sqrt(h/phi)))/pi,0))
        if(model == "cubic") el$rho.h = expression(ifelse(h >= 0 & h < phi,1-(7*(h/phi)**2-8.75*(h/phi)**3+3.5*(h/phi)**5-0.75*(h/phi)**7),0))
        if(model == "gaussian") el$rho.h = expression(exp(-(h/phi)**2))
        if(model == "exponential") el$rho.h = expression(exp(-(h/phi)))
        if(model == "identity") el$rho.h = expression(ifelse(h==0,1,0))
        if(model == "pure.nugget") el$rho.h = expression(ifelse(h>=0 & h<tol.nugget,1,0))
        if(!is.na(match(model,model.RF))){
            requireNamespace("RandomFields", quietly = TRUE)
            l.args = formals(get(model))
            add.args = add.args[!is.na(match(names(add.args),names(l.args)))]
            if(any(grepl("Aniso",names(add.args)))) add.args[[grepl("Aniso",names(add.args))]] <- m.c[["Aniso"]]
            if(any(grepl("proj",names(add.args)))) add.args[[grepl("proj",names(add.args))]] <- m.c[["proj"]]
            plug.r = ifelse(length(add.args)>0,paste(",",paste(names(add.args),add.args,sep="=",collapse=","),sep=""),"")
            plug.f = ifelse(length(add.args)>0,paste(",",paste(names(add.args),names(add.args),sep="=",collapse=","),sep=""),"")
            r.args = paste("phi",ifelse(any(!is.na(match(model,c("RMmatern","RMwhittle","RMgengneiting","RMcauchy","RMgencauchy","RMstable","RMpoweredexp")))),",kappa",""),plug.r,sep="")
            f.args = paste("var=1,scale=phi",
                    ifelse(any(!is.na(match(model,c("RMmatern","RMwhittle")))),",nu=kappa",
                           ifelse(any(!is.na(match(model,c("RMgengneiting")))),",mu=kappa[1],kappa=kappa[2]",
                                  ifelse(any(!is.na(match(model,c("RMcauchy")))),",gamma=kappa",
                                         ifelse(any(!is.na(match(model,c("RMgencauchy")))),",beta=kappa[1],alpha=kappa[2]",
                                                ifelse(any(!is.na(match(model,c("RMstable","RMpoweredexp")))),",alpha=kappa",""))))),
                    plug.f,sep="")
        }
        if(any(!is.na(match(model,model.RF)))){
            el$rho.h = eval(parse(text=paste("function(",r.args,"){",paste("RandomFields::",model,sep=""),"(",f.args,")}",sep="")))
            environment(el$rho.h) <- parent.frame()
        }
    }
    if(what == "bounds" | what == "all"){
        if(!is.na(match(model,c("matern","cauchy","RMmatern","RMwhittle","RMcauchy")))){
            el$lower = c(phi=0,kappa=0)
            el$upper = c(phi=Inf,kappa=Inf)
        }
        if(!is.na(match(model,c("RMgengneiting")))){
            el$lower = c(phi=0,kappa1=1,kappa2=0)
            el$upper = c(phi=Inf,kappa1=Inf,kappa2=3)
        }
        if(!is.na(match(model,c("gencauchy","RMgencauchy")))){
            el$lower = c(phi=0,kappa1=0,kappa2=0)
            el$upper = c(phi=Inf,kappa1=Inf,kappa2=2)
        }
        if(!is.na(match(model,c("powered.exponential","RMstable","RMpoweredexp")))){
            el$lower = c(phi=0,kappa=0)
            el$upper = c(phi=Inf,kappa=2)
        }
        if(!is.na(match(model,c("identity","pure.nugget","spherical","wave","exponential","gaussian","cubic","circular","RMnugget","RMgneiting","RMexp","RMgauss","RMspheric")))){#"RMnugget"?
            el$lower = c(phi=0)
            el$upper = c(phi=Inf)
        }
    }
    if(what == "expr"){
        return(el$rho.h)
    } else {
        return(el)
    }
}
.semiVar = function(model = c("matern","gaussian","exponential","power","cubic","penta.spherical","spherical","wave","sin.hole","identity","pure.nugget"), what = c("expr","bounds","all"), tol.nugget = 1.0e-15){
 .defaultArg("model","matern")
 .defaultArg("what","expr")
 sl = list()
 if(what == "expr" | what == "all"){
  if(model == "matern") sl$gamma.h = expression(ifelse(h > tol.nugget,nugget+sill*(1-(2**(1-kappa))/gamma(kappa))*((h/phi)**kappa)*besselK(h/phi, nu = kappa, expon.scaled = FALSE),nugget+0))
  if(model == "gaussian") sl$gamma.h = expression(ifelse(h > tol.nugget,nugget+sill*(1-exp(-(h/phi)**2)),nugget+0))
  if(model == "exponential") sl$gamma.h = expression(ifelse(h > tol.nugget,nugget+sill*(1-exp(-(h/phi))),nugget+0))
  if(model == "power") sl$gamma.h = expression(ifelse(h > tol.nugget,nugget+sill*(h**phi),nugget+0))
  if(model == "cubic") sl$gamma.h = expression(ifelse((h > tol.nugget & h <= phi),nugget+sill*(7*((h/phi)**2)-(35/4)*((h/phi)**3)+(7/2)*((h/phi)**5)-(3/4)*((h/phi)**7)),ifelse(tol.nugget<phi & phi<h,nugget+sill,nugget+0)))
  if(model == "penta.spherical") sl$gamma.h = expression(ifelse((h > tol.nugget & h <= phi),nugget+sill*((15/8)*(h/phi)-(5/4)*((h/phi)**3)+(3/8)*((h/phi)**5)),ifelse(tol.nugget<phi & phi < h,nugget+sill,nugget+0)))
  if(model == "spherical") sl$gamma.h = expression(ifelse((h > tol.nugget & h <= phi),nugget+sill*((3/2)*(h/phi)+(1/2)*((h/phi)**3)),ifelse(tol.nugget<phi & phi < h,nugget+sill,nugget+0)))
  if(model == "wave") sl$gamma.h = expression(ifelse(h > tol.nugget,nugget+sill*(1-((h/phi)**(-1))*sin(h/phi)),nugget+0))
  if(model == "sin.hole") sl$gamma.h = expression(ifelse(h > tol.nugget,nugget+sill*(1-((pi*h/phi)**(-1))*sin(pi*h/phi)),nugget+0))
  if(model == "pure.nugget") sl$gamma.h = expression(ifelse(h > tol.nugget,nugget+0,0))
  if(model == "identity") sl$gamma.h = expression(ifelse(h > tol.nugget,nugget+sill,nugget+0))
 }
 if(what == "bounds" | what == "all"){
  if(model == "matern"){
   sl$lower = c(sill=0, phi=0, kappa=0, nugget=0)
   sl$upper = c(sill=Inf, phi=Inf, kappa=Inf, nugget=Inf)
  }
  if(!is.na(match(model, c("gaussian","exponential","power","cubic","penta.spherical","spherical","wave","sin.hole")))){
   sl$lower = c(sill=0, phi=0, nugget=0)
   sl$upper = c(sill=Inf, phi=Inf, nugget=Inf)
  }
  if(model == "identity"){
   sl$lower = c(sill=0, nugget=0)
   sl$upper = c(sill=Inf, nugget=Inf)
  }
  if(model == "pure.nugget"){
   sl$lower = c(nugget=0)
   sl$upper = c(nugget=Inf)
  }
 }
 if(what == "expr"){
  return(sl$gamma.h)
 } else {
  return(sl)
 }
}
.autoCorF = function(h, phi = 1, kappa = NULL, model = c("matern","identity","pure.nugget","powered.exponential","spherical","wave","exponential","gaussian","cubic","circular","gencauchy","cauchy","RMmatern","RMwhittle","RMgneiting","RMgengneiting","RMnugget","RMcauchy","RMexp","RMgencauchy","RMgauss","RMspheric","RMstable","RMpoweredexp"), ...){
    m.c <- match.call(expand.dots=TRUE)
    if(missing(phi) | is.null(phi)){phi = 1}
    .defaultArg("model","RMwhittle")
    add.args = list(...)
    if(is.null(add.args$tol.nugget)) tol.nugget <- 1.0e-15 else tol.nugget <- add.args$tol.nugget
    model.RF = c("RMmatern","RMwhittle","RMgneiting","RMgengneiting","RMnugget","RMcauchy","RMexp","RMgencauchy","RMgauss","RMspheric","RMstable","RMpoweredexp")
    if(!is.na(match(model,model.RF))){
        requireNamespace("RandomFields", quietly = TRUE)
    }
    rho = .autoCor(model = model, what = "all", ...)
    if(any(grepl("kappa",names(rho$lower))) & is.null(kappa)){
        stop("Parameter kappa was not defined!")
    }
    theta = c(phi=phi,kappa=kappa)
    for(i in 1:length(rho$lower)){
        if(!(theta[i]>=rho$lower[i] && theta[i]<=rho$upper[i])){
            stop(paste("Parameter ",names(theta)[i]," lies outside allowed region for its family!\nEnter ",names(theta)[i]," value in ",ifelse(rho$lower[i]==0,"(","["),rho$lower[i],",",rho$upper[i],ifelse(rho$upper[i]==Inf,").","]."),sep="",collapse=""))
        }
    }
    if(!any(grepl("kappa",names(rho$lower))) & !is.null(kappa)){
        message("kappa not needed for the family, so not used!")
    }
    h.mat = (is.matrix(h) || is.data.frame(h)) && nrow(h)==ncol(h)
    if(h.mat){
        if(!isSymmetric(h) || nrow(h)!=ncol(h)) stop("h must be a vector or a symmetric n*n matrix of distances!")
        Rho.h = matrix(1,nrow=nrow(h),ncol=ncol(h))
        low.h = lower.tri(Rho.h,diag=FALSE)
        h = h[low.h]
    }
    if(any(grepl(model,c("matern","identity","pure.nugget","powered.exponential","spherical","wave","exponential","gaussian","cubic","circular","gencauchy","cauchy")))){
        aux = eval(rho$rho.h)
    } else {
        call.cf = eval(rho$rho.h)
        fargs <- formals(call.cf)
        nargs <- names(formals(call.cf))
        largs <- list()
        for(i in 1:length(nargs)){
            if(any(grepl(nargs[i],names(theta)))){
                largs[[i]] <- theta[grepl(nargs[i],names(theta))]
            } else {
                largs[[i]] <- fargs[[i]]
            }
            names(largs)[i] <- nargs[i]
        }
        d.c <- call("call.cf")
        for(i in 1:length(largs)){
            if(names(largs)[i]=="Aniso"){
                d.c[[names(largs)[i]]] <- m.c[["Aniso"]]
            } else if(names(largs)[i]=="proj") {
                d.c[[names(largs)[i]]] <- m.c[["proj"]]
            } else {
                d.c[[names(largs)[i]]] <- largs[[i]]
            }
        }
        e.c <- eval(d.c)
        environment(e.c) <- parent.frame()
        aux = RandomFields::RFcov(e.c, distances=h, dim=2)
    }
    if(h.mat){
        Rho.h[low.h] = aux
        Rho.h = t(Rho.h)
        Rho.h[low.h] = aux
        return(Rho.h)
    } else {
        return(aux)
    }
}
.semiVarF = function(h, phi = NULL, sill = NULL, kappa = NULL, nugget = NULL, model = c("matern","gaussian","exponential","power","powered.exponential","cubic","circular","penta.spherical","spherical","wave","sin.hole","cauchy","gencauchy","identity","pure.nugget","RMmatern","RMwhittle","RMgneiting","RMgengneiting","RMnugget","RMcauchy","RMexp","RMgencauchy","RMgauss","RMspheric","RMstable","RMpoweredexp"), use.cor = FALSE, tol.nugget = 1.0e-15, ...){
    m.c = match.call(expand.dots = TRUE)
    cor.model = paste(formals(.autoCorF)$model)[-1]
    only.var = paste(formals(.semiVar)$model)[-1]
    if(missing(use.cor)) use.cor <- NULL
    if(is.null(use.cor)){
        in.cor <- any(grepl(model,cor.model))
        in.var <- any(grepl(model,only.var))
        use.cor <- ifelse((in.cor & in.var) | (in.cor & !in.var), TRUE,ifelse(!in.cor & in.var, FALSE, NA))
    }
    if(missing(phi) | is.null(phi)){phi = NULL}
    if(missing(sill) | is.null(sill)){sill = NULL}
    if(missing(kappa) | is.null(kappa)){kappa = NULL}
    if(missing(nugget) | is.null(nugget)){nugget = 0}
    model = match.arg(model)
    model.RF = c("RMmatern","RMwhittle","RMgneiting","RMgengneiting","RMnugget","RMcauchy","RMexp","RMgencauchy","RMgauss","RMspheric","RMstable","RMpoweredexp")
    if(!is.na(match(model,model.RF))){
        requireNamespace("RandomFields", quietly = TRUE)
    }
    if(use.cor){
        rho.v = .autoCorF(h = h, phi = phi, kappa = kappa, model = model, ...)
        return(as.vector(nugget)*(1-(h<tol.nugget)) + as.vector(sill)*(1 - rho.v))
    } else {
        gamma.h = .semiVar(model = model, what = "all", tol.nugget = tol.nugget)
        if(any(grepl("kappa",names(gamma.h$lower))) & is.null(kappa)){
            stop("Parameter kappa was not defined!")
        }
        theta.all = c(sill=as.vector(sill),phi=as.vector(phi),kappa=as.vector(kappa),nugget=as.vector(nugget))
        isnt = is.na(match(names(theta.all),names(gamma.h$lower)))
        theta = theta.all[!isnt]
        not.given = is.na(match(names(gamma.h$lower),names(theta)))
        if(any(not.given)){
            stop(paste("Parameter",paste(names(gamma.h$lower)[not.given],sep="",collapse=", "),"not given but needed for the covariance family!"))
        }
        if(any(isnt)){
            message(paste("Parameter",paste(names(theta.all)[isnt],sep="",collapse=", "),"not needed for the family, so not used!"))
        }
        for(i in 1:length(gamma.h$lower)){
            if(!(theta[i]>=gamma.h$lower[i] && theta[i]<=gamma.h$upper[i])){
                stop(paste("Parameter ",names(theta)[i]," lies outside allowed region for its family!\nEnter ",names(theta)[i]," value in ",ifelse(gamma.h$lower[i]==0,"(","["),gamma.h$lower[i],",",gamma.h$upper[i],ifelse(gamma.h$upper[i]==Inf,").","]."),sep="",collapse=""))
            }
        }
        h.mat = (is.matrix(h) || is.data.frame(h)) && nrow(h)==ncol(h)
        if(h.mat){
            if(!isSymmetric(h) || nrow(h)!=ncol(h)) stop("h must be a vector or a symmetric n*n matrix of distances!")
            G.h = matrix(nugget,nrow=nrow(h),ncol=ncol(h))
            low.h = lower.tri(G.h,diag=FALSE)
            h = h[low.h]
        }
        if(any(grepl(model,c("matern","gaussian","exponential","power","cubic","penta.spherical","spherical","wave","sin.hole","identity","pure.nugget")))){
            aux = eval(gamma.h$gamma.h)
        } else {
            call.cf = eval(gamma.h$gamma.h)
            fargs <- formals(call.cf)
            nargs <- names(formals(call.cf))
            largs <- list()
            for(i in 1:length(nargs)){
                if(any(grepl(nargs[i],names(theta)))){
                    largs[[i]] <- theta[grepl(nargs[i],names(theta))]
                } else {
                    largs[[i]] <- fargs[[i]]
                }
                names(largs)[i] <- nargs[i]
            }
            d.c <- call("call.cf")
            for(i in 1:length(largs)){
                if(names(largs)[i]=="Aniso"){
                    d.c[[names(largs)[i]]] <- m.c[["Aniso"]]
                } else if(names(largs)[i]=="proj") {
                    d.c[[names(largs)[i]]] <- m.c[["proj"]]
                } else {
                    d.c[[names(largs)[i]]] <- largs[[i]]
                }
            }
            e.c <- eval(d.c)
            environment(e.c) <- parent.frame()
            aux = RandomFields::RFvariogram(e.c, distances=h, dim=2)
        }
        if(h.mat){
            G.h[low.h] = aux
            G.h = t(G.h)
            G.h[low.h] = aux
            return(G.h)
        } else {
            return(aux)
        }
    }
}
.h2D <- function(x, angle = 0, ratio = 1){
	if(angle != 0 | ratio != 1){
		H.m = matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),2,2,byrow = TRUE)
		D.m = diag(c(1,1/ratio))
		A.m = D.m%*%H.m
	} else {
		A.m = diag(1,2)
	}
	.hNormRow <- function(coord,r.loc, A = NULL){
		d.ro = NULL
		for(r.other in 1:nrow(coord)){
			h.diff = coord[r.loc,] - coord[r.other,]
			if(!is.null(A)){
				h.diff = A%*%h.diff
			}
			d.ro[r.other] = sqrt(t(h.diff)%*%h.diff)
		}
		return(as.numeric(d.ro))
	}
	hs.row = NULL
	for(r.row in 1:nrow(x)){
		hs.row = cbind(hs.row,.hNormRow(x,r.row,A=A.m))
	}
	list(coord = x, h = hs.row)
}
.h2D <- function(x,angle=0,ratio=1) {
 if (!is.numeric(x))
 stop("argument x must be numeric")
 if (!(nrow(x)>1) | !(ncol(x)==2))
 stop("argument x must be of dimension n*2")
 if (!(ratio > 0 & ratio <= 1) | !(angle >= 0 & angle <= 360))
 stop("ratio must belong to (0,1] and angle to [0,360]")
 n = as.integer(nrow(x))
 s = apply(x,2,as.double)
 a = as.double(angle)
 r = as.double(ratio)
 h = matrix(as.double(rep(0,n*n)),n,n)
 out <- .Fortran("h2d", n=n, s=s, a=a, r=r, h=h, PACKAGE = "scpm")
 list(coord=x,h=out$h)
}
.createCov = function(dist.m = NULL, coords.m = NULL, to.Inf = 1.0e+24, rhos = NULL, ...){
    force.cor <- FALSE
    if(missing(rhos)) rhos <- NULL
	if(missing(to.Inf)) to.Inf = 1.0e+24
	if(is.null(to.Inf)) to.Inf = 1.0e+24
	add.list = list(...)
	if(!is.null(add.list$phi)){ phi = add.list$phi } else { phi = NULL }
	if(!is.null(add.list$sill)){ sill = add.list$sill } else {sill = NULL }
	if(!is.null(add.list$kappa)){ kappa = add.list$kappa } else { kappa = NULL }
	if(!is.null(add.list$nugget)){ nugget = add.list$nugget } else { nugget = NULL }
	if(!is.null(add.list$model)){ model = add.list$model } else { model = NULL }
	if(is.null(add.list$use.cor)){ use.cor = FALSE } else { use.cor = add.list$use.cor }
	if(is.null(add.list$tol.nugget)){ tol.nugget = 1.0e-15 } else { tol.nugget = add.list$tol.nugget }
	if(is.null(nugget) & is.null(sill) & is.null(rhos)){
        stop(".createCov: nugget,sill or rhos must be defined!")
    } else if(!is.null(nugget) & !is.null(sill) & is.null(rhos)){
        rhos <- sill/(nugget + sill)
    } else if(!is.null(nugget) & is.null(sill) & !is.null(rhos)){
        sill <- nugget/(1-rhos) - nugget
    } else if(is.null(nugget) & !is.null(sill) & !is.null(rhos)){
        nugget <- sill/rhos - sill
    } else if(is.null(nugget) & is.null(sill) & !is.null(rhos)){
        force.cor <- TRUE
    } else {
        stop(".createCov: nugget and sill, only rhos or rhos with nugget or sill must be specified!")
    }
    if(missing(dist.m)) dist.m <- NULL
    if(missing(coords.m)) coords.m <- NULL
	if(is.null(add.list$angle)) { angle = 0 } else { angle = add.list$angle } 
	if(is.null(add.list$ratio)) { ratio = 1 } else { ratio = add.list$ratio }
	lcor = paste(formals(.autoCorF)$model)[-1]
	lvar = paste(formals(.semiVarF)$model)[-1]
    if(is.null(dist.m) & is.null(coords.m)){
        stop("matrix of distances or coordinates must be given!")
    } else if(is.null(dist.m) & !is.null(coords.m)){
        dist.m <- .h2D(as.matrix(coords.m), angle = angle, ratio = ratio)$h
    }
	in.cor = any(!is.na(match(add.list$model,lcor)))
	in.var = any(!is.na(match(add.list$model,lvar)))
	use.var = ifelse(((in.cor & in.var) | (!in.cor & in.var)) & !force.cor,TRUE,ifelse((in.cor & in.var) & force.cor,FALSE,ifelse((!in.cor & in.var) & force.cor,TRUE,ifelse(in.cor & !in.var,FALSE,NA))))
	if(!is.na(use.var)){
		if(use.var){
			if(ifelse(in.var & in.cor & !use.cor,TRUE,FALSE)){
				use.cor = !use.cor
			}
			if(ifelse(in.var & !in.cor & use.cor,TRUE,FALSE)){
				use.cor = !use.cor
			}
            if(force.cor) stop(".createCov: model does not allow only rhos. Additional nugget or sill must be specified!")
			c.gI = call(".semiVarF",h=to.Inf,phi=phi,sill=sill,kappa=kappa,nugget=0,model=model,use.cor=use.cor)
			gamma.Inf = eval(c.gI)
			c.gh = call(".semiVarF",h=dist.m,phi=phi,sill=sill,kappa=kappa,nugget=0,model=model,use.cor=use.cor)
			gamma.h = eval(c.gh)
			c.mat = gamma.Inf - gamma.h
			gamma.n = call(".semiVarF",h=dist.m,phi=phi,sill=0,kappa=kappa,nugget=nugget,model="pure.nugget",tol.nugget=tol.nugget,use.cor=FALSE)
			c.mat = eval(gamma.n) + c.mat
			c.mat[c.mat==Inf] = to.Inf
			c.mat[c.mat==-Inf] = -to.Inf
			return(c.mat)
		} else {
            if(!force.cor){
                sills <- nugget + sill
    			c.gh = call(".semiVarF",h=dist.m,nugget=nugget,model="pure.nugget",tol.nugget=tol.nugget,use.cor=FALSE)
    			c.rh = call(".autoCorF",h=dist.m,phi=phi,kappa=kappa,model=model)
    			c.mat = (1/sills)%x%eval(c.gh) + rhos%x%eval(c.rh)
            } else {
    			c.rh = call(".autoCorF",h=dist.m,phi=phi,kappa=kappa,model=model)
                c.mat <- (1-rhos)%x%(dist.m<tol.nugget) + rhos%x%eval(c.rh)
            }
			c.mat[c.mat==Inf] = to.Inf
			c.mat[c.mat==-Inf] = -to.Inf
			return(c.mat)
		}
	} else {
		stop("model not in the list of available models!")
	}
}
.scpLoglik <- function(par, beta.par = FALSE, is.log.rhos = TRUE,
                              obs.z, dist.m, W, Xs, Z, V, fix.nugget = NULL, fix.kappa = NULL,
                              pty = "none", Q.base = NULL, Ztransf = "RV", model, ml = TRUE, use.profile = TRUE,
                              inverse.method = "solve", itol = .Machine$double.eps, createA = "eigen",
                              diagA = TRUE, nugget.tol = 1.0e-15, print.pars = TRUE, ...){
    ac <- paste(formals(.autoCor)$model)[-1]
    sv <- paste(formals(.semiVar)$model)[-1]
    ul <- unique(c(ac,sv))
    fl <- sv[is.na(match(sv,ac))]
    non.stat <- any(grepl(paste("^",model,"$",sep="",collapse=""),fl))
    nugget.par <- ifelse(non.stat & is.null(fix.nugget),TRUE,FALSE)
    fit.nugget = ifelse(is.null(fix.nugget),TRUE,FALSE)
    fit.rhos = ifelse(fit.nugget,TRUE,FALSE)
    fit.kappa = ifelse(is.null(fix.kappa),TRUE,FALSE)
    drop.nugget = ifelse(!non.stat | (non.stat & !fit.nugget),TRUE,FALSE)
    drop.rhos = ifelse((!non.stat & !fit.nugget) | (non.stat & !fit.nugget),TRUE,FALSE)
    n.xy = length(obs.z)
    if(!beta.par) add.beta = 0 else add.beta = ncol(Xs)
    if(beta.par) beta.s = par[1:add.beta] else beta.s <- NULL
    if(pty=="none") smooth <- FALSE else smooth <- TRUE
    if(smooth){
        if(pty=="cs" | pty=="ps"){
            add.alpha <- 2
        } else if(pty=="tps"){
            add.alpha <- 1
        } else if(pty=="none"){
            add.alpha <- 0
        } else { stop("Non valid smooth type definition!") }
    } else {
        add.alpha <- 0
    }
    if(smooth){
        alpha.xy = par[(add.beta+1):(add.beta+add.alpha)]
        alpha <- sum(alpha.xy)
        rho.xy = alpha.xy/alpha
    } else {
        alpha.xy = NULL
        alpha <- NULL
        rho.xy = NULL
    }
    if(!drop.rhos){
        next.pos = add.beta + add.alpha + 2
        if(is.log.rhos){
            log.rhos <- c(log.rhos = par[add.beta + add.alpha + 1]) #log(sigma2/sigma2s)
            rhos <- c(rhos = exp(log.rhos)) #sigma2/sigma2s
        } else {
            rhos <- c(rhos = par[add.beta + add.alpha + 1]) #sigma2/sigma2s
        }
    } else {
        if(!is.null(fix.nugget)){
            next.pos = add.beta + add.alpha + 1
            rhos <- c(rhos = 1 - fix.nugget)
            log.rhos <- c(log.rhos = log(rhos))
        } else {
            stop("A value for fix.nugget must be specified!")
        }
    }
    isI = ifelse(model=="identity",TRUE,FALSE)
    if(!isI){
        phi = par[next.pos]
        if(any(model==c("matern", "powered.exponential", "cauchy", "gencauchy", "RMmatern", "RMwhittle", "RMgengneiting", "RMcauchy", "RMgencauchy", "RMstable", "RMpoweredexp"))){
            if(fit.kappa){ kappa = par[(next.pos+1):ifelse(fit.nugget & !drop.nugget,length(par)-1,length(par))] } else { kappa = fix.kappa }
        } else {
            kappa = NULL
        }
    } else {
        phi = NULL
        kappa = NULL
    }
    diff.length <- length(par)-length(c(beta.s,alpha.xy,if(fit.rhos & !drop.rhos){rhos},phi,if(is.null(fix.kappa)){kappa}))
    if(sum((fit.nugget & !drop.nugget))!=sum(diff.length)){
        if(sum((fit.nugget & !drop.nugget))>sum(diff.length)) stop(paste("Length of parameter defined different than required. Perhaps parameter tau2 must be included in par for model ",model,".",sep=""))
        if(sum((fit.nugget & !drop.nugget))<sum(diff.length)) stop(paste("Length of parameter defined different than required. Perhaps parameter tau2 should not be included in par for model ",model,".",sep=""))
    }
    if(nugget.par){
        tau2 <- c(tau2 = par[length(par)])
    } else if(!is.null(fix.nugget)){
        tau2 <- c(tau2 = fix.nugget)
    } else {
        tau2 <- 0
    }
    if(!non.stat){
        R.s = .createCov(dist.m = dist.m, phi = phi, sill = NULL, kappa = kappa, nugget = NULL, rhos = rhos,
                         model = model, tol.nugget = nugget.tol, ...)
        R.u = W%*%R.s%*%t(W)
    } else {
        sigma2s <- c(sigma2s = tau2/rhos)
        sigma2 <- c(sigma2 = sigma2s - tau2)
        Sigma.s = .createCov(dist.m = dist.m, phi = phi, sill = NULL, kappa = kappa, nugget = tau2, rhos = rhos,
                             model = model, tol.nugget = nugget.tol, ...)
        R.s = (1/sigma2s)%x%Sigma.s
        Sigma.u = W%*%Sigma.s%*%t(W)
        R.u = (1/sigma2s)%x%Sigma.u
    }
    if(smooth){
        if(pty=="tps"){
            l.Q <- Q.base
            Z <- Q.base$K
            G.m = (1/alpha)%x%l.Q$A
            Gi.m = alpha%x%l.Q$Ai
        } else if(pty=="cs" | pty=="ps"){
            Q.rho = .tensorPenalty(obj = Q.base, alpha = rho.xy, type = "P")
			l.Q <- list(V = V)
			if(Ztransf=="V" | Ztransf=="VA"){
				Z <- Z
			} else {
				Z <- .generateZs(R.s, V, itol = itol)
			}
			l.Q$Ai <- .generateAi(Z, Q.rho, itol = itol)
			l.Q$A <- solve(l.Q$Ai , tol = itol)
            if(Ztransf=="VA" | Ztransf=="RVA"){
				el.A <- eigen(l.Q$A)
				l.Q$Asr <- el.A$vectors%*%diag(sqrt(el.A$values))%*%t(el.A$vectors)
				Z <- Z%*%l.Q$Asr
                G.m = (1/alpha)%x%diag(1,nrow(l.Q$A))
                Gi.m = alpha%x%diag(1,nrow(l.Q$A))
            } else {
                G.m = (1/alpha)%x%l.Q$A
                Gi.m = alpha%x%l.Q$Ai
            }
        }
        Zs <- W%*%Z
        LL = Zs%*%G.m%*%t(Zs) + R.u
    } else {
        LL = R.u
    }
    LL.i = .whichInverse(inverse.method, LL, tol = itol, Re.values = TRUE)
    if(!beta.par | !ml) XLX = t(Xs)%*%LL.i%*%Xs
    if(!beta.par){
        XLXi = .whichInverse(inverse.method,XLX, tol = itol, Re.values = TRUE)
        XLz = t(Xs)%*%LL.i%*%obs.z
        beta.s = XLXi%*%XLz
    }
    e.c = (obs.z - Xs%*%beta.s)
    RSS = t(e.c)%*%LL.i%*%e.c
    den = ifelse(ml,n.xy,n.xy-ncol(Xs))
    sigma2s = RSS/den
    tau2 <- (1-rhos)*sigma2s
    sigma2 = sigma2s - tau2
    if(ml){
        if(!use.profile){
            logL = -0.5*den*log(2*pi) -0.5*determinant(sigma2s%x%LL)$modulus -0.5*RSS
        }
        plogL = -0.5*den*log(RSS) -0.5*determinant(LL)$modulus
    } else {
        if(!use.profile){
            logL = -0.5*den*log(2*pi) -0.5*determinant(sigma2s%x%LL)$modulus -0.5*RSS -0.5*determinant(sigma2s%x%XLX)$modulus
        }
        plogL = -0.5*den*log(RSS) -0.5*determinant(LL)$modulus -0.5*determinant(XLX)$modulus
    }
    if(print.pars){
        message("##########################################")
        message("####### estimated parameter values #######")
        message("##########################################")
        message("##########################################")
        if(use.profile){
            message(paste("profile-logLik =",round(plogL,6)))
        } else {
            message(paste("logLik =",round(logL,6)))
        }
        message(paste(paste("beta.s[",1:length(beta.s),"] = ",sep=""),round(beta.s,6),sep="",collapse=" "))
        if(exists("alpha.xy")) { if(!is.null(alpha.xy)) message(paste(paste("alpha[",1:length(alpha.xy),"] = ",sep=""),round(alpha.xy,6),sep="",collapse=" ")) }
        if(exists("alpha")) { if(!is.null(alpha)) message(paste("alpha =",round(alpha,6))) }
        if(exists("rhos")) { if(!is.null(rhos)) message(paste("rhos =",round(rhos,6))) }
        if(exists("sigma2s")) { if(!is.null(sigma2s)) message(paste("sigma2s =",round(sigma2s,6))) }
        if(exists("sigma2")) { if(!is.null(sigma2)) message(paste("sigma2 =",round(sigma2,6))) }
        if(exists("phi")) { if(!is.null(phi)) message(paste("phi =",round(phi,6))) }
        if(exists("kappa")){
        if(class(kappa)!="function" & !is.null(kappa)){
            message(paste(paste("kappa[",1:length(kappa),"] =",sep=""),round(kappa,6),sep="",collapse=" "))
        }
    }
    if(exists("tau2")) { if(!is.null(tau2)) message(paste("tau2 =",round(tau2,6))) }
        message("##########################################")
    }
    if(use.profile){
        return(as.numeric(-plogL))
    } else {
        return(as.numeric(-logL))
    }
}
.randomParsGrid = function(pars.bounds, control = list()){
 if(is.null(control$is.zero)){ control$is.zero = 1.0e-6 }
 if(is.null(control$is.inf)){ control$is.inf = 10000 }
 if(is.null(control$size.each.grid)){ control$size.each.grid = 30 }
 if(is.null(control$size.sample)){ control$size.sample = 25 }
 c.zero = control$is.zero; c.Inf = control$is.inf; c.one = 1-control$is.zero
 c.ninits = control$size.each.grid; c.lsize = control$size.sample
 c.grid.0toInf = c(seq(c.zero,1,l=0.4*c.ninits), seq(1,10,l=0.4*c.ninits),
                   seq(10,100,l=0.15*c.ninits), seq(100,1000,l=0.025*c.ninits), seq(1000,c.Inf,l=0.025*c.ninits))
 c.grid.0toInf = c.grid.0toInf[c.grid.0toInf>=c.zero & c.grid.0toInf<=c.Inf]
 c.grid.mInfto0 = rev(-c.grid.0toInf)
 c.grid.mInftoInf = c(seq(c.zero,1,l=0.4*c.ninits*0.5),seq(1,10,l=0.4*c.ninits*0.5),
                      seq(10,100,l=0.15*c.ninits*0.5),seq(100,1000,l=0.025*c.ninits*0.5),seq(1000,c.Inf,l=0.025*c.ninits*0.5))
 c.grid.mInftoInf = c.grid.mInftoInf[c.grid.mInftoInf > c.zero & c.grid.mInftoInf < c.Inf]
 c.grid.mInftoInf = c(-rev(c.grid.mInftoInf),c.grid.mInftoInf)
 c.grid.0to2 = c(seq(c.zero,1,l=0.5*c.ninits),seq(1,2,l=0.5*c.ninits))
 c.grid.0to1 = c(seq(c.zero,1,l=c.ninits))
 c.pars.grid = list()
 for(i in 1:length(pars.bounds$lower)){
  if(pars.bounds$lower[i]==-Inf & pars.bounds$upper[i]==Inf) c.pars.grid[[i]] = unique(c.grid.mInftoInf)
  if(pars.bounds$lower[i]==0 & pars.bounds$upper[i]==Inf) c.pars.grid[[i]] = unique(c.grid.0toInf)
  if(pars.bounds$lower[i]==-Inf & pars.bounds$upper[i]==0) c.pars.grid[[i]] = unique(c.grid.mInfto0)
  if(pars.bounds$lower[i]==0 & pars.bounds$upper[i]==2) c.pars.grid[[i]] = unique(c.grid.0to2)
  if(pars.bounds$lower[i]==0 & pars.bounds$upper[i]==1) c.pars.grid[[i]] = unique(c.grid.0to1)
 }
 c.pars.grid = expand.grid(c.pars.grid)
 colnames(c.pars.grid) = names(pars.bounds$lower)
 c.ids = sample(1:nrow(c.pars.grid),min(nrow(c.pars.grid),c.lsize),replace=FALSE)
 c.min = as.matrix(c.pars.grid[c.ids,])
 colnames(c.min) = colnames(c.pars.grid)
 rownames(c.min) = rownames(c.pars.grid)[c.ids]
 return(c.min)
}
.getBounds = function(pars = list(fix.nugget = NULL, fix.kappa = NULL, is.log.rhos = TRUE, model = "matern"), smooth.type = c("none","cs","ps","tps")){
    .defaultArg("smooth.type","none")
    if(is.null(pars)) stop("pars list must be specified!")
    if(is.null(pars$is.log.rhos)) pars$is.log.rhos = TRUE
    ac <- paste(formals(.autoCor)$model)[-1]
    sv <- paste(formals(.semiVar)$model)[-1]
    ul <- unique(c(ac,sv))
    fl <- sv[is.na(match(sv,ac))]
    non.stat <- any(grepl(paste("^",pars$model,"$",sep="",collapse=""),fl))
    nugget.par <- ifelse(non.stat & is.null(pars$fix.nugget),TRUE,FALSE)
    fit.nugget = ifelse(is.null(pars$fix.nugget),TRUE,FALSE)
    fit.rhos = ifelse(fit.nugget,TRUE,FALSE)
    fit.kappa = ifelse(is.null(pars$fix.kappa),TRUE,FALSE)
    c.bo.sv = .tryCatchWE(.semiVar(pars$model,"bounds"))
    c.bo.ac = .tryCatchWE(.autoCor(pars$model,"bounds"))
    c.use.var = ifelse((!c.bo.sv$error & !c.bo.ac$error) | (!c.bo.sv$error & c.bo.ac$error),
                        TRUE,
                        ifelse((c.bo.sv$error & !c.bo.ac$error),FALSE,NA))
    if(!is.na(c.use.var)){
        if(c.use.var){
            c.lower = c.bo.sv$value$lower
            c.upper = c.bo.sv$value$upper
        } else {
            c.lower = c.bo.ac$value$lower
            c.upper = c.bo.ac$value$upper
        }
    } else {
        stop("Model not in the list of models!")
    }
    c.names = names(c.lower)
    if(any(grepl("sill",c.names))){
        c.id <- grepl("sill",c.names)
        c.names <- c.names[!c.id]
        c.lower <- c.lower[!c.id]
        c.upper <- c.upper[!c.id]
    }
    c.names = c(ifelse(pars$is.log.rhos,"log.rhos","rhos"),c.names)
    c.lower = c(ifelse(pars$is.log.rhos,-Inf,0),c.lower)
    c.upper = c(ifelse(pars$is.log.rhos,0,1),c.upper)
    if(smooth.type=="cs" | smooth.type=="ps"){
        c.names = c("alpha1","alpha2",c.names)
        c.lower = c(0,0,c.lower)
        c.upper = c(Inf,Inf,c.upper)
    } else if(smooth.type=="tps"){
        c.names = c("alpha",c.names)
        c.lower = c(0,c.lower)
        c.upper = c(Inf,c.upper)
    }
    drop.nugget = ifelse(!non.stat | (non.stat & !fit.nugget),TRUE,FALSE)
    if(any(grepl("nugget",c.names))){
        if(!drop.nugget){
            c.names[grepl("nugget",c.names)] = "tau2"
        } else {
            c.id = grepl("nugget",c.names)
            c.names = c.names[!c.id]
            c.lower = c.lower[!c.id]
            c.upper = c.upper[!c.id]
        }
    } else {
        if(!drop.nugget){
            c.names = c(c.names,"tau2")
            c.lower = c(c.lower,0)
            c.upper = c(c.upper,Inf)
        }
    }
    drop.rhos = ifelse((!non.stat & !fit.nugget) | (non.stat & !fit.nugget),TRUE,FALSE)
    if(any(grepl("(log.rhos|rhos)",c.names))){
        if(drop.rhos){
            c.id = grepl("(log.rhos|rhos)",c.names)
            c.names = c.names[!c.id]
            c.lower = c.lower[!c.id]
            c.upper = c.upper[!c.id]
        }
    }
    drop.kappa = ifelse(is.null(pars$fix.kappa),FALSE,ifelse(is.logical(pars$fix.kappa),ifelse(pars$fix.kappa,TRUE,FALSE),TRUE))
    if(any(grepl("kappa",c.names))){
        if(drop.kappa){
            c.id = grepl("kappa",c.names)
            c.names = c.names[!c.id]
            c.lower = c.lower[!c.id]
            c.upper = c.upper[!c.id]
        }
    }
    names(c.lower) = c.names
    names(c.upper) = c.names
    list(lower = c.lower, upper = c.upper)
}
.gs = function(zV, XM, ZM, VM, hM, Qel, W,
                          use = list(penalty = "cs", ps.order = 2, Ztr = "RVA",
                                     inverse = ".invEigen", tol = .Machine$double.neg.eps,
                                     Amethod = "eigen", Adiag = TRUE),
                          pars = list(addbeta = FALSE, is.log.rhos = TRUE,
                                      fix.nugget = NULL, fix.kappa = NULL,
                                      model = "matern", nugget.tol = 1.0e-15),
                          method = list(use.reml = FALSE, use.profile = TRUE),
                          control = list(is.zero = 1.0e-6, is.inf = 10000, size.each.grid = 30,
                                         size.sample = 40),
                          print.pars = FALSE){
    if(is.null(use)) stop("List use must be specified!")
    if(is.null(pars)) stop("List pars must be specified!")
    if(is.null(method)) stop("List method must be specified!")
    if(is.null(control)) stop("List control must be specified!")
    if(is.null(use$penalty)) stop("Element penalty must be defined in use list!")
    if(is.null(use$ps.order) & use$penalt=="ps") stop("Element ps.order must be defined in use list!")
    if(is.null(use$Ztr)) stop("Element Ztr must be defined in use list!")
    if(is.null(use$inverse)) stop("Element inverse must be defined in use list!")
    if(is.null(use$tol)) stop("Element tol must be defined in use list!")
    if(is.null(use$Amethod)) stop("Element Amethod must be defined in use list!")
    if(is.null(use$Adiag)) stop("Element Adiag must be defined in use list!")
    if(is.null(pars$addbeta)) stop("Element addbeta must be defined in pars list!")
    if(is.null(pars$is.log.rhos)) stop("Element is.log.rhos must be defined in pars list!")
    if(is.null(pars$model)) stop("Element model must be defined in pars list!")
    if(is.null(pars$nugget.tol)) stop("Element nugget.tol must be defined in pars list!")
    if(is.null(method$use.reml)) stop("Element use.reml must be defined in method list!")
    if(is.null(method$use.profile)) stop("Element use.profile must be defined in method list!")
    if(is.null(control$is.zero)) stop("Element is.zero must be defined in control list!")
    if(is.null(control$is.inf)) stop("Element is.inf must be defined in control list!")
    if(is.null(control$size.each.grid)) stop("Element size.each.grid must be defined in control list!")
    if(is.null(control$size.sample)) stop("Element size.sample must be defined in control list!")
    use$Adiag = ifelse(use$penalty=="ps", TRUE, use$Adiag)
    c.theta = .getBounds(pars = pars, smooth.type=use$penalty)
    c.npars = length(c.theta$lower)
    c.nat = control$size.each.grid**c.npars
    if(c.nat > 35**5){
        c.base.size = 10
        while(c.base.size**c.npars < 35**5){
            c.base.size = c.base.size + 1
        }
        control$size.each.grid = c.base.size - 1
        message(paste("The value size.each.grid**length(pars) cannot be greated than 50**5.\nChanging value to size.each.grid = ",c.base.size-1,".\nTo change manually modify control and pars lists!",sep=""))
    }
    c.start = c.theta
    if(pars$addbeta){
        c.start$lower = c(beta=rep(-Inf,ncol(XM)),c.theta$lower)
        c.start$upper = c(beta=rep(+Inf,ncol(XM)),c.theta$upper)
    }
    c.init.pars = .randomParsGrid(c.theta, control = control)
    c.ilL = NULL
    for(i in 1:nrow(c.init.pars)){
		message(paste("Step",i,"of",nrow(c.init.pars)))
        aux = .tryCatchWE(.scpLoglik(
            par = c.init.pars[i,],
            beta.par = FALSE,
            is.log.rhos = pars$is.log.rhos,
            obs.z = zV,
            dist.m = hM,
            W = W,
            Xs = XM,
			Z = ZM,
			V = VM,
            fix.nugget = pars$fix.nugget,
            fix.kappa = pars$fix.kappa,
            pty = use$penalty,
            Q.base = Qel,
            Ztransf = use$Ztr,
            model = pars$model,
            ml = !method$use.reml,
            use.profile = method$use.profile,
            inverse.method = use$inverse,
            itol = use$tol,
            createA = use$Amethod,
            diagA = use$Adiag,
            nugget.tol = pars$nugget.tol,
            print.pars = print.pars
        ))
        c.ilL[i] = ifelse(!aux$error, aux$value, NA)
    }
    c.start$pars = c.init.pars[match(min(c.ilL, na.rm = TRUE),c.ilL),]
    if(length(c.start$pars)==1)
        names(c.start$pars) <- colnames(c.init.pars)
    if(pars$addbeta){
        c.beta = .betaFromPars(pars = c.start$pars, zV = zV, XM = XM, ZM = ZM, VM = VM, hM = hM, Qel = Qel, W,
                              use = use, pars.str = pars, use.reml = method$use.reml, only.beta = TRUE)
        c.start$pars = c(c.beta,c.start$pars)
    }
    return(c.start)
}
.rightDistance = function(z,psi){
    aux = cbind((z-psi)*(z>psi),-(z>psi))
    colnames(aux) = paste(c("U(","V("),psi,")",sep="")
    return(aux)
}
.matrixUV = function(z=NULL,psi){
    zMatrix = matrix(0,nrow=length(z),ncol=2*length(psi))
    cnamesz = rep("",ncol=2*length(psi))
    for(i in 1:length(psi)){
        aux = .rightDistance(z,psi[i])
        cnamesz[(2*i-1):(2*i)] = colnames(aux)
        zMatrix[,(2*i-1):(2*i)] = aux
    }
    colnames(zMatrix) = cnamesz
    return(zMatrix)
}
.psiUpdate = function(beta, psi0, UVcols){
	if(length(psi0)*2!=length(UVcols)) stop("Length of psi must match the number of columns of UV matrix!")
    return(psi0 + (beta[seq(min(UVcols)+1,max(UVcols),2)])/(beta[seq(min(UVcols),max(UVcols)-1,2)]))
}
.psiIterative = function(cp.call, psi0, UVcols, obj, control.cp = list(maxiter = 20, maxskip = 10)){
	psi0.v = NULL
	UVcols.v = NULL
	for(i in 1:length(psi0)){
		psi0.v = c(psi0.v,psi0[[i]])
		UVcols.v = c(UVcols.v,UVcols[[i]])
	}
	beta.step = matrix(0,nrow=control.cp$maxiter,ncol=ncol(obj$XM))
	psi.step = matrix(0,nrow=control.cp$maxiter,ncol=length(psi0.v))
	beta.step[1,] = rep(0,ncol(obj$XM))
	psi.step[1,] = psi0.v
	num.skip = 0
	ii = 2
	while(ii < (control.cp$maxiter + 1) && num.skip < control.cp$maxskip){
		if(ii>2){
			acc.psi = 0
			newUV = NULL
			for(i in 1:length(cp.call)){
				mp = match("psi",names(cp.call[[i]]))
				cp.call[[i]][[mp]] = psi.step[ii-1,(acc.psi+1):(acc.psi+length(psi0[[i]]))]
				newUV[[i]] = eval(cp.call[[i]], parent.frame())
				obj$XM[,UVcols[[i]]] = newUV[[i]]$UV
				acc.psi = acc.psi + length(psi0[[i]])
			}
		}
		cat(".")
		beta.aux = .betaFromRSPZ(obj)
		beta.step[ii,] = beta.aux$beta
		psi.step[ii,] = .psiUpdate(beta.step[ii,], psi.step[ii-1,], UVcols.v)
		ii = ii + 1
	}
	psi = psi.step[nrow(psi.step),]
	acc.psi = 0
	psi.l = NULL
	for(i in 1:length(cp.call)){
		psi.l[[i]] = psi[(acc.psi+1):(acc.psi+length(psi0[[i]]))]
		acc.psi = acc.psi + length(psi0[[i]])
	}
	list(psi = psi.l, beta = beta.step[nrow(beta.step),], XM = obj$XM, UVcols = UVcols,
		 elements = beta.aux, history = list(beta = beta.step, psi = psi.step))
}
.RSPZfromPars = function(pars, zV, XM, ZM, VM, hM, Qel, W,
                              use = list(penalty = "cs", ps.order = 2, Ztr = "RVA",
                                         inverse = ".invEigen", tol = .Machine$double.neg.eps,
                                         Amethod = "eigen", Adiag = TRUE),
                              pars.str = list(fix.nugget = NULL, fix.kappa = NULL, is.log.rhos = TRUE,
                                              model = "matern", nugget.tol = 1.0e-15)){
    if(missing(Qel)) Qel <- NULL
    if(missing(use)) use <- NULL
    if(is.null(use)){
        use = list(penalty = "none", inverse = ".invEigen", tol = .Machine$double.neg.eps)
    }
    if(is.null(use$penalty)) use$penalty <- "none"
    if(is.null(use$intercept)) use$intercept <- TRUE
    if(is.null(pars.str)){
        stop("pars.str must be defined!")
    } else {
        bounds = .getBounds(pars = list(fix.nugget = pars.str$fix.nugget,
                                fix.kappa = pars.str$fix.kappa,
                                is.log.rhos = pars.str$is.log.rhos,
                                model = pars.str$model), smooth.type = use$penalty)
        c.req = names(bounds$lower)
    }
    c.id.el = match(c.req,names(pars))
    if(all(!is.na(c.id.el))){
        pars = pars[c.id.el]
    } else {
        stop(paste("One or more required elements of pars were not specified!\nEnter ",paste(c.req[is.na(c.id.el)],sep="",collapse=", "),".",sep=""))
    }
    ac <- paste(formals(.autoCor)$model)[-1]
    sv <- paste(formals(.semiVar)$model)[-1]
    ul <- unique(c(ac,sv))
    fl <- sv[is.na(match(sv,ac))]
    non.stat <- any(grepl(paste("^",pars.str$model,"$",sep="",collapse=""),fl))
    nugget.par <- ifelse(non.stat & is.null(pars.str$fix.nugget),TRUE,FALSE)
    fit.nugget = ifelse(is.null(pars.str$fix.nugget),TRUE,FALSE)
    fit.rhos = ifelse(fit.nugget,TRUE,FALSE)
    fit.kappa = ifelse(is.null(pars.str$fix.kappa),TRUE,FALSE)
    drop.nugget = ifelse(!non.stat | (non.stat & !fit.nugget),TRUE,FALSE)
    drop.rhos = ifelse((!non.stat & !fit.nugget) | (non.stat & !fit.nugget),TRUE,FALSE)
    if(use$penalty=="none") smooth <- FALSE else smooth <- TRUE
    if(smooth){
        if(use$penalty=="cs" | use$penalty=="ps"){
            c.id.a = match(c("alpha1","alpha2"),names(pars))
        } else if(use$penalty=="tps"){
            c.id.a = match(c("alpha"),names(pars))
        } else if(use$penalty!="none"){ stop("Non valid smooth type definition!") }
        if(all(!is.na(c.id.a))) c.alpha.xy = pars[c.id.a] else stop("One or more elements of alpha not specified!")
        c.alpha <- sum(c.alpha.xy)
        c.rho.xy = c.alpha.xy/c.alpha
    } else {
        c.alpha.xy = NULL
        c.alpha <- NULL
        c.rho.xy = NULL
    }
    c.id.p = match("phi",names(pars))
    if(all(!is.na(c.id.p))) c.phi = pars[c.id.p] else c.phi = NULL
    if(is.null(pars.str$fix.kappa)){
        c.id.k = grepl("kappa",names(pars))
        if(any(c.id.k)) c.kappa = pars[c.id.k] else c.kappa = NULL
    } else {
        c.kappa = pars.str$fix.kappa
    }
    if(fit.nugget & !drop.nugget){
        c.id.t = match("tau2",names(pars))
        if(all(!is.na(c.id.t))) c.tau2 = pars[c.id.t] else c.tau2 = NULL
    } else if(!is.null(pars.str$fix.nugget)){
        c.tau2 <- NULL
        c.rhos <- 1 - pars.str$fix.nugget
        c.log.rhos <- log(c.rhos)
    } else {
        c.tau2 <- NULL
    }
    if(fit.rhos & !drop.rhos){
        c.id.s = match(ifelse(pars.str$is.log.rhos,"log.rhos","rhos"),names(pars))
        if(all(!is.na(c.id.s))){
            if(pars.str$is.log.rhos){
                c.log.rhos = pars[c.id.s]
                c.rhos = exp(c.log.rhos)
            } else {
                c.rhos = pars[c.id.s]
                c.log.rhos = log(c.rhos)
            }
        }
    } else if(!is.null(pars.str$fix.nugget)){
        c.tau2 <- NULL
        c.rhos <- 1 - pars.str$fix.nugget
        c.log.rhos <- log(c.rhos)
    } else {
        c.tau2 <- 0
        c.rhos <- 1
        c.log.rhos <- 0
    }
    if(!non.stat){
        c.R.s = .createCov(dist.m = hM, phi = c.phi, sill = NULL, kappa = c.kappa, nugget = NULL, rhos = c.rhos,
                         model = pars.str$model, tol.nugget = pars.str$nugget.tol)
        c.R.u = W%*%c.R.s%*%t(W)
    } else {
        c.sigma2s <- c.tau2/c.rhos
        c.sigma2 <- c.sigma2s - c.tau2
        c.Sigma.s = .createCov(dist.m = hM, phi = c.phi, sill = NULL, kappa = c.kappa, nugget = c.tau2, rhos = c.rhos,
                             model = pars.str$model, tol.nugget = pars.str$nugget.tol)
        c.R.s = (1/c.sigma2s)%x%c.Sigma.s
        c.Sigma.u = W%*%c.Sigma.s%*%t(W)
        c.R.u = (1/c.sigma2s)%x%c.Sigma.u
    }
    c.H.u = .whichInverse(use$inverse, c.R.u, tol = use$tol, Re.values = TRUE)
    if(smooth){
        if(use$penalty=="tps"){
            c.Ps <- Qel$K
            c.l.Ps <- Qel
            c.Z <- Qel$K
            c.Gm = (1/c.alpha)%x%c.l.Ps$A
            c.Gi = c.alpha%x%c.l.Ps$Ai
        } else if(use$penalty=="cs" | use$penalty=="ps"){
            c.Ps = .tensorPenalty(obj = Qel, alpha = c.rho.xy, type = "P")
			c.l.Ps <- list(V = VM)
			if(use$Ztr=="V" | use$Ztr=="VA"){
				c.Z <- ZM
			} else {
				c.Z <- .generateZs(c.R.s, VM, itol = use$tol)
			}
			c.l.Ps$Ai <- .generateAi(c.Z, c.Ps, itol = use$tol)
			c.l.Ps$A <- solve(c.l.Ps$Ai , tol = use$tol)
            if(use$Ztr=="VA" | use$Ztr=="RVA"){
				el.A <- eigen(c.l.Ps$A)
				c.l.Ps$Asr <- el.A$vectors%*%diag(sqrt(el.A$values))%*%t(el.A$vectors)
				c.Z <- c.Z%*%c.l.Ps$Asr
                c.Gm = (1/c.alpha)%x%diag(1,nrow(c.l.Ps$A))
                c.Gi = c.alpha%x%diag(1,nrow(c.l.Ps$A))
            } else {
                c.Gm = (1/c.alpha)%x%c.l.Ps$A
                c.Gi = c.alpha%x%c.l.Ps$Ai
            }
        }
        c.Zs <- W%*%c.Z
        c.LL = c.Zs%*%c.Gm%*%t(c.Zs) + c.R.u
        c.ZHZa = t(c.Zs)%*%c.H.u%*%c.Zs + c.Gi
        c.ZHZai = .whichInverse(use$inverse, c.ZHZa, tol = use$tol, Re.values = TRUE)
        c.tZ = c.ZHZai%*%t(c.Zs)%*%c.H.u
    } else {
        c.LL = c.R.u
        c.ZHZa = NULL
        c.ZHZai = NULL
        c.tZ = NULL
    }
    c.In = diag(1,nrow(XM))
    c.LLi = .whichInverse(use$inverse, c.LL, tol = use$tol, Re.values = TRUE)
    list(alpha = if(exists("c.alpha.xy")){c.alpha.xy}, sigma2s=if(exists("c.sigma2s")){c.sigma2s}, sigma2s=if(exists("c.sigma2")){c.sigma2}, rhos = c.rhos, log.rhos = c.log.rhos, phi = c.phi, kappa = c.kappa, tau2 = c.tau2, R.s = c.R.s, R.u = c.R.u, H.u = c.H.u, Sigma.s = if(exists("c.Sigma.s")){c.Sigma.s}, Sigma.u = if(exists("c.Sigma.u")){c.Sigma.u}, W = W, XM = XM, Z = if(exists("c.Z")){c.Z}, Zs = if(exists("c.Zs")){c.Zs}, G = if(exists("c.Gm")){c.Gm}, Gi = if(exists("c.Gi")){c.Gi}, L = c.LL, Li = c.LLi, TZ = if(exists("c.tZ")){c.tZ}, ZHZa = if(exists("c.ZHZa")){c.ZHZa}, ZHZai = if(exists("c.ZHZai")){c.ZHZai}, zV = zV, use = use, pars.str = pars.str, penalty = if(use$penalty!="none"){list(P = c.Ps, decompose = c.l.Ps)}, fit.nugget = fit.nugget, Qel = Qel)
}
.betaFromPars = function(pars, zV, XM, ZM, VM, hM, Qel, W,
                              use = list(penalty = "cs", ps.order = 2, Ztr = "RVA",
                                         inverse = ".invEigen", tol = .Machine$double.neg.eps,
                                         Amethod = "eigen", Adiag = TRUE),
                              pars.str = list(fix.nugget = NULL, fix.kappa = NULL, is.log.rhos = TRUE,
                                              model = "matern", nugget.tol = 1.0e-15),
                              use.reml = FALSE, only.beta = TRUE){
    if(missing(Qel)) Qel <- NULL
    if(missing(use)) use <- NULL
    if(is.null(use)){
        use = list(penalty = "none", inverse = ".invEigen", tol = .Machine$double.neg.eps)
    }
    if(is.null(use$penalty)) use$penalty <- "none"
    if(is.null(use$intercept)) use$intercept <- TRUE
    if(is.null(pars.str)){
        stop("pars.str must be defined!")
    } else {
        bounds = .getBounds(pars = list(fix.nugget = pars.str$fix.nugget,
                                fix.kappa = pars.str$fix.kappa,
                                is.log.rhos = pars.str$is.log.rhos,
                                model = pars.str$model), smooth.type = use$penalty)
        c.req = names(bounds$lower)
    }
    c.id.el = match(c.req,names(pars))
    if(all(!is.na(c.id.el))){
        pars = pars[c.id.el]
    } else {
        stop(paste("One or more required elements of pars were not specified!\nEnter ",paste(c.req[is.na(c.id.el)],sep="",collapse=", "),".",sep=""))
    }
    ac <- paste(formals(.autoCor)$model)[-1]
    sv <- paste(formals(.semiVar)$model)[-1]
    ul <- unique(c(ac,sv))
    fl <- sv[is.na(match(sv,ac))]
    non.stat <- any(grepl(paste("^",pars.str$model,"$",sep="",collapse=""),fl))
    nugget.par <- ifelse(non.stat & is.null(pars.str$fix.nugget),TRUE,FALSE)
    fit.nugget = ifelse(is.null(pars.str$fix.nugget),TRUE,FALSE)
    fit.rhos = ifelse(fit.nugget,TRUE,FALSE)
    fit.kappa = ifelse(is.null(pars.str$fix.kappa),TRUE,FALSE)
    drop.nugget = ifelse(!non.stat | (non.stat & !fit.nugget),TRUE,FALSE)
    drop.rhos = ifelse((!non.stat & !fit.nugget) | (non.stat & !fit.nugget),TRUE,FALSE)
    if(use$penalty=="none") smooth <- FALSE else smooth <- TRUE
    if(smooth){
        if(use$penalty=="cs" | use$penalty=="ps"){
            c.id.a = match(c("alpha1","alpha2"),names(pars))
        } else if(use$penalty=="tps"){
            c.id.a = match(c("alpha"),names(pars))
        } else if(use$penalty!="none"){ stop("Non valid smooth type definition!") }
        if(all(!is.na(c.id.a))) c.alpha.xy = pars[c.id.a] else stop("One or more elements of alpha not specified!")
        c.alpha <- sum(c.alpha.xy)
        c.rho.xy = c.alpha.xy/c.alpha
    } else {
        c.alpha.xy = NULL
        c.alpha <- NULL
        c.rho.xy = NULL
    }
    c.id.p = match("phi",names(pars))
    if(all(!is.na(c.id.p))) c.phi = pars[c.id.p] else c.phi = NULL
    c.fit.nugget <- fit.nugget
    c.fit.kappa <- fit.kappa
    if(is.null(pars.str$fix.kappa)){
        c.id.k = grepl("kappa",names(pars))
        if(any(c.id.k)) c.kappa = pars[c.id.k] else c.kappa = NULL
    } else {
        c.kappa = pars.str$fix.kappa
    }
    if(fit.nugget & !drop.nugget){
        c.id.t = match("tau2",names(pars))
        if(all(!is.na(c.id.t))) c.tau2 = pars[c.id.t] else c.tau2 = NULL
    } else if(!is.null(pars.str$fix.nugget)){
        c.tau2 <- NULL
        c.rhos <- 1 - pars.str$fix.nugget
        c.log.rhos <- log(c.rhos)
    } else {
        c.tau2 <- 0
        c.rhos <- 1
        c.log.rhos <- 0
    }
    if(fit.rhos & (!drop.rhos)){
        c.id.s = match(ifelse(pars.str$is.log.rhos,"log.rhos","rhos"),names(pars))
        if(all(!is.na(c.id.s))){
            if(pars.str$is.log.rhos){
                c.log.rhos = pars[c.id.s]
                c.rhos = exp(c.log.rhos)
            } else {
                c.rhos = pars[c.id.s]
                c.log.rhos = log(c.rhos)
            }
        }
    } else if(!is.null(pars.str$fix.nugget)){
        c.tau2 <- NULL
        c.rhos <- 1 - pars.str$fix.nugget
        c.log.rhos <- log(c.rhos)
    } else {
        c.tau2 <- NULL
        c.rhos <- 1
        c.log.rhos <- 0
    }
    if(!non.stat){
        c.R.s = .createCov(dist.m = hM, phi = c.phi, sill = NULL, kappa = c.kappa, nugget = NULL, rhos = c.rhos,
                         model = pars.str$model, tol.nugget = pars.str$nugget.tol)
        c.R.u = W%*%c.R.s%*%t(W)
    } else {
        c.sigma2s <- c.tau2/c.rhos
        c.sigma2 <- c.sigma2s - c.tau2
        c.Sigma.s = .createCov(dist.m = hM, phi = c.phi, sill = NULL, kappa = c.kappa, nugget = c.tau2, rhos = c.rhos,
                             model = pars.str$model, tol.nugget = pars.str$nugget.tol)
        c.R.s = (1/c.sigma2s)%x%c.Sigma.s
        c.Sigma.u = W%*%c.Sigma.s%*%t(W)
        c.R.u = (1/c.sigma2s)%x%c.Sigma.u
    }
    c.H.u = .whichInverse(use$inverse, c.R.u, tol = use$tol, Re.values = TRUE)
    c.In = diag(1,nrow(XM))
    if(smooth){
        if(use$penalty=="tps"){
            c.Ps <- Qel$K
            c.l.Ps <- Qel
            c.Z <- Qel$K
            c.Gm = (1/c.alpha)%x%c.l.Ps$A
            c.Gi = c.alpha%x%c.l.Ps$Ai
        } else if(use$penalty=="cs" | use$penalty=="ps"){
            c.Ps = .tensorPenalty(obj = Qel, alpha = c.rho.xy, type = "P")
			c.l.Ps <- list(V = VM)
			if(use$Ztr=="V" | use$Ztr=="VA"){
				c.Z <- ZM
			} else {
				c.Z <- .generateZs(c.R.s, VM, itol = use$tol)
			}
			c.l.Ps$Ai <- .generateAi(c.Z, c.Ps, itol = use$tol)
			c.l.Ps$A <- solve(c.l.Ps$Ai , tol = use$tol)
            if(use$Ztr=="VA" | use$Ztr=="RVA"){
				el.A <- eigen(c.l.Ps$A)
				c.l.Ps$Asr <- el.A$vectors%*%diag(sqrt(el.A$values))%*%t(el.A$vectors)
				c.Z <- c.Z%*%c.l.Ps$Asr
                c.Gm = (1/c.alpha)%x%diag(1,nrow(c.l.Ps$A))
                c.Gi = c.alpha%x%diag(1,nrow(c.l.Ps$A))
            } else {
                c.Gm = (1/c.alpha)%x%c.l.Ps$A
                c.Gi = c.alpha%x%c.l.Ps$Ai
            }
        }
        c.Zs <- W%*%c.Z
        c.LL = c.Zs%*%c.Gm%*%t(c.Zs) + c.R.u
        c.ZHZa = t(c.Zs)%*%c.H.u%*%c.Zs + c.Gi
        c.ZHZai = .whichInverse(use$inverse, c.ZHZa, tol = use$tol, Re.values = TRUE)
        c.tZ = c.ZHZai%*%t(c.Zs)%*%c.H.u
        c.pZ = c.Zs%*%c.tZ
        c.IpZ = c.In - c.pZ
    } else {
        c.Zs <- NULL
        c.LL = c.R.u
        c.ZHZa = NULL
        c.ZHZai = NULL
        c.tZ = NULL
        c.pZ = NULL
        c.IpZ = c.In
    }
    c.LLi = .whichInverse(use$inverse, c.LL, tol = use$tol, Re.values = TRUE)
    c.XLX = t(XM)%*%c.LLi%*%XM
    c.XLXi = .whichInverse(use$inverse,c.XLX, tol = use$tol, Re.values = TRUE)
    c.tXs = c.XLXi%*%t(XM)%*%c.LLi
    c.beta.s = c.tXs%*%zV
    names(c.beta.s) = paste("beta",1:length(c.beta.s),sep="")
    if(only.beta){
        names(c.beta.s) <- colnames(XM)
        return(c.beta.s)
    } else {
        if(smooth){
            c.r = c.tZ%*%(zV-XM%*%c.beta.s)
            c.yf = XM%*%c.beta.s + c.Zs%*%c.r
        } else {
            c.r <- NULL
            c.yf = XM%*%c.beta.s
        }
        if(use$is.X!="none"){
            c.Vdim <- ifelse(use$is.X=="tensor",3+ifelse(use$intercept,1,0),2+ifelse(use$intercept,1,0))
            if(ncol(XM)>c.Vdim){
                c.C = XM[,-((ncol(XM)-(c.Vdim-1)):ncol(XM))]
                c.X = XM[,(ncol(XM)-(c.Vdim-1)):ncol(XM)]
                if(smooth){
                    c.pZ = c.Zs%*%c.tZ
                    c.IpZ = c.In - c.pZ
                } else {
                    c.pZ = NULL
                    c.IpZ = c.In
                }
                c.XHX = t(c.X)%*%c.H.u%*%c.IpZ%*%c.X
                c.XHXi = .whichInverse(use$inverse, c.XHX, tol = use$tol, Re.values = TRUE)
                c.tX = c.XHXi%*%t(c.X)%*%c.H.u%*%c.IpZ
                c.pX = c.X%*%c.tX
                c.IpX <- c.In - c.pX
                if(smooth){
                    c.Sxz = c.pX + c.pZ%*%c.IpX
                } else {
                    c.Sxz <- c.pX
                }
                c.ISxz = c.In - c.Sxz
                c.CHISC = t(c.C)%*%c.H.u%*%c.ISxz%*%c.C
                c.CHISCi = .whichInverse(use$inverse, c.CHISC, tol = use$tol, Re.values = TRUE)
                c.tC = c.CHISCi%*%t(c.C)%*%c.H.u%*%c.ISxz
                c.pC = c.C%*%c.tC
                c.IpC = c.In-c.pC
                c.Mxz = c.pC + c.Sxz%*%c.IpC
                if(smooth){
                    c.gf = c.X%*%(c.beta.s[(ncol(XM)-(c.Vdim-1)):ncol(XM)]) + c.Zs%*%c.r
                } else {
                    c.gf = c.X%*%(c.beta.s[(ncol(XM)-(c.Vdim-1)):ncol(XM)])
                }
            } else {
                c.X <- XM
                if(smooth){
                    c.pZ = c.Zs%*%c.tZ
                    c.IpZ = c.In - c.pZ
                } else {
                    c.pZ = NULL
                    c.IpZ = c.In
                }
                c.XHX = t(c.X)%*%c.H.u%*%c.IpZ%*%c.X
                c.XHXi = .whichInverse(use$inverse, c.XHX, tol = use$tol, Re.values = TRUE)
                c.tX = c.XHXi%*%t(c.X)%*%c.H.u%*%c.IpZ
                c.pX = c.X%*%c.tX
                c.IpX <- c.In - c.pX
                if(smooth){
                    c.Sxz = c.pX + c.pZ%*%c.IpX
                } else {
                    c.Sxz <- c.pX
                }
                c.Mxz = c.Sxz
                c.ISxz = c.In - c.Sxz
                c.gf = c.yf
            }
        } else {
        	c.Mxz = XM%*%c.tXs
            c.Sxz = NULL
            c.gf = NULL
        }
        c.epsilon = zV - XM%*%c.beta.s
        c.RSS = t(c.epsilon)%*%c.LLi%*%(c.epsilon)
        if(!use.reml) c.df = nrow(XM) else c.df = (nrow(XM)-ncol(XM))
        c.sigma2s = c.RSS/c.df
        c.tau2 <- (1-c.rhos)*c.sigma2s
        c.sigma2 = c.sigma2s - c.tau2
        c.Sigma.u = c.sigma2s%x%c.R.u
        colnames(c.XLX) <- colnames(XM)
        rownames(c.XLX) <- colnames(XM)
        colnames(c.XLXi) <- colnames(XM)
        rownames(c.XLXi) <- colnames(XM)
        names(c.beta.s) <- colnames(XM)
        list(alpha = if(exists("c.alpha.xy")){c.alpha.xy}, sigma2s = c.sigma2s, sigma2 = c.sigma2, rhos = c.rhos, log.rhos = c.log.rhos, phi = c.phi, kappa = c.kappa, tau2 = c.tau2, beta = c.beta.s, r = c.r, R = c.R.u, g = if(exists("c.gf")){c.gf}, fitted.values = if(exists("c.yf")){c.yf}, S = if(exists("c.Sxz")){c.Sxz}, M = if(exists("c.Mxz")){c.Mxz}, X = XM, Z = if(exists("c.Z")){c.Z}, Zs = if(exists("c.Zs")){c.Zs}, W = W, XLX = c.XLX, XLXi = c.XLXi, var.beta = c.sigma2s%x%c.XLXi, conf.error = if(!only.beta){sqrt(diag(c.Mxz%*%c.Sigma.u%*%t(c.Mxz)))}, pred.error = if(!only.beta){sqrt(diag(c.Sigma.u + c.Mxz%*%c.Sigma.u%*%t(c.Mxz)))}, penalty = if(use$penalty!="none"){list(P = c.Ps, decompose = c.l.Ps)})
    }
}
.betaFromRSPZ = function(obj){
    c.In = diag(1,nrow(obj$XM))
    c.XLX = t(obj$XM)%*%obj$Li%*%obj$XM
    c.XLXi = .whichInverse(obj$use$inverse, c.XLX, tol = obj$use$tol, Re.values = TRUE)
    c.tXs = c.XLXi%*%t(obj$XM)%*%obj$Li
    c.beta.s = c.tXs%*%obj$zV
    names(c.beta.s) = paste("beta",1:length(c.beta.s),sep="")
    list(beta = c.beta.s, XLX = c.XLX, XLXi = c.XLXi, TX = c.tXs)
}
.fitFromRSPZ = function(RZ, IP, use.reml = FALSE){
	c.fit.nugget = RZ$fit.nugget
	c.Ps = RZ$penalty$P
	c.l.Ps = RZ$penalty$decompose
	use = RZ$use
    if(use$penalty=="none") smooth <- FALSE else smooth <- TRUE
	pars.str = RZ$pars.str
    c.rhos = c(rhos = RZ$rhos)
    c.log.rhos = c(log.rhos = RZ$log.rhos)
	c.tau2 = RZ$tau2
	XM = IP$XM
	zV = RZ$zV
	c.LLi = RZ$Li
	c.XLX = t(XM)%*%c.LLi%*%XM
	c.XLXi = .whichInverse(use$inverse, c.XLX, tol = use$tol, Re.values = TRUE)
	c.R.u = RZ$R.u
	c.H.u = RZ$H.u
	c.Sigma.u = RZ$Sigma.u
	c.Z = RZ$Z
    c.Zs = RZ$Zs
	c.tZ = RZ$TZ
	c.tXs = c.XLXi%*%t(XM)%*%c.LLi
	c.beta.s = IP$beta
	c.In = diag(1,nrow(XM))
    if(smooth){
        c.r = c.tZ%*%(zV-XM%*%c.beta.s)
        c.yf = XM%*%c.beta.s + c.Zs%*%c.r
    } else {
        c.r <- NULL
        c.yf = XM%*%c.beta.s
    }
    if(use$is.X!="none"){
        c.Vdim <- ifelse(use$is.X=="tensor",3+ifelse(use$intercept,1,0),2+ifelse(use$intercept,1,0))
    	if(ncol(XM)>c.Vdim){
            c.C = XM[,-((ncol(XM)-(c.Vdim-1)):ncol(XM))]
            c.X = XM[,(ncol(XM)-(c.Vdim-1)):ncol(XM)]
            if(smooth){
                c.pZ = c.Zs%*%c.tZ
                c.IpZ = c.In - c.pZ
            } else {
                c.pZ = NULL
                c.IpZ = c.In
            }
            c.XHX = t(c.X)%*%c.H.u%*%c.IpZ%*%c.X
            c.XHXi = .whichInverse(use$inverse, c.XHX, tol = use$tol, Re.values = TRUE)
            c.tX = c.XHXi%*%t(c.X)%*%c.H.u%*%c.IpZ
            c.pX = c.X%*%c.tX
            c.IpX <- c.In - c.pX
            if(smooth){
                c.Sxz = c.pX + c.pZ%*%c.IpX
            } else {
                c.Sxz <- c.pX
            }
            c.ISxz = c.In - c.Sxz
            c.CHISC = t(c.C)%*%c.H.u%*%c.ISxz%*%c.C
            c.CHISCi = .whichInverse(use$inverse, c.CHISC, tol = use$tol, Re.values = TRUE)
            c.tC = c.CHISCi%*%t(c.C)%*%c.H.u%*%c.ISxz
            c.pC = c.C%*%c.tC
            c.IpC = c.In-c.pC
            c.Mxz = c.pC + c.Sxz%*%c.IpC
            if(smooth){
                c.gf = c.X%*%(c.beta.s[(ncol(XM)-(c.Vdim-1)):ncol(XM)]) + c.Zs%*%c.r
            } else {
                c.gf = c.X%*%(c.beta.s[(ncol(XM)-(c.Vdim-1)):ncol(XM)])
            }
    	} else {
            c.X <- XM
            if(smooth){
                c.pZ = c.Zs%*%c.tZ
                c.IpZ = c.In - c.pZ
            } else {
                c.pZ = NULL
                c.IpZ = c.In
            }
            c.XHX = t(c.X)%*%c.H.u%*%c.IpZ%*%c.X
            c.XHXi = .whichInverse(use$inverse, c.XHX, tol = use$tol, Re.values = TRUE)
            c.tX = c.XHXi%*%t(c.X)%*%c.H.u%*%c.IpZ
            c.pX = c.X%*%c.tX
            c.IpX <- c.In - c.pX
            if(smooth){
                c.Sxz = c.pX + c.pZ%*%c.IpX
            } else {
                c.Sxz <- c.pX
            }
            c.Mxz = c.Sxz
            c.ISxz = c.In - c.Sxz
            c.gf = c.yf
    	}
    } else {
    	c.Mxz = XM%*%c.tXs
        c.Sxz = NULL
        c.gf = NULL
    }
	c.epsilon = zV - XM%*%c.beta.s
	c.RSS = t(c.epsilon)%*%c.LLi%*%(c.epsilon)
	if(!use.reml) c.df = nrow(XM) else c.df = (nrow(XM)-ncol(XM))
	c.sigma2s = c.RSS/c.df
    c.tau2 <- (1-c.rhos)*c.sigma2s
    c.sigma2 = c.sigma2s - c.tau2
    c.Sigma.u = c.sigma2s%x%c.R.u
    colnames(c.XLX) <- colnames(XM)
    rownames(c.XLX) <- colnames(XM)
    colnames(c.XLXi) <- colnames(XM)
    rownames(c.XLXi) <- colnames(XM)
    names(c.beta.s) <- colnames(XM)
	list(alpha = RZ$alpha, sigma2s = c.sigma2s, sigma2 = c.sigma2, phi = RZ$phi, kappa = RZ$kappa, tau2 = c.tau2, rhos = c.rhos, log.rhos = c.log.rhos, beta = c.beta.s, r = c.r, psi = IP$psi, history.psi = IP$history, R = c.R.u, g = if(exists("c.gf")){c.gf}, fitted.values = if(exists("c.yf")){c.yf},	S = if(exists("c.Sxz")){c.Sxz}, M = if(exists("c.Mxz")){c.Mxz}, X = XM, Z = if(exists("c.Z")){c.Z}, Zs = if(exists("c.Zs")){c.Zs}, W = RZ$W, XLX = c.XLX, XLXi = c.XLXi, var.beta = c.sigma2s%x%c.XLXi, conf.error = sqrt(diag(c.Mxz%*%c.Sigma.u%*%t(c.Mxz))), pred.error = sqrt(diag(c.Sigma.u + c.Mxz%*%c.Sigma.u%*%t(c.Mxz))), penalty = if(use$penalty!="none"){list(P = c.Ps, decompose = c.l.Ps)})
}
.scpInitialize = function(zV, XM, ZM, VM, hM, Qel, W, ch, use, pars, initial, method,
                       control.optim, control.gs, control.ch, print.pars){
    if(missing(use)){
        use = list(penalty = "none", is.X = "none", ps.order = 2, Ztr = "RV",
             inverse = "solve", tol = .Machine$double.neg.eps,
             Amethod = "eigen", Adiag = TRUE)
    } else {
        if(is.null(use)){
            stop("List use must be specified!")
        } else {
            if(is.null(use$penalty)){
                use$penaly = "none"
            }
            if(use$penalty=="cs" | use$penalty=="ps"){
                use$is.X = "tensor"
            }
            if(use$penalty=="tps"){
                use$is.X = "tps"
            }
            if(is.null(use$is.X) & use$penalty=="none"){
                use$is.X = "none"
            }
            if(is.null(use$ps.order) & use$penalty=="ps"){
                use$ps.order = 2
            }
            if(is.null(use$Ztr)){
                use$Ztr = "RV"
            }
            if(is.null(use$inverse)){
                use$inverse = "solve"
            }
            if(is.null(use$tol)){
                use$tol = .Machine$double.neg.eps
            }
            if(is.null(use$Amethod)){
                use$Amethod = "eigen"
            }
            if(is.null(use$Adiag)){
                use$Adiag = ifelse(use$penalty=="ps", TRUE, TRUE)
            }
        }
    }
    if(missing(pars)){
        pars = list(fix.nugget = NULL, fix.kappa = NULL, is.log.rhos = TRUE,
                  model = "matern", nugget.tol = 1.0e-15)
    } else {
        if(is.null(pars)){
            stop("List pars must be specified!")
        } else {
            if(is.null(pars$model)){
                pars$model = "matern"
            }
            if(is.null(pars$nugget.tol)){
                pars$nugget.tol = 1.0e-15
            }
        }
    }
    use.gs = FALSE
    alt.gs = FALSE
    if(missing(initial)){
        use.gs = TRUE
        initial = NULL
    }
    if(is.null(initial)) message("Initial values not specified. Using internal search!")
        if(missing(method)){
            method = list(use.reml = FALSE, use.profile = TRUE, linearize = "eigen")
        } else {
        if(is.null(method)){
            stop("List method must be specified!")
        } else {
            if(is.null(method$use.reml)){
                method$use.reml = FALSE
            }
            if(is.null(method$use.profile)){
                method$use.profile = TRUE
            }
            if(is.null(method$linearize)){
                method$linearize = "eigen"
            }
        }
    }
    if(missing(control.optim)){
        control.optim = list(trace = 5, maxit = 300)
    } else {
        if(is.null(control.optim)){
            control.optim = list()
        }
    }
    if(missing(control.gs)){
        control.gs = list(is.zero = 1.0e-6, is.inf = 10000, size.each.grid = 30, size.sample = 40, num.attempts = 5)
    } else {
        if(is.null(control.gs)){
            stop("List control.gs must be specified!")
        } else {
            if(is.null(control.gs$is.zero)){
                control.gs$is.zero = 1.0e-6
            }
            if(is.null(control.gs$is.inf)){
                control.gs$is.inf = 10000
            }
            if(is.null(control.gs$size.each.grid)){
                control.gs$size.each.grid = 30
            }
            if(is.null(control.gs$size.sample)){
                control.gs$size.sample = 40
            }
            if(is.null(control.gs$num.attempts)){
                control.gs$num.attempts = 5
            }
            if(is.null(control.gs$replace.optim)){
                control.gs$replace.optim = FALSE
            }
        }
    }
    if(missing(ch)){
        control.ch = list(maxiter = 20, maxskip = 10)
    } else {
        if(is.null(control.ch)){
            stop("List control.ch must be specified!")
        } else {
            if(is.null(control.ch$maxiter)){
                control.ch$maxiter = 20
            }
            if(is.null(control.ch$maxskip)){
                control.ch$maxskip = 10
            }
        }
    }
    if(missing(print.pars)){ print.pars = FALSE }
    use$Adiag = ifelse(use$penalty=="ps", TRUE, use$Adiag)
    cycle <- 1
    while(cycle <= control.gs$num.attempts){
        cat(".")
		message(paste("Computing: cycle",cycle))
        c.initial = unlist(initial)
        c.theta = .getBounds(pars = pars, smooth.type = use$penalty)
        c.req = names(c.theta$lower)
        c.start = c.theta
        c.id.el = match(c.req,names(c.initial))
        if(all(!is.na(c.id.el))){
            c.start$pars = c.initial[c.id.el]
        } else {
            use.gs = TRUE
            c.start$pars = c.initial[!is.na(match(names(c.initial),c.req))]
        }
        if(!is.null(initial)){
            if(any(grepl("beta",names(initial)))){
                if(all(!is.na(match(names(initial$beta),paste("beta",1:ncol(XM),sep=""))))){
                    pars$addbeta = TRUE
                    c.start$pars = c(initial$beta,c.start$pars)
                    c.start$lower = c(beta=rep(-Inf,ncol(XM)),c.start$lower)
                    c.start$upper = c(beta=rep(+Inf,ncol(XM)),c.start$upper)
                } else {
                    message("Length of starting values for beta does not match the required!\nUsing internal search for starting values not given!")
                    pars$addbeta = TRUE
                    use.gs = TRUE
                }
            } else {
                pars$addbeta = FALSE
            }
        } else {
            pars$addbeta = FALSE
        }
        if(use.gs | control.gs$replace.optim){
			message("One or more starting values were not provided.\nUsing internal search!")
            c.gstart = .gs(zV = zV, XM = XM, ZM = ZM, VM = VM, hM = hM, Qel = Qel, W = W,
                                    use = use, pars = pars,
                                    method = method, control = control.gs,
                                    print.pars = FALSE)
            c.gs = c.gstart
            if(use.gs){
                c.id.given = match(names(c.gstart$pars),names(c.start$pars))
                if(sum(!is.na(c.id.given))>0){
                    c.gstart$pars[!is.na(c.id.given)] = c.start$pars[c.id.given[!is.na(c.id.given)]]
                } else {
                    alt.gs = TRUE
                }
                c.start = c.gstart
            }
        }
        message(paste("Estimating ",paste(names(c.start$pars),sep="",collapse=", "),sep="",collapse=""))
        message(paste("Starting at ",paste(round(c.start$pars,5),sep="",collapse=", "),sep="",collapse=""))
        requireNamespace("stats", quietly = TRUE)
        message(paste(ifelse(method$use.reml,"REML","ML"),"estimation.",ifelse(method$use.profile,"Use profiling.\n","Not profiling.\n"),sep=" "))
        aux = .tryCatchWE(
            stats::optim(
            par = c.start$pars,
            .scpLoglik,
            gr = NULL,
            beta.par = pars$addbeta,
            is.log.rhos = pars$is.log.rhos,
            obs.z = zV,
            dist.m = hM,
            W = W,
            Xs = XM,
			Z = ZM,
			V = VM,
            fix.nugget = pars$fix.nugget,
            fix.kappa = pars$fix.kappa,
            pty = use$penalty,
            Q.base = Qel,
            Ztransf = use$Ztr,
            model = pars$model,
            ml = !method$use.reml,
            use.profile = method$use.profile,
            inverse.method = use$inverse,
            itol = use$tol,
            createA = use$Amethod,
            diagA = use$Adiag,
            nugget.tol = pars$nugget.tol,
            print.pars = print.pars,
            method = "L-BFGS-B", control = control.optim,
            hessian = TRUE, lower = c.start$lower, upper = c.start$upper
            )
        )
        cycle <- cycle + 1
        if(!aux$error){
            if(aux$value$convergence==0){
                cycle <- control.gs$num.attempts + 1
            }
        }
    }
    if(!aux$error){
        c.fit = aux$value
    } else {
        if((use.gs & alt.gs) | control.gs$replace.optim){
            c.fit = list()
            if(use.gs & alt.gs){
                c.fit$par = c.start$pars
            }
            if(control.gs$replace.optim){
                c.fit$par = c.gs$pars
            }
			message("Computing the Hessian")
            c.fit$hessian = stats::optimHess(
                par = c.fit$par,
                .scpLoglik,
                gr = NULL,
                beta.par = pars$addbeta,
                is.log.rhos = pars$is.log.rhos,
                obs.z = zV,
                dist.m = hM,
                W = W,
                Xs = XM,
				Z = ZM,
				V = VM,
                fix.nugget = pars$fix.nugget,
                fix.kappa = pars$fix.kappa,
                pty = use$penalty,
                Q.base = Qel,
                Ztransf = use$Ztr,
                model = pars$model,
                ml = !method$use.reml,
                use.profile = method$use.profile,
                inverse.method = use$inverse,
                itol = use$tol,
                createA = use$Amethod,
                diagA = use$Adiag,
                nugget.tol = pars$nugget.tol,
                print.pars = print.pars,
                control = control.optim
            )
        } else {
            message("Error in optim procedure!")
            stop(aux$value)
        }
    }
    if(is.null(ch)){
		message("Obtaining mean estimates")
        c.out = .betaFromPars(pars = c.fit$par, zV = zV, XM = XM, ZM = ZM, VM = VM, hM = hM, Qel = Qel, W = W,
                                 use = use, pars.str = pars, use.reml = method$use.reml, only.beta = FALSE)
    } else {
		message("Obtaining mean estimates")
        c.elem = .RSPZfromPars(pars = c.fit$par, zV = zV, XM = XM, ZM = ZM, VM = VM, hM = hM, Qel = Qel, W = W,
                                  use = use, pars.str = pars)
		message("Obtaining change-points")
        c.psi.l = .psiIterative(cp.call = ch$cp.call, psi0 = ch$psi0, UVcols = ch$UVcols,
                                  obj = c.elem, control.cp = control.ch)
		message("Obtaining fitted model")
        c.out = .fitFromRSPZ(RZ = c.elem, IP = c.psi.l, use.reml = method$use.reml)
    }
    if(!aux$error){
        c.out$optim = c.fit
    } else {
        c.out$gs = c.fit
    }
	cat("\n")
    return(c.out)
}
landim1 <- load(ifelse(file.exists("data/landim1.rda"),"data/landim1.rda",file.path(getwd(), "data/landim1.rda")))
