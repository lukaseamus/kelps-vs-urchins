library(nnls)
# Gaussian peak coefficient file, all pigments
gaussian.peaks <- read.table("gaussian.peak.parameters.txt", header = TRUE) 

# Specific absorption coefficient (sac) file including alias pigments (zea, cryp, anth)
sac.table <- read.table("specific absorption coefficients.txt", header =T)
sac <- sac.table$L.g.cm
names(sac) <- sac.table$pigment

# Individual Gaussian peak spectrum as function of wavelength (nm):
#
#	a exp(-((wavelength - xp - m.bias) / (s.bias * s))^2 / 2)
#
# where xp is peak wavelength (nm)
#	   s is peak halfwidth (nm), and
#	   a is peak weight (dimensionless)

gaussian.peak.function <- function(g, x, m.bias=0, s.bias=1) {

	peaks <- function(h) {
		h["a"] * exp(-0.5 * ((x - h["xp"] - m.bias) / (s.bias * h["s"]))^2)
	}

	# Gaussian spectral components
	p <- apply(g, 1, peaks)	

	return(apply(p, 1, sum))
}


# core <- c("Chl.a", "Chl.b", "Chl.c1", "Chl.c2", "Phe.a", "Phe.b", "Allo", "bb.Car",  
# 	"Cantha", "Diadino", "Diato", "Dino", "Echin", "Fuco", "Lut", "Myxo",     
# 	"c.Neo", "Peri", "Viola")

pigment.basis <- function(w, pigments=core, m.bias= 0, s.bias=1) {

	# Make Gaussian peak spectra
	X.list <- by(gaussian.peaks[,-1], gaussian.peaks[, 1], 
		function(x) { gaussian.peak.function(x, w, m.bias, s.bias) })

	X.pig <- matrix(unlist(X.list), ncol = length(unique(gaussian.peaks$pigment)))
	colnames(X.pig) <- names(X.list)

	# NULL means make basis of the whole pigment list
	if (length(pigments) > 0) {		
		X.pig <- X.pig[, colnames(X.pig) %in% pigments]
	}

	return(X.pig)
}



# w <- 400:700
# X <- pigment.basis(w)
# matplot(w, X, type="l", lty=1)

background.basis <- function(w, k=6) {

	# Polynomial background basis 
	# (decreasing to keep non-negativity)

	u <- (max(w) - w) / (max(w) - min(w))

	X.bg <- rep(1, length(w))
	for (i in 1:k) {
		X.bg <- cbind(X.bg, u^i)
	}

	return(X.bg)
}

# w <- 400:700
# X <- background.basis(w)
# matplot(w, X, type="l", lty=1)

# Estimating pigment weigths using NNLS
pigment.fit <- function(w, y, pigments=core, k=6, m.bias= 0, s.bias= 1) {

	# Basis for fitting
	X <- cbind(background.basis(w, k), pigment.basis(w, pigments, m.bias, s.bias))

	# Fit each spectrum (cloumn of y) by nnls (non-negative least squares)
	
	if (is.vector(y) == TRUE) {
	m <- vector("list", 1)
		m <- nnls(X, y)
	
	} else {

	m <- vector("list", ncol(y))
	for (i in 1:ncol(y)) {
		m[[i]] <- nnls(X, y[[i]])
	
	}
}

	return(list(m=m, k=k, y=y, X=X, w=w, p=pigments))
}

# Calculating fitted spectrum (background + pigment)
fitted.spectrum <- function(m) { 	
	
	if (is.vector(m$y) == TRUE) {
	y.fit <- fitted(m$m)

	} else {
		
	y.fit <- sapply(m$m, fitted)
	colnames(y.fit) <- colnames(m$y)
}
	return(y.fit)

}

# Calculating fitted background absorbance spectra
background.spectrum <- function(m) { 	
	
	if (is.vector(m$y) == TRUE) {
	coefficients <- coef(m$m)

	# Row 1:k+1 is background	
	bg.coefficients <- coefficients[1:(m$k + 1)] 
	bg.fitted.spectrum <- m$X[ , 1:(m$k + 1)] %*% bg.coefficients 

	
	} else {
	
	coefficients <- sapply(m$m, coef)	
	bg.coefficients <- coefficients[1:(m$k + 1), ] 
	bg.fitted.spectrum <- m$X[ ,1:(m$k + 1)] %*% bg.coefficients 
	
	colnames(bg.fitted.spectrum) <- colnames(m$y)
}
	return(bg.fitted.spectrum)
}

# Calculating fitted pigment absorbance spectra
pigment.spectrum <- function(m) {		
	
	if (is.vector(m$y) == TRUE) {
	coefficients <- coef(m$m)

	# Row k+2 to last row are pigments	
	pigment.coefficients <- coefficients[(m$k + 2):length(coefficients)] 	
	pigment.fitted.spectrum <- m$X[ , (m$k + 2):ncol(m$X)] %*% pigment.coefficients 
	
	} else {

	coefficients <- sapply(m$m, coef)

	# Row k+2 to last row are pigments	
	pigment.coefficients <- coefficients[(m$k + 2):nrow(coefficients), ] 	
	pigment.fitted.spectrum <- m$X[ , (m$k + 2):ncol(m$X)] %*% pigment.coefficients 
	
	colnames(pigment.fitted.spectrum) <- colnames(m$y)
}
	return(pigment.fitted.spectrum)

}

# Calculating pigment concentrations (mg/L) in extract
pigment.concentration <- function(m, pathl = 1) {

	if (is.vector(m$y) == TRUE) {
	coefficients <- coef(m$m)
	
	# divide by optical pathlength to get absorption 
	pigment.absorption <- coefficients[(m$k + 2):length(coefficients)] / pathl

	pigment.absorption <- pigment.absorption / apply((m$X[, (m$k + 2):length(coefficients)]), 2, max)

	# Subset the specific absorption coefficient table if pigment subset is used	
	if (length(m$p) > 0){
		sac <- sac[names(sac) %in% m$p] 
	}

	mg.L.extract <- pigment.absorption  * (1000 / sac[colnames(m$X)[(m$k + 2):length(coefficients)]])
	
	} else {

coefficients <- sapply(m$m, coef)
	
	# divide by optical pathlength to get absorption 
	pigment.absorption <- coefficients[(m$k + 2):nrow(coefficients), ] / pathl 

	pigment.absorption <- sweep(pigment.absorption, 1, 
		apply((m$X[, (m$k + 2):nrow(coefficients)]), 2, max), "/")


	rownames(pigment.absorption) <- colnames(m$X)[(m$k + 2):nrow(coefficients)]
	colnames(pigment.absorption) <- colnames(m$y)
	
	# Subset the specific absorption coefficient table if pigment subset is used	
	if (length(m$p) > 0){
		sac <- sac[names(sac) %in% m$p] 
	}

	mg.L.extract <- sweep(pigment.absorption, 1, 1000 / sac[colnames(m$X)[(m$k + 2):nrow(coefficients)]], "*")
}
	return(mg.L.extract)	
} 


# Example 
#w <- 400:700
#y <- read.table("COMSAT abs spec.txt", header = T)
#y <- subset(y, wl >= 400 & wl <= 700)[,-1]

# M <- pigment.fit(w, y)
# names(M)

#fit.spec <- fitted.spectrum(M)
#matplot(w, fit.spec, type = "l")

#bg.fitted.spec <- background.spectrum(M)
#matplot(w, bg.fitted.spec, type = "l")

#pigment.fitted.spec <- pigment.spectrum(M)
#matplot(w, pigment.fitted.spec, type = "l")

#pigment.concentration(M, pathl = 1)



