# Rscript TOM_Vegan-Ordination.R (clean)

# Example Run : 
# Rscript TOM_Vegan-Ordination.R --slave --inomics=input.xlsx --osheet=genera --inpheno=input.xlsx --psheet=pheno --contrast=#BBC_BaCo --log=1 --output=ResultsDir


library(factoextra)
library(ggplot2)
library(pracma)
library(readxl)
library(vegan)

# Part 1 : command line argument parsing
	
	run.arguments <- commandArgs(TRUE)
	
	valid.run.parameters <- c( "inomics", "osheet", "inpheno", "psheet", "contrast", "log", "output" )
	
	for ( i in 1:length( run.arguments ) ) {

		if ( strcmpi( substr( run.arguments[i], 1, 2 ), "--" ) & grepl( "=", run.arguments[i], fixed = TRUE) ) {
			key.pair <- strsplit( run.arguments[i], split="=")
			run.parameter <- gsub( "--", "", key.pair[[1]][1] )
			run.argument <- key.pair[[1]][2]
			
			if ( ( run.parameter %in% valid.run.parameters ) & ( run.parameter == "inomics" ) ) {
				inomics <- run.argument
			}
			if ( ( run.parameter %in% valid.run.parameters ) & ( run.parameter == "osheet" ) ) {
				osheet <- run.argument
			}
			if ( ( run.parameter %in% valid.run.parameters ) & ( run.parameter == "inpheno" ) ) {
				inpheno <- run.argument
			}
			if ( ( run.parameter %in% valid.run.parameters ) & ( run.parameter == "psheet" ) ) {
				psheet <- run.argument
			}
			if ( ( run.parameter %in% valid.run.parameters ) & ( run.parameter == "contrast" ) ) {
				contrast <- run.argument
			}
			if ( ( run.parameter %in% valid.run.parameters ) & ( run.parameter == "log" ) ) {
				logtransform <- run.argument
			}			
			if ( ( run.parameter %in% valid.run.parameters ) & ( run.parameter == "output" ) ) {
				outdir <- run.argument
			}
		}
	}


# Part 2 : retrieve and transform data

	# Import from Excel ( two different data sheets : with samples in rows, and characteristics or features in columns )
	# It is important that samples in both data sheets are in the exact identical order !


	# Import omics df
	# output : features in rows, samples in columns

	omics <- read_excel(inomics, sheet = osheet) # imports as a Tibble, but these don't allow for row names [ https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html ]
	omics <- as.data.frame(omics) # convert Tibble to conventional DataFrame
	
	rownames(omics) <- omics[,1] # get the row names
	omics[,1] <- NULL # remove the row name column

	
	# Import metadata df
	# output : samples in rows, characteristics in columns

	meta <- read_excel(inpheno, sheet = psheet) # imports as a Tibble, but these don't allow for row names [ https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html ]
	meta <- as.data.frame(meta) # convert Tibble to conventional DataFrame
	
	rownames(meta) <- meta[,1] # get the row names
	meta[,1] <- NULL # remove the row name column


	# Transformations of omics df
	# log-transform, center
	# Y' = log ( A * Y + 1 ) ; where A is the 'strength' of the log transformation : 1, 10, 100, 1000, etc., default = 1

	logtransform = as.integer(logtransform)
	omics.log <- ( logtransform * omics ) + 1 
	omics.log <- log( omics.log ) 
	omics.sc <- scale(omics.log, center = TRUE, scale = FALSE) 


# Part 3 : ordination

	# PCA with Vegan
	# example on the omics df, with the 'Niche' column from the metadata df as phenotypic information 

	pca.out <- vegan::rda( omics.sc, scale = FALSE ) 
	pca.out.sum <- summary(pca.out)


	# RDA with Vegan
	# example on the omics df, with the 'Niche' column from the metadata df as explanatory variable
	
	rda.out <- vegan::rda( omics.sc ~ get(contrast, meta) + Condition(NULL), data = meta, scale = FALSE, na.action = na.fail, subset = NULL )	# previously : ~ `#BBCBaCo` +
	rda.out.sum <- summary(rda.out)
	
	# calculate RDA p-value by permutation test with Vegan 'permutest'
	
	rda.out.pval <- permutest ( rda.out , permutations = 999, pairwise = FALSE )
	
	
	# print SUM PCA/RDA results and RDA p-value to screen
	
	print ( "# SUMMARY Report of PCA/RDA results :" )
	cat("\n")
	print ( rda.out )
	print ( "# SUMMARY Report of RDA p-value :" )
	print ( rda.out.pval )
	cat("\n")
	
	# print SUM PCA/RDA results and RDA p-value to outfile
	
	sumoutfile <- paste0(outdir,"/PCA-RDA-contrast_",contrast,".screenprint.summary.log")
	sink(sumoutfile)
	print ( "# SUMMARY Report of PCA/RDA results :" )
	cat("\n")
	print ( rda.out )
	print ( "# SUMMARY Report of RDA p-value :" )
	print ( rda.out.pval )
	sink()
	
	# print ALL PCA/RDA results to outfile
	
	alloutfile <- paste0(outdir,"/PCA-RDA-contrast_",contrast,".screenprint.report.log")
	sink(alloutfile)
	print ( "# FULL Report of PCA/RDA results :" )
	cat("\n")
	print ( rda.out.sum )
	sink()
	
	
# Part 4: plotting graphs

	# Graphical print, from command line to device
	
	# Create groups

	mygroups <- as.numeric(as.factor( get(contrast, meta) )) # note that 'Niche' currently is hard-coded, idem for the number of groups and their colors !

	lab.mygroups = ( get(contrast, meta) )	# previously : lab.mygroups = ( meta$`#BBCBaCo` )
	
	pch.mygroups <- c(21, 22, 23, 24, 25)[mygroups]
	
	col.mygroups <- c("blue", "red", "green", "yellow", "pink")[mygroups]


mypdf <- paste0(outdir,"/PCA-RDA-contrast_",contrast,".plots.pdf")

pdf(file=mypdf)

	# Plot PCA biplot with Vegan in 'biplot'
	biplot(pca.out, type = 'n') # optional : display = "sites" or, display = "species"
	
		points(pca.out, display = "sites", cex=0.7, pch=pch.mygroups, col=col.mygroups, bg="white")
		# text(pca.out, display = "sites", cex=0.5, col="black")
	
		points(pca.out, display = "species", cex=0.7, pch=21, col="blue", bg="white")
		text(pca.out, display = "species", cex=0.5, col="blue")
		
		for ( i in unique (col.mygroups) ) ( ordihull(pca.out, groups = col.mygroups, show.group = i, col = i, draw = 'polygon', alpha = 50) ) # label = TRUE

		for ( i in unique (lab.mygroups) ) ( ordispider(pca.out, groups = col.mygroups, show.group = col.mygroups[ match(i,lab.mygroups) ], col = 'black', lty = 'blank', label = TRUE, labels = i) )

		taxloadings <- scores(pca.out, display = "species")
		arrowlen <- 0.9
		arrows(0, 0, arrowlen * taxloadings[, 1],  arrowlen * taxloadings[, 2], length = 0.05, col = "darkcyan")
		
	
	# Plot RDA triplot with Vegan in 'ordiplot'
	ordiplot(rda.out, type = 'n') |> # optional : display = "sites" or, display = "species"
	
		points("sites", cex=0.7, pch=pch.mygroups, col=col.mygroups, bg="white") |>
		# text("sites", cex=0.5, col="black") |>
	
		points("species", cex=0.7, pch=21, col="blue", bg="white") |>
		text("species", cex=0.5, col="blue", arrows = TRUE)

		for ( i in unique (col.mygroups) ) ( ordihull(rda.out, groups = col.mygroups, show.group = i, col = i, draw = 'polygon', alpha = 50) ) # label = TRUE
		
		for ( i in unique (lab.mygroups) ) ( ordispider(rda.out, groups = col.mygroups, show.group = col.mygroups[ match(i,lab.mygroups) ], col = 'black', lty = 'blank', label = TRUE, labels = i) )
		
dev.off()

# run from command line with : Rscript RDA.R --slave
