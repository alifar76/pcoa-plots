rm(list=ls())
require(grid)
require(vegan)
require(ggplot2)
require(gridExtra)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


dataframe_filter <- function(namearrange,meta,MYdata){
	rearrlist <- lapply(namearrange,function(x) eval(parse(text = paste("subset(MYdata,",meta, "==x)"))))
	newdataf <- do.call("rbind", rearrlist)
	return (newdataf)
}



PCOA_calcs <- function(dminp,mapfile,var,plottitle,colours){
	dmat <- read.table(dminp,header = T, sep = "\t", check.names = F, comment.char= "",row.names=1)
	MYmetaEF <- read.table(mapfile,header = T, sep = "\t", check.names = T,comment.char= "")
	namearrange <- colnames(dmat)
	meta <- "X.SampleID"
	MYmetaEF <- dataframe_filter(namearrange,meta,MYmetaEF)
	set.seed(123)
	MYpcoa <- cmdscale(as.dist(dmat), k = 2, eig=TRUE, add=TRUE) # Earlier had used dist function instead of vegdist, so was getting different values. dist did not support bray. 
	eig2 <- eigenvals(MYpcoa)
	print(eig2 / sum(eig2)[1])
	pc1 <- round((eig2 / sum(eig2)[1])[1]*100,2)
	pc2 <- round((eig2 / sum(eig2)[1])[2]*100,2)
	xaxis <- paste("PCoA1"," - ",pc1,"%",sep="")
	yaxis <- paste("PCoA2"," - ",pc2,"%",sep="")
	set.seed(123)
	adon <- tryCatch(adonis(as.dist(dmat) ~ MYmetaEF[,var], data=MYmetaEF, permutations = 1000, autotransform=T),error=function(e) NA)
	pvalue <- ifelse(is.na(adon$aov.tab[1,5]),"p-value = NA, variation = NA",paste("p-value = ",round(adon$aov.tab[1,6],3),", variation = ",round(adon$aov.tab[1,5],3)*100,"%",sep=""))
	PCOA <- data.frame(MYpcoa$points[,1],MYpcoa$points[,2], MYmetaEF[,var][!is.na(MYmetaEF[,var])])
	colnames(PCOA) <- c("PCoA1", "PCoA2" ,"comparison_groups")
	print (PCOA)
	plot.new()
  	#this is used to generate the confidence circles around the 95% confidence intervals
  	set.seed(123)
	ord<-ordiellipse(MYpcoa, PCOA$comparison_groups, display = "sites", kind = "se", conf = 0.95, label = T)
	df_ell <- data.frame()
	for(g in levels(PCOA$comparison_groups)){
		if (table(PCOA$comparison_groups)[[g]] > 2){
			df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCOA[PCOA$comparison_groups==g,],
      		vegan:::veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale))),comparison_groups=g))
    	}
	}
	plab <- gsub("^\\s+|\\s+$", "", strsplit(pvalue, ",")[[1]][1])
	varlab <- gsub("^\\s+|\\s+$", "", strsplit(pvalue, ",")[[1]][2])
	ppos <- max(MYpcoa$points[,2])
	varpos <- ppos - 0.025
	xpos <- min(MYpcoa$points[,1])+0.05
	PCOA <- cbind(PCOA,rownames(PCOA))
	colnames(PCOA) <- c(colnames(PCOA)[1:length(colnames(PCOA))-1],"samplenames")
	print (PCOA)
	# labs(color="Clinical Outcome") changes legend label. 
	p1 <- ggplot(data = PCOA, aes(x=PCoA1, y=PCoA2,)) + 
	geom_point(aes(color = comparison_groups),size=6) +				# Play around with size parameter to adjust dot sizes in PCoA plot
	scale_color_manual(values = colours)+ 
	geom_text(aes(label=samplenames),hjust=0, vjust=0) +			# Play-around with hjust and vjust values to your liking to best make plot
	xlab(xaxis) + 
	ylab(yaxis) + 
	geom_path(data=df_ell, aes(x=Dim1, y=Dim2,colour=comparison_groups), size=1, linetype=2) + 
	labs(title = plottitle, color="Phenotype") + 
	theme(plot.title = element_text(size = 25),axis.text=element_text(size=15,color="black"),axis.title=element_text(size=20,color="black"), legend.key.size = unit(1,"cm"),legend.text=element_text(size=15),legend.title=element_text(size=15)) + 
	annotate("text",x=c(xpos,xpos),y=c(ppos,varpos),label=c(plab,varlab))
	return (p1)
}



## Change these values:
setwd('/Users/alifaruqi/Desktop/Projects/Development_Tools/Github_Scripts/pcoa-plots/src')					# Working directory
outputname <- "beta_diversity_figures"																		# Name of output figure
dminp <- c("bray_curtis_dm.txt","canberra_dm.txt","unweighted_unifrac_dm.txt","weighted_unifrac_dm.txt")  # Vector of input DM names
## The Distance Matrices input names and plot titles must have same indices in their corresponding dimp and plottitle vector. 
compname <- "Abnormal/Normal"
mapfile <- "mapfile.txt"			# Name of mapping file
var <- "Status"						# Name of metadata variable on which adonis calculations are to be performed
colours <- c("yellow","blue")			# Colour for dots. If more than two level of metavariable, passs on more colours


### These function calls do not change
dmtest <- c("Comparison (Bray-Curtis)","Comparison (Canberra)","Comparison (Unweighted UniFrac)","Comparison (Weighted UniFrac)")  # Plot titles for each of DM used
plottitle <- sapply(dmtest,function(x){
    paste(compname,x,sep=" ")
   })
plottile <- as.vector(plottitle)
b1 <- PCOA_calcs(dminp[1],mapfile,var,plottitle[1],colours)
b2 <- PCOA_calcs(dminp[2],mapfile,var,plottitle[2],colours)
b3 <- PCOA_calcs(dminp[3],mapfile,var,plottitle[3],colours)
b4 <- PCOA_calcs(dminp[4],mapfile,var,plottitle[4],colours)

pdf(paste(outputname,".pdf",sep=""),height=15, width=20, useDingbats=FALSE)
multiplot(b1,b2,b3,b4,cols=2)
dev.off()
