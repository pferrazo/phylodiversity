###load libraries
require(classInt)
require(maps)
require(nlme)
require(picante)
require(raster)
require(RColorBrewer)
require(segmented)
require(vegan)

###load community (genus-by-site) matix
comm <- read.csv(file.choose(), header = T, row.names = 1, sep = ',')
###load phylogeny
phy <- read.tree('phy.tre')
##matching objects
tmp <- match.phylo.comm(phy, comm)
comm <- tmp$comm
phy <- tmp$phy

###compute lineage diversity
LD <- ses.pd(comm, phy, null.model = 'taxa.labels', runs = 999)

###load metadata (md) matrix comprising explanatory variables
#md and community matrix have the same row order 
md <- read.csv(file.choose(), header = T, row.names = 1, sep = ',', dec = '.')
md <- cbind(LD, md)

###regressions
#linear
summary(lm(md$LD ~ md$MAP))#LD = Lineage Diversity; MAP = Mean Annual Precipitation
#quadratic
summary(lm(md$LD ~ poly(md$MAP, 2, raw = T))
#piecewise
summary(segmented(lm(md$LD ~ md$MAP), seg.Z = ~ MAP, data = md))

###subset matrix following break-point in piecewise regression
low.MAP <- md[which(md$MAP<=1490),]
hi.MAP <- md[which(md$MAP>1490),]

###Generalised Least Squares (controlling for spatial AC)
##low.MAP
attach(low.MAP)
f1 <- formula(LD ~ MAP)

low.null <- gls(f1, data = low.MAP)
low.ac <- gls(f1, correlation = corExp(form = ~ 'Long10' + 'Lat10', nugget = TRUE), data = low.MAP)
#check available correlation structures at https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/corClasses.html

##hi.MAP
detach(low.MAP)
attach(hi.MAP)

hi.null <- gls(f1, data = hi.MAP)
hi.ac <- gls(f1, correlation = corExp(form = ~ 'Long10' + 'Lat10', nugget = TRUE), data = hi.MAP)
#extracting AIC values
AIC(low.null)
AIC(low.ac)
AIC(hi.null)
AIC(hi.ac)

###Figure 2
latlong <- cbind(md$Long10, md$Lat10)
par(mfrow=c(1,2))

plot(md$LD ~ md$MAP, col = "white", xlab = "Mean Annual Precipitation (mm)", ylab = "Lineage Diversity (sesPD)",, 
     	cex.axis = 1.5, cex.lab = 1.5)
main = title('a', cex.main = 1.5, adj = 0)

#setting color ramp
var <- md$LD
nclas <- brewer.pal(4, "Spectral")
nclas <- nclas[4:1]
class <- classIntervals(var)
colcode <- findColours(class, nclas)	

points(md$MAP, md$LD, pch = c(17, 16)[as.numeric(md$MAP_beta)], col = colcode)
#extracting coefficients of correlation for legend
low.r = bquote(r^2 == .(format(cor(low.MAP$LD,predict(low.ac))^2, digits = 2)))
hi.r = bquote(r^2 == .(format(cor(hi.MAP$LD,predict(hi.ac))^2, digits = 2)))
legend("topleft", legend = low.r, bty = "n", cex = 1.5)
legend("topright", legend = hi.r, bty = "n", cex = 1.5)

#extracting ranges for curves
x.low <- c(min(md$MAP), 1490)
y.low <- x.low* 0.0034780-5.3610
x.hi <- c(1490,max(md$MAP))
y.hi <- x.hi*(-0.0006022)+ 0.7165
low.ci <- lm(sesPD ~ PrecAnn, data = low.MAP)
myPredict <- predict(low.ci, interval="confidence", level = 0.99)
ix <- sort(low.MAP$PrecAnn,index.return=T)$ix
polygon(c(rev(low.MAP$PrecAnn[ix]), low.MAP$PrecAnn[ix]),
	c(rev(myPredict[ ix,3]),
	myPredict[ ix,2]), border = NA, col = gray(0.2,alpha = 0.5))
lines(x_low, y_low, lwd = 5, col = gray(0.2,alpha = 0.9))
hi.ci <- lm(sesPD ~ PrecAnn, data = hi.MAP)
myPredict <- predict(hi.ci, interval="confidence", level = 0.99)
ix <- sort(hi.MAP$PrecAnn,index.return=T)$ix
polygon(c(rev(hi.MAP$PrecAnn[ix]), hi.MAP$PrecAnn[ix]),
	c(rev(myPredict[ ix,3]),
	myPredict[ ix,2]), border = NA, col = gray(0.2,alpha = 0.5))
lines(x_hi, y_hi, lwd = 5, col = gray(0.2,alpha = 0.9))
abline(v = 1490, lty = 2, lwd = 5, col = gray(0.2,alpha = 0.9))

dev.off()

###Figure 3

PAs <- read.csv('PAs_final.csv', header = T, row.names = 1, sep = ',')
wdpa <- shapefile('wdpa_crop')
south.amer <- shapefile('amer_sul')
brasil <- shapefile('Brazil')

#e <- c(-77.8533, -34.8608, -24.9044, 8.4711)
#wdpa <- crop(wdpa, e)

map(xlim = c(-77.8533, -34.8608), ylim = c(-24.9044, 8.4711), col = 'white')
map.axes()

plot(wdpa, col = 'gray90', border = 'gray90', add = T)
plot(south.amer, add = T)
points(PAs$Long10, PAs$Lat10, pch = 16,
	col = c('black',
	rgb(t(col2rgb("blue"))/255, alpha = 0.5),
	rgb(t(col2rgb("red"))/255, alpha = 0.5))[as.numeric(PAs$PA)], cex = 0.5)
legend('topright', c('high LD unprotected', 'low LD unprotected',
	'protected'), pch = 16,
	col = c(rgb(t(col2rgb("red"))/255, alpha = 0.5),
	rgb(t(col2rgb("blue"))/255, alpha = 0.5), 'black'), bty = 'n')

dev.off()

