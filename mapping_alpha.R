require(vegan)
require(segmented)
require(raster)
require(maps)
require(classInt)
require(RColorBrewer)
require(mgcv)

mapping <- read.csv('mapping_alpha.csv', header = T,
	row.names = 1, sep = ',')

attach(mapping)
fit1 <- lm(sunplin ~ PrecAnn)
summary(segmented(fit1, seg.Z = ~ PrecAnn, data = mapping))

low.MAP <- mapping[which(mapping$PrecAnn<=1480),]
hi.MAP <- mapping[which(mapping$PrecAnn>1480),]

##gls
attach(low.MAP)

###low.ac
f1 <- formula(PDmean ~ PrecAnn)
f1.null <- formula(PDmean ~ 1)

low.ac <- gls(f1, correlation = corExp(form = ~ 'Long10' + 'Lat10',
	nugget = TRUE), data = low.MAP)

vario1 <- Variogram(low.null, form = ~'Long10' + 'Lat10', resType = "pearson")
plot(vario1, smooth = T, ylim = c(0, 1.1))
vario2 <- Variogram(low.ac, form = ~'Long10' + 'Lat10',
	resType = 'normalized')
plot(vario2, smooth = F, ylim = c(0, 1.1))

cor(low.MAP$sesPD,predict(low.ac))^2 #0.4883892

##AIC low
#Spherical = 2958.449
#Linear = 2961.12
#Ratio = 2950.353
#Gaussian = 2980.033
#Exponential = 2930.691
#NULL = 3387.047
#delta exponential = -456.356

detach(low.MAP)
attach(hi.MAP)

hi.null <- gls(f1, data = hi.MAP)
hi.ac <- gls(f1, correlation = corExp(form = ~ 'Long10' + 'Lat10',
	nugget = TRUE), data = hi.MAP)

vario3 <- Variogram(hi.null, form = ~'Long10' + 'Lat10',
	resType = "pearson")
vario4 <- Variogram(hi.ac, form = ~'Long10' + 'Lat10',
	resType = 'normalized')

tiff('spatialAC1.tiff', height = 15, width = 15, unit = 'cm', res = 600)
plot(vario1, smooth = T, ylim = c(0, 1.3),
	xlab = 'Distance (km)', main = 'a')
dev.off()
tiff('spatialAC2.tiff', height = 15, width = 15, unit = 'cm', res = 600)
plot(vario2, smooth = F, ylim = c(0, 1.3),
	xlab = 'Distance (km)', main = 'b')
dev.off()
tiff('spatialAC3.tiff', height = 15, width = 15, unit = 'cm', res = 600)
plot(vario3, smooth = T, ylim = c(0, 1.3),
	xlab = 'Distance (km)', main = 'c')
dev.off()
tiff('spatialAC4.tiff', height = 15, width = 15, unit = 'cm', res = 600)
plot(vario4, smooth = F, ylim = c(0, 1.3),
	xlab = 'Distance (km)', main = 'd')
dev.off()

##AIC hi
#Spherical = 2344.146
#Linear = 2344.146
#Ratio = 2209.484
#Gaussian = 2215.874
#Exponential = 2204.931
#NULL = 2413.85
#delta exponential = -208.919

cor(hi.MAP$sesPD,predict(hi.ac))^2 #0.1084856

Vario1D <- Variogram(B1D, form = ~ Long10 + Lat10,
	robust = TRUE, maxDist = 5000, resType = "normalized")

plot(Vario1D, smooth = FALSE)

#maps
detach(hi.MAP)
attach(mapping)
ll <- cbind(Long10, Lat10)

tiff('mapping_alpha.tiff', height = 15, width = 30, unit = "cm", res = 600)
par(mfrow=c(1,2))

plot(log(SR) ~ PrecAnn, type = 'n',
	xlab = 'Mean Annual Precipitation (mm)',
	ylab = 'log(Species Richness)', main = '(a)',
	cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)

var <- SR
nclas <- brewer.pal(4, "Spectral")
nclas <- nclas[4:1]
class <- classIntervals(var)
colcode <- findColours(class, nclas)

points(PrecAnn, log(SR),
	pch = c(17, 16)[as.numeric(mapping$MAP)], col = colcode)
SR.lm <- summary(lm(log(SR) ~ poly(PrecAnn, 2, raw = T)))
SR.r <- bquote(r^2 == .(format(SR.lm$adj.r.squared, digits = 2)))
legend("topleft", legend = SR.r, bty = "n", cex = 1.5)
SR.poly <- lm(log(SR) ~ poly(PrecAnn, 2, raw = T))
myPredict <- predict(SR.poly, interval="confidence", level = 0.99)
ix <- sort(PrecAnn,index.return=T)$ix
polygon(c(rev(PrecAnn[ix]), PrecAnn[ix]), c(rev(myPredict[ ix,3]),
	myPredict[ ix,2]), border = NA, col = gray(0.2,alpha = 0.5))
curve.dat = data.frame(PrecAnn, predict(SR.poly))
curve.dat = curve.dat[order(curve.dat$PrecAnn),]
lines(curve.dat, lwd = 5, col = gray(0.2,alpha = 0.9))


map(xlim = c(-77.8533,-34.8608), ylim = c(-24.9044,8.4711), lty = 'dashed')
map.axes(cex.axis = 1.5)
points(mapping$Long10, mapping$Lat10,
	pch = c(17, 16)[as.numeric(mapping$MAP)], col = colcode, cex = 0.5)
ordisurf(ll, mapping$PrecAnn, col = "black", main = "", add = T,
	level = c(1200, 1800, 2400, 2800), labcex = 0.8)
main = title('(b)', cex.main = 1.5)

dev.off()


#sesPD
tiff('mapping_sesPD_NEE.tiff', height = 15, width = 15, unit = "cm", res = 600)
#par(mfrow=c(1,2))

plot(mapping$sesPD ~ mapping$PrecAnn, col = "white",
	xlab = "Mean Annual Precipitation (mm)",
	ylab = "Lineage Diversity (sesPD)", cex.lab = 1.5, cex.axis = 1.5)
main = title('a', cex.main = 1.5, adj = 0)
low.MAP <- mapping[which(mapping$PrecAnn<1490),]
hi.MAP <- mapping[which(mapping$PrecAnn>=1490),]

var <- mapping$sesPD
nclas <- brewer.pal(4, "Spectral")
nclas <- nclas[4:1]
class <- classIntervals(var)
colcode <- findColours(class, nclas)

points(mapping$PrecAnn, mapping$sesPD,
	pch = c(17, 16)[as.numeric(mapping$MAP)], col = colcode)
low.r = bquote(r^2 == .(format(cor(low.MAP$sesPD,predict(low.ac))^2,
	digits = 2)))
hi.r = bquote(r^2 == .(format(cor(hi.MAP$sesPD,predict(hi.ac))^2,
	digits = 2)))
legend("topleft", legend = low.r, bty = "n", cex = 1.5)
legend("topright", legend = hi.r, bty = "n", cex = 1.5)
x_low <- c(min(mapping$PrecAnn), 1490)
y_low <- x_low* 0.0034780-5.3610
x_hi <- c(1490,max(mapping$PrecAnn))
y_hi <- x_hi*(-0.0006022)+ 0.7165
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

tiff('mapping_sesPD_NEE2.tiff', height = 15, width = 15, unit = "cm", res = 600)

map(xlim = c(-77.8533,-34.8608), ylim = c(-24.9044,8.4711), lty = 'dashed')
map.axes(cex.axis = 1.5)
points(mapping$Long10, mapping$Lat10,
	pch = c(17, 16)[as.numeric(mapping$MAP)], col = colcode, cex = 0.5)
ordisurf(ll, mapping$PrecAnn, col = "black", main = "", add = T,
	level = c(1200, 1800), labcex = 0.8)
ordisurf(ll, mapping$PrecAnn, col = "gray60", main = "", add = T,
	level = c(2400, 2800), labcex = 0.8)
main = title('b', cex.main = 1.5, adj = 0)

dev.off()


#endemism
tiff('mapping_endemism_final3.tiff', height = 15, width = 30, unit = "cm", res = 600)
par(mfrow=c(1,2))

plot(endemism ~ PrecAnn, log = 'y', type = 'n',
	xlab = 'Mean Annual Precipitation (mm)',
	ylab = 'Phylogenetic Endemism (myrs)',
	cex.lab = 1.5, cex.axis = 1.5)
main = title('c', cex.main = 1.5, adj = 0)

var <- endemism
nclas <- brewer.pal(4, "Spectral")
nclas <- nclas[4:1]
class <- classIntervals(var)
colcode <- findColours(class, nclas)

points(PrecAnn, endemism,
	pch = c(17, 16)[as.numeric(mapping$MAP)], col = colcode)
endemism.lm <- summary(lm(log(endemism) ~ poly(PrecAnn, 2, raw = T)))
endemism.r <- bquote(r^2 == .(format(endemism.lm$adj.r.squared, digits = 2)))
legend("topleft", legend = endemism.r, bty = "n", cex = 1.5)
endemism.poly <- lm(endemism ~ poly(PrecAnn, 2, raw = T))
myPredict <- predict(endemism.poly, interval="confidence", level = 0.99)
ix <- sort(PrecAnn,index.return=T)$ix
polygon(c(rev(PrecAnn[ix]), PrecAnn[ix]), c(rev(myPredict[ ix,3]),
	myPredict[ ix,2]), border = NA, col = gray(0.2,alpha = 0.5))
curve.dat = data.frame(PrecAnn, predict(endemism.poly))
curve.dat = curve.dat[order(curve.dat$PrecAnn),]
lines(curve.dat, lwd = 5, col = gray(0.2,alpha = 0.9))


map(xlim = c(-77.8533,-34.8608), ylim = c(-24.9044,8.4711), lty = 'dashed')
map.axes(cex.axis = 1.5)
points(mapping$Long10, mapping$Lat10,
	pch = c(17, 16)[as.numeric(mapping$MAP)], col = colcode, cex = 0.5)
ordisurf(ll, mapping$PrecAnn, col = "black", main = "", add = T,
	level = c(1200, 1800), labcex = 0.8)
ordisurf(ll, mapping$PrecAnn, col = "gray60", main = "", add = T,
	level = c(2400, 2800), labcex = 0.8)
main = title('d', cex.main = 1.5, adj = 0)

dev.off()

##supp maps
tiff('mapping_MAT2.tiff', height = 15, width = 30, unit = "cm", res = 600)
par(mfrow=c(1,2))
plot(sesPD ~ TempAnn, col = "white",
	xlab = "Mean Annual Temperature (°C)",
	ylab = "Lineage Diversity (sesPD)", cex.lab = 1.5, cex.axis = 1.5)
main = title('a', cex.main = 1.5, adj = 0)

var <- sesPD
nclas <- brewer.pal(4, "Spectral")
nclas <- nclas[4:1]
class <- classIntervals(var)
colcode <- findColours(class, nclas)

points(TempAnn, sesPD, pch = 16, col = colcode)
legend("bottomleft", 'ns', bty = "n", cex = 1.5)

map(xlim = c(-77.8533,-34.8608), ylim = c(-24.9044,8.4711), col = 'darkgray')
map.axes(cex.axis = 1.5)
points(mapping$Long10, mapping$Lat10, pch = 16, col = colcode, cex = 0.5)
ordisurf(ll, mapping$TempAnn, col = "black", main = "", add = T,
	levels = c(21, 23, 25, 26, 27), labcex = 0.8)
main = title('b', cex.main = 1.5, adj = 0)

dev.off()

tiff('mapping_TempMin2.tiff', height = 15, width = 30, unit = "cm", res = 600)

par(mfrow=c(1,2))
plot(sesPD ~ TempMin, col = "white",
	xlab = "Mean Minimum Temperature (°C)",
	ylab = "Lineage Diversity (sesPD)", cex.lab = 1.5, cex.axis = 1.5)
main = title('c', cex.main = 1.5, adj = 0)

var <- sesPD
nclas <- brewer.pal(4, "Spectral")
nclas <- nclas[4:1]
class <- classIntervals(var)
colcode <- findColours(class, nclas)

points(TempMin, sesPD, pch = 16, col = colcode)
legend("bottomleft", 'ns', bty = "n", cex = 1.5)

map(xlim = c(-77.8533,-34.8608), ylim = c(-24.9044,8.4711), col = 'darkgray')
map.axes(cex.axis = 1.5)
points(mapping$Long10, mapping$Lat10, pch = 16, col = colcode, cex = 0.5)
ordisurf(ll, mapping$TempMin, col = "black", main = "", add = T,
	nlevels = 5, labcex = 0.8)
main = title('d', cex.main = 1.5, adj = 0)

dev.off()

#PAs
require(raster)
require(maps)

PAs <- read.csv('PAs_final.csv', header = T, row.names = 1, sep = ',')
wdpa <- shapefile('wdpa_crop')
south.amer <- shapefile('amer_sul')
brasil <- shapefile('Brazil')
ucpi <- shapefile('UCPI')
ucus <- shapefile('UCUS')

#e <- c(-77.8533, -34.8608, -24.9044, 8.4711)
#wdpa <- crop(wdpa, e)

tiff('PAs_final.tiff', height = 15, width = 15, unit = 'cm', res = 600)

map(xlim = c(-77.8533, -34.8608), ylim = c(-24.9044, 8.4711), col = 'white')
map.axes()

plot(wdpa, col = 'gray90', border = 'gray90', add = T)
#plot(ucpi, col = 'gray90', border = 'gray90', add = T)
plot(south.amer, add = T)
#plot(brasil, add = T)
points(PAs$Long10, PAs$Lat10, pch = 16,
	col = c('black',
	rgb(t(col2rgb("blue"))/255, alpha = 0.5),
	rgb(t(col2rgb("red"))/255, alpha = 0.5))[as.numeric(PAs$PA)], cex = 0.5)
legend('topright', c('high LD unprotected', 'low LD unprotected',
	'protected'), pch = 16,
	col = c(rgb(t(col2rgb("red"))/255, alpha = 0.5),
	rgb(t(col2rgb("blue"))/255, alpha = 0.5), 'black'), bty = 'n')

dev.off()

#amazon
attach(amazon)
ll <- cbind(Long10, Lat10)

tiff('amazon_TempMin.tiff', height = 15, width = 30, unit = 'cm', res = 600)
par(mfrow = c(1, 2))
plot(sesPD ~ TempMin, type = 'n',
	xlab = 'Mean Mininum Temperature (°C)',
	ylab = 'Lineage Diversity (sesPD)', main = '(a)',
	cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)

var <- sesPD
nclas <- brewer.pal(4, "Spectral")
nclas <- nclas[4:1]
class <- classIntervals(var)
colcode <- findColours(class, nclas)

points(TempMin, sesPD, pch = 17, col = colcode)
TM.lm <- summary(lm(sesPD ~ TempMin))
TM.r <- bquote(r^2 == .(format(TM.lm$adj.r.squared, digits = 2)))
legend("bottomleft", legend = TM.r, bty = "n", cex = 1.5)
TM.lm <- lm(sesPD ~ TempMin)
myPredict <- predict(TM.lm, interval="confidence", level = 0.99)
ix <- sort(TempMin,index.return=T)$ix
polygon(c(rev(TempMin[ix]), TempMin[ix]), c(rev(myPredict[ ix,3]),
	myPredict[ ix,2]), border = NA, col = gray(0.2,alpha = 0.5))
curve.dat = data.frame(TempMin, predict(TM.lm))
curve.dat = curve.dat[order(curve.dat$TempMin),]
lines(curve.dat, lwd = 5, col = gray(0.2,alpha = 0.9))

map(xlim = c(min(Long10), max(Long10)),
	ylim = c(min(Lat10), (max(Lat10))), col = 'darkgray')
map.axes(cex.axis = 1.5)
points(Long10, Lat10, pch = 17, col = colcode, cex = 0.5)
ordisurf(ll, TempMin, col = "black", main = "", add = T,
	levels = c(15, 18, 20, 22), labcex = 0.8)
main = title('(b)', cex.main = 1.5)

dev.off()

#amazon

basal <- read.csv('mapping_amazon.csv', header = T, row.names = 1, sep = ',')
attach(basal)


tiff('longitude_basal.tiff', height = 30, width = 15, unit = "cm", res = 600)
par(mfrow=c(2,1))

plot(basal_proportion ~ Long10, pch = 16, col = gray(0.2,alpha = 0.5),
	xlab = 'Longitude (°)',
	ylab = 'Proportion of basal lineages',
	cex.lab = 1.5, cex.axis = 1.5)
main = title('a', cex.main = 1.5, adj = 0)

long.lm <- summary(lm(basal_proportion ~ poly(Long10, 2, raw = T)))
long.r <- bquote(r^2 == .(format(long.lm$adj.r.squared, digits = 2)))
legend("topright", legend = long.r, bty = "n", cex = 1.5)
long.poly <- lm(basal_proportion ~ poly(Long10, 2, raw = T))
myPredict <- predict(long.poly, interval="confidence", level = 0.99)
ix <- sort(Long10,index.return=T)$ix
polygon(c(rev(Long10[ix]), Long10[ix]), c(rev(myPredict[ ix,3]),
	myPredict[ ix,2]), border = NA, col = gray(0.2,alpha = 0.5))
curve.dat = data.frame(Long10, predict(long.poly))
curve.dat = curve.dat[order(curve.dat$Long10),]
lines(curve.dat, lwd = 5, col = gray(0.2,alpha = 0.9))

plot(sesPD ~ basal_proportion, pch = 16, col = gray(0.2,alpha = 0.5),
	xlab = 'Proportion of basal lineages',
	ylab = 'Lineage Diversity (sesPD)',
	cex.lab = 1.5, cex.axis = 1.5)
main = title('b', cex.main = 1.5, adj = 0)
legend('topright', legend = 'ns')

dev.off()

