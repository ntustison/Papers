library(MASS)
library(bestglm)
md<-read.csv('meandiffusion.csv')
vol<-read.csv('volume.csv')
th<-read.csv('thickness.csv')
demog<-read.csv('demographics.csv')
# remove some bad data
whold<-( demog$AGE > 50 & demog$AGE > 1  )
whmid<-( demog$AGE > 40 & demog$AGE < 60  )
whyoung<-( demog$AGE < 41 & demog$AGE > 1  )
wh<-whold
print(sum(wh)/nrow(demog))
demog<-demog[wh,]
vol<-vol[wh,]
th<-th[wh,]
md<-md[wh,]
# training / testing data
subs1<-c(1:(nrow(demog)/2))
subs1<-subs1*2
subs2<-subs1-1
subs2<-subs2[1:(length(subs2)-1)] ; length(subs2) ; length(subs1)
# make data frame with all training data
Xy1<-data.frame(vol=vol[subs1,],th=th[subs1,],sex=demog$SEX[subs1],age=demog$AGE[subs1])
Xy2<-data.frame(vol=vol[subs2,],th=th[subs2,],sex=demog$SEX[subs2],age=demog$AGE[subs2])
Xy1<-data.frame(th=th[subs1,],sex=demog$SEX[subs1],age=demog$AGE[subs1])
names(Xy1)<-c(paste("X",1:(ncol(Xy1)-1),sep=""),"y")
Xy2<-data.frame(th=th[subs2,],sex=demog$SEX[subs2],age=demog$AGE[subs2])
names(Xy2)<-c(paste("X",1:(ncol(Xy2)-1),sep=""),"y")
Xy<-Xy1
ans<-bestglm(Xy, IC="CV" , t=125 )
# ans<-bestglm(Xy, IC="CV", t=10,CVArgs=list(Method="d", K=2, REP=10))
# now apply to test data
Xy<-Xy2
predictedage<-predict(ans$BestModel,newdata=Xy)
realage<-Xy$y
print(paste('Mean:',mean(abs(predictedage-realage))))
ranger<-c( min( c(realage,predictedage) ),  max( c(realage,predictedage) ) )
pdf( file = "brian_results.pdf", 8, 8 );
plot( ranger, ranger, type = "n",
    xlab = expression( paste( bold( "Real Age (years)" ) ) ),
    ylab = expression( paste( bold( "Predicted Age (years)" ) ) ),
    frame.plot = FALSE  );
  title( main = paste( "Cortical thickness prediction of age "), col.main = "black", font.main = 2 );
  points( realage, predictedage, pch = 22, lwd = 1.5, cex = 0.75, col = "red" );
  axis( 1, tick = TRUE, lwd = 2, font = 2 );
  axis( 2, tick = TRUE, lwd = 2, font = 2 );
dev.off();
print(cor.test(predictedage,realage))
