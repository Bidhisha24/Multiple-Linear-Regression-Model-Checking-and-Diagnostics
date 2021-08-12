# Read the data and construct the scatterplot matrix


library(readxl)
file="C:/Users/USER/Desktop/project/data.xlsx"
data1=read_excel(file)
n=nrow(data1)
p=ncol(data1)-2

# Check for linearity
correlation=cor(data1[,3],data1[4:(p+2)])
correlation


#Fit linear model to the data
data1.model=lm(LifeExpectancy~InfantDeaths+BMI+UnderFiveDeaths+Polio+Diphtheria+IncomeCompositionOfResources,data=data1)
data1.model




#Find if the model fits the data well
summary(data1.model)


# Check for model inadequacies

unit=seq(1,1,length.out=n)
x1=data1$InfantDeaths
x2=data1$BMI
x3=data1$UnderFiveDeaths
x4=data1$Polio
x5=data1$Diphtheria
x6=data1$IncomeCompositionOfResources
y=data1$LifeExpectancy
x=cbind(unit,x1,x2,x3,x4,x5,x6)
betahat=coefficients(data1.model)
yhat=x%*%betahat
ehat=residuals(data1.model)
ehat
r=rstudent(data1.model)
r

windows()
plot(yhat,r,ylab="R-Student residuals",xlab="Fitted values", main="Graph of R-student residuals against fitted values")
#windows()
#plot(y,r,main="graph of Life Expectancy against R-student residuals",ylab="R-student residuals",xlab="Life Expectancy in age")
windows()
plot(x[,2],r,main="graph of R-student residuals against Number of Infant Deaths",ylab="R-student residuals",xlab="Number of Infant Deaths per 10000 population")
windows()
plot(x[,3],r,main="graph of R-student residuals against Body Mass Index(BMI)",ylab="R-student residuals",xlab="Average Body Mass Index(BMI) of the entire population")
windows()
plot(x[,4],r,main="graph of R-student residuals against Under 5 Deaths",ylab="R-student residuals",xlab="Number of Under 5 Deaths per 10000 population")
windows()
plot(x[,5],r,main="graph of R-student residuals against Percentage of Polio immunisation",ylab="R-student residuals",xlab="Percentage of Polio (Pol3) immunisation coverage among 1 year olds")
windows()
plot(x[,6],r,main="graph of R-student residuals against Percentage DTP3 immunisation",ylab="R-student residuals",xlab="Percentage of DTP3 immunisation coverage among 1 year olds")
windows()
plot(x[,7],r,main="graph of R-student residuals against Income composition of resources",ylab="R-student residuals",xlab="Income composition of resources ranging from 0 to 1")



# Plot the normal probability plot
sorted_r=sort(r)
prob=(seq(1:16)-0.5)/n
windows()
plot(sorted_r,prob,main="Normal Probablity Plot",xlab="R-Student residuals",ylab="Probability")


# Find the hat matrix and its diagonal elements

h=hatvalues(data1.model)
h



# Find the ouliers or the leverage points

outlier=h[c(h>(2*p/n))]
outlier



# Cook's Distance

cd=cooks.distance(data1.model)
cd
influence=cd[c(cd>(4/(n-p)))]
influence

windows()
with(data1,plot(LifeExpectancy,cooks.distance(data1.model), xlab="Life Expectancy", ylab="Cook's Distance", main="Graph of Cook's Distance against Life Expectancy"))
identify(data1$LifeExpectancy,cd)




# Compute DFBETAS

dfb=dfbetas(data1.model)
dfb
num=seq(1:n)
beta0=num[c(abs(dfb[,1])>2/sqrt(n))]
beta1=num[c(abs(dfb[,2])>2/sqrt(n))]
beta2=num[c(abs(dfb[,3])>2/sqrt(n))]
beta3=num[c(abs(dfb[,4])>2/sqrt(n))]
beta4=num[c(abs(dfb[,5])>2/sqrt(n))]
beta5=num[c(abs(dfb[,6])>2/sqrt(n))]
beta6=num[c(abs(dfb[,7])>2/sqrt(n))]
beta0
beta1
beta2
beta3
beta4
beta5
beta6




# Compute DFFITS

dff=dffits(data1.model)
dff
outfits=dff[c(abs(dff)>2*sqrt(p/n))]
outfits



# Compute COVRATIO

covr=covratio(data1.model)
covr
outcov=covr[c(abs(covr-1)>3*p/n)]
outcov

improve=covr[c(covr>1)]
degrade=covr[c(covr<1)]
improve
degrade



# Find the mean squared residuals
msres=t(ehat)%*%ehat/(n-p)
msres=msres[1,1]
msres



# Estimate the true value of standard deviation of pure error

D=c(0)
E=c(0)
k=1
for(i in 1:(4*n-p))
{
	if(i<=(n-1))
	{
		if(i==1)
		{
			j=1
		}
		temp=betahat*(x[j,]-x[j+1,])
		temp1=(t(temp)%*%temp)/msres
		D[k]=temp1[1]
		E[k]=abs(ehat[j]-ehat[j+1])
		k=k+1
		j=j+1

	}
	if(i>(n-1) && i<=(n-1+n-2))
	{
		if(i==16)
		{
			j=1
		}
		temp=betahat*(x[j,]-x[j+2,])
		temp1=(t(temp)%*%temp)/msres
		D[k]=temp1[1]
		E[k]=abs(ehat[j]-ehat[j+2])
		k=k+1
		j=j+1

	}
	if(i>(n-1+n-2) && i<=(n-1+n-2+n-3))
	{
		if(i==30)
		{
			j=1
		}

		temp=betahat*(x[j,]-x[j+3,])
		temp1=(t(temp)%*%temp)/msres
		D[k]=temp1[1]
		E[k]=abs(ehat[j]-ehat[j+3])
		k=k+1
		j=j+1

	}
	if(i>(n-1+n-2+n-3) && i<=(n-1+n-2+n-3+n-4))
	{
		if(i==43)
		{
			j=1
		}

		temp=betahat*(x[j,]-x[j+4,])
		temp1=(t(temp)%*%temp)/msres
		D[k]=temp1[1]
		E[k]=abs(ehat[j]-ehat[j+4])
		k=k+1
		j=j+1

	}
}
mat=cbind(D,E)
matrix=mat[order(mat[,1]),]
Eu=matrix[1:15,2]
sigmahat=0.886*sum(Eu)/15
sigmahat
var=sigmahat*sigmahat
var



# Eliminating outliers and checking model

library(dplyr)

remove1=data1[-1,]
remove1.model=lm(LifeExpectancy~InfantDeaths+BMI+UnderFiveDeaths+Polio+Diphtheria+IncomeCompositionOfResources,data=remove1)
remove1.model
ehat1=residuals(remove1.model)
msres1=t(ehat1)%*%ehat1/(n-1-p)
msres1=msres1[1,1]
msres1
correlation1=cor(remove1[,3],remove1[4:9])


remove12=data1[-12,]
remove12.model=lm(LifeExpectancy~InfantDeaths+BMI+UnderFiveDeaths+Polio+Diphtheria+IncomeCompositionOfResources,data=remove12)
remove12.model
ehat12=residuals(remove12.model)
msres12=t(ehat12)%*%ehat12/(n-1-p)

msres12=msres12[1,1]
msres12
correlation12=cor(remove12[,3],remove12[4:9])


remove16=data1[-16,]
remove16.model=lm(LifeExpectancy~InfantDeaths+BMI+UnderFiveDeaths+Polio+Diphtheria+IncomeCompositionOfResources,data=remove16)
remove16.model
ehat16=residuals(remove16.model)
msres16=t(ehat16)%*%ehat16/(n-1-p)
msres16=msres16[1,1]
msres16
correlation16=cor(remove16[,3],remove16[4:9])


remove1_12=data1[-c(1,12),]
remove1_12.model=lm(LifeExpectancy~InfantDeaths+BMI+UnderFiveDeaths+Polio+Diphtheria+IncomeCompositionOfResources,data=remove1_12)
remove1_12.model
ehat1_12=residuals(remove1_12.model)
msres1_12=t(ehat1_12)%*%ehat1_12/(n-2-p)
msres1_12=msres1_12[1,1]
msres1_12
correlation1_12=cor(remove1_12[,3],remove1_12[4:9])


remove1_16=data1[-c(1,16),]
remove1_16.model=lm(LifeExpectancy~InfantDeaths+BMI+UnderFiveDeaths+Polio+Diphtheria+IncomeCompositionOfResources,data=remove1_16)
remove1_16.model
ehat1_16=residuals(remove1_16.model)
msres1_16=t(ehat1_16)%*%ehat1_16/(n-2-p)
msres1_16=msres1_16[1,1]
msres1_16
correlation1_16=cor(remove1_16[,3],remove1_16[4:9])


remove12_16=data1[-c(12,16),]
remove12_16.model=lm(LifeExpectancy~InfantDeaths+BMI+UnderFiveDeaths+Polio+Diphtheria+IncomeCompositionOfResources,data=remove12_16)
remove12_16.model
ehat12_16=residuals(remove12_16.model)
msres12_16=t(ehat12_16)%*%ehat12_16/(n-2-p)
msres12_16=msres12_16[1,1]
msres12_16
correlation12_16=cor(remove12_16[,3],remove12_16[4:9])


removeAll=data1[-c(1,12,16),]
removeAll.model=lm(LifeExpectancy~InfantDeaths+BMI+UnderFiveDeaths+Polio+Diphtheria+IncomeCompositionOfResources,data=removeAll)
removeAll.model
ehatAll=residuals(removeAll.model)
msresAll=t(ehatAll)%*%ehatAll/(n-3-p)
msresAll=msresAll[1,1]
msresAll
correlationAll=cor(removeAll[,3],removeAll[4:9])

betahat=coefficients(data1.model)
betahat1=coefficients(remove1.model)
betahat12=coefficients(remove12.model)
betahat16=coefficients(remove16.model)
betahat1_12=coefficients(remove1_12.model)
betahat1_16=coefficients(remove1_16.model)
betahat12_16=coefficients(remove12_16.model)
betahatAll=coefficients(removeAll.model)

data.frame(betahat,betahat1,betahat12,betahat16,betahat1_12,betahat1_16,betahat12_16,betahatAll)
data.frame(msres,msres1,msres12,msres16,msres1_12,msres1_16,msres12_16,msresAll,var)


diff1=abs((betahat-betahat1)/betahat)*100
diff12=abs((betahat-betahat12)/betahat)*100
diff16=abs((betahat-betahat16)/betahat)*100
diff1_12=abs((betahat-betahat1_12)/betahat)*100
diff1_16=abs((betahat-betahat1_16)/betahat)*100
diff12_16=abs((betahat-betahat12_16)/betahat)*100
diffAll=abs((betahat-betahatAll)/betahat)*100
data.frame(diff1,diff12,diff16,diff1_12,diff1_16,diff12_16,diffAll)






