n=50
a=0.001
b=0.001

s=10000
B=matrix(0,4,1)
S=s*diag(4)
m=10000
Xsource=matrix(NA,50,4)
y=matrix(NA,50,1)

Xsource=as.matrix(read.table("C:\\Users\\somayehyoussefi\\Desktop\\HW2.txt")[1:50,2:5])
y=as.matrix(read.table("C:\\Users\\somayehyoussefi\\Desktop\\HW2.txt")[1:50,1])
Xsource=t(Xsource)
y=t(y)
x=matrix(0,4,50)

Like=matrix(NA,15,1)
Likei=matrix(NA,m,15)

sigma2=matrix(1,m,15)
beta=array(0,dim=c(15,m,4))
M1=matrix(NA,4,1)
M2=matrix(NA,4,4)

for (count in 1:15)
{
if (count==1){
            x=Xsource
        
} else if (count==2){
		x=matrix(0,4,50)
            x[2,]=Xsource[2,];x[3,]=Xsource[3,];x[4,]=Xsource[4,]
           
} else if (count==3){
		x=matrix(0,4,50)
	      x[1,]=Xsource[1,];x[3,]=Xsource[3,];x[4,]=Xsource[4,]
         
} else if (count==4){
		x=matrix(0,4,50)
	    x[1,]=Xsource[1,];x[2,]=Xsource[2,];x[4,]=Xsource[4,]

} else if (count==5){
		x=matrix(0,4,50)
	    x[1,]=Xsource[1,];x[2,]=Xsource[2,];x[3,]=Xsource[3,]

} else if (count==6){
		x=matrix(0,4,50)
	    x[3,]=Xsource[3,];x[4,]=Xsource[4,]

} else if (count==7){
		x=matrix(0,4,50)
	    x[2,]=Xsource[2,];x[4,]=Xsource[4,]

} else if (count==8){
		x=matrix(0,4,50)
	   x[2,]=Xsource[2,];x[3,]=Xsource[3,]

} else if (count==9){
		x=matrix(0,4,50)
	    x[1,]=Xsource[1,];x[2,]=Xsource[2,]

} else if (count==10){
		x=matrix(0,4,50)
	   x[1,]=Xsource[1,];x[3,]=Xsource[3,]

} else if (count==11){
		x=matrix(0,4,50)
	   x[1,]=Xsource[1,];x[4,]=Xsource[4,]

} else if (count==12){
		x=matrix(0,4,50)
	   x[1,]=Xsource[1,]

} else if (count==13){
		x=matrix(0,4,50)
	   x[2,]=Xsource[2,]

} else if (count==14){
		x=matrix(0,4,50)
	    x[3,]=Xsource[3,]

} else {
		x=matrix(0,4,50)
            x[4,]=Xsource[4,]
	
}
sum=0
for(i in 2:m)
{
Error=matrix(NA,1,50)
invS=solve(S)
xyt=x%*%(t(y))
M2=solve(1/sigma2[i-1,count]*(x%*%t(x))+invS)
M1=M2%*%(1/sigma2[i-1,count]*xyt+invS%*%B)
rr=matrix(0,4,1)
rr[1,1]=rnorm(1,0,1)
rr[2,1]=rnorm(1,0,1)
rr[3,1]=rnorm(1,0,1)
rr[4,1]=rnorm(1,0,1)
beta[count,i,]=M1+chol(M2)%*%rr

terma=a+n/2
termb=1/2*((y-t(beta[count,i,])%*%x)%*%t(y-t(beta[count,i,])%*%x))+b
sigma2[i,count]=1/rgamma(1,terma,termb)

Error=y-t(beta[count,i,])%*%x
Error2=Error%*%t(Error)

Likei[i,count]=exp(-Error2/(2*sigma2[i,count]))/(sqrt(2*pi*sigma2[i,count]))^50
sum=sum+1/Likei[i,count]
}
Like[count,1]=m/sum

}

BayesFactor=matrix(NA,15,15)
for (k in 1:15)
{
for (l in 1:15)
{
BayesFactor[k,l]=Like[k,1]/Like[l,1]
}
}
write.table(BayesFactor,"Bayes.csv", sep=",")
write.table(Likei,"Likei.csv", sep=",")
MEAN=matrix(NA,15,4)
for (i in 1:15)
{
for (j in 1:4)
{
MEAN[i,j]=mean(beta[i,,j])
}
}
write.table(MEAN,"mmmm.csv", sep=",")


par(mfrow = c(4, 4))
hist(beta[1,,1])
hist(beta[1,,2])
hist(beta[1,,3])
hist(beta[1,,4])
hist(beta[2,,1])
hist(beta[2,,2])
hist(beta[2,,3])
hist(beta[2,,4])

hist(beta[3,,1])
hist(beta[3,,2])
hist(beta[3,,3])
hist(beta[3,,4])

hist(beta[4,,1])
hist(beta[4,,2])
hist(beta[4,,3])
hist(beta[4,,4])

par(mfrow = c(4, 4))
hist(beta[5,,1])
hist(beta[5,,2])
hist(beta[5,,3])
hist(beta[5,,4])
hist(beta[6,,1])
hist(beta[6,,2])
hist(beta[6,,3])
hist(beta[6,,4])

hist(beta[7,,1])
hist(beta[7,,2])
hist(beta[7,,3])
hist(beta[7,,4])

hist(beta[8,,1])
hist(beta[8,,2])
hist(beta[8,,3])
hist(beta[8,,4])

par(mfrow = c(4, 4))
hist(beta[9,,1])
hist(beta[9,,2])
hist(beta[9,,3])
hist(beta[9,,4])
hist(beta[10,,1])
hist(beta[10,,2])
hist(beta[10,,3])
hist(beta[10,,4])

hist(beta[11,,1])
hist(beta[11,,2])
hist(beta[11,,3])
hist(beta[11,,4])

hist(beta[12,,1])
hist(beta[12,,2])
hist(beta[12,,3])
hist(beta[12,,4])

par(mfrow = c(3, 4))
hist(beta[13,,1])
hist(beta[13,,2])
hist(beta[13,,3])
hist(beta[13,,4])
hist(beta[14,,1])
hist(beta[14,,2])
hist(beta[14,,3])
hist(beta[14,,4])

hist(beta[15,,1])
hist(beta[15,,2])
hist(beta[15,,3])
hist(beta[15,,4])


C=matrix(0,1,1)
index=matrix(1,m,1)
C<-Likei[2,1]
n=m-1
for (i in 2:n)
{
for (j in 1:15)
{
compare<-Likei[i,j]
if (compare>=C) 
{
C<-compare
index[i,1]<-j
}
mm<-index[i,1]
}
C<-Likei[i+1,mm] 
}


sum=matrix(0,15,1)

for (i in 1:m)
{
for (j in 1:15)
{
if (index[i,1]==j)
{ sum[j,1]=sum[j,1]+1
}
}
}

write.table(sum,"number.csv", sep=",")
