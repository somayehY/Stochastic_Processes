n=50
a=0.001
b=0.001

s=10000
B=matrix(0,4,1)
S=s*diag(4)
m=5000
x=matrix(NA,50,4)
y=matrix(NA,50,1)
x=as.matrix(read.table("C:\\Users\\somayehyoussefi\\Desktop\\HW2.txt")[1:50,2:5])
y=as.matrix(read.table("C:\\Users\\somayehyoussefi\\Desktop\\HW2.txt")[1:50,1])
x=t(x)
y=t(y)

sigma2=matrix(1,m)
beta=matrix(0,m,4)
M1=matrix(NA,4,1)
M2=matrix(NA,4,4)

for(i in 2:m)
{

invS=solve(S)
xyt=x%*%(t(y))
M2=solve(1/sigma2[i-1]*(x%*%t(x))+invS)
M1=M2%*%(1/sigma2[i-1]*xyt+invS%*%B)
rr=matrix(0,4,1)
rr[1,1]=rnorm(1,0,1)
rr[2,1]=rnorm(1,0,1)
rr[3,1]=rnorm(1,0,1)
rr[4,1]=rnorm(1,0,1)
beta[i,]=M1+chol(M2)%*%rr
terma=a+n/2
termb=1/2*((y-t(beta[i,])%*%x)%*%t(y-t(beta[i,])%*%x))+b
sigma2[i]=1/rgamma(1,terma,termb)
}
hist(beta[,1])
hist(beta[,2])
hist(beta[,3])
hist(beta[,4])
hist(sigma2)



