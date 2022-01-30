n=50
a=0.001
b=0.001

s=10000
B=matrix(0,4,1)
S=s*diag(4)
m=5000

y=matrix(NA,50,1)
x=as.matrix(read.table("C:\\Users\\somayehyoussefi\\Desktop\\HW2.txt")[1:50,2:5])
y=as.matrix(read.table("C:\\Users\\somayehyoussefi\\Desktop\\HW2.txt")[1:50,1])
x=t(x)
y=t(y)

sigma2=matrix(1,m)
beta=matrix(0,m,4)
M1=matrix(NA,4,1)
M2=matrix(NA,4,4)
rr=matrix(0,4,1)
Like=matrix(NA,15,1)
Likei=matrix(NA,m,1)
count=1

	

for(i in 2:m)
sum=0
{

for (count in 1:15)
{
if (count==1){
            rr[1,1]=rnorm(1,0,1); rr[2,1]=rnorm(1,0,1); rr[3,1]=rnorm(1,0,1); rr[4,1]=rnorm(1,0,1);
} else if (count==2){
            rr[1,1]=0; rr[2,1]=rnorm(1,0,1); rr[3,1]=rnorm(1,0,1); rr[4,1]=rnorm(1,0,1);S
} else if (count==3){
		rr[1,1]=rnorm(1,0,1); rr[2,1]=0; rr[3,1]=rnorm(1,0,1); rr[4,1]=rnorm(1,0,1);
} else if (count==4){
		rr[1,1]=rnorm(1,0,1); rr[2,1]=rnorm(1,0,1); rr[3,1]=0; rr[4,1]=rnorm(1,0,1);
} else if (count==5){
		rr[1,1]=rnorm(1,0,1); rr[2,1]=rnorm(1,0,1); rr[3,1]=rnorm(1,0,1); rr[4,1]=0;
} else if (count==6){
		rr[1,1]=0; rr[2,1]=0; rr[3,1]=rnorm(1,0,1); rr[4,1]=rnorm(1,0,1);
} else if (count==7){
		rr[1,1]=0; rr[2,1]=rnorm(1,0,1); rr[3,1]=0; rr[4,1]=rnorm(1,0,1);
} else if (count==8){
		rr[1,1]=0; rr[2,1]=rnorm(1,0,1); rr[3,1]=rnorm(1,0,1); rr[4,1]=0;
} else if (count==9){
		rr[1,1]=rnorm(1,0,1); rr[2,1]=rnorm(1,0,1); rr[3,1]=0; rr[4,1]=0;
} else if (count==10){
		rr[1,1]=rnorm(1,0,1); rr[2,1]=0; rr[3,1]=rnorm(1,0,1); rr[4,1]=0;
} else if (count==11){
		rr[1,1]=rnorm(1,0,1); rr[2,1]=0; rr[3,1]=0; rr[4,1]=rnorm(1,0,1);
} else if (count==12){
		rr[1,1]=rnorm(1,0,1); rr[2,1]=0; rr[3,1]=0; rr[4,1]=0;
} else if (count==13){
		rr[1,1]=0; rr[2,1]=rnorm(1,0,1); rr[3,1]=0; rr[4,1]=0;
} else if (count==14){
		rr[1,1]=0; rr[2,1]=0; rr[3,1]=rnorm(1,0,1); rr[4,1]=0;
} else {
		rr[1,1]=0; rr[2,1]=0; rr[3,1]=0; rr[4,1]=rnorm(1,0,1);
}


Error=matrix(NA,1,50)

invS=solve(S)
xyt=x%*%(t(y))
M2=solve(1/sigma2[i-1]*(x%*%t(x))+invS)
M1=M2%*%(1/sigma2[i-1]*xyt+invS%*%B)
beta[i,]=M1+chol(M2)%*%rr
Error=y-t(beta[i,])%*%x
Error2=Error%*%t(Error)

terma=a+n/2
termb=1/2*((y-t(beta[i,])%*%x)%*%t(y-t(beta[i,])%*%x))+b
sigma2[i]=1/rgamma(1,terma,termb)

Likei[i,1]=exp(-Error2/(2*sigma2[i]))/(sqrt(2*pi*sigma2[i]))^50
sum=sum+1/Likei[i,1]


Like[count,1]=m/sum
}



hist(beta[,1])
hist(beta[,2])
hist(beta[,3])
hist(beta[,4])
hist(sigma2)



