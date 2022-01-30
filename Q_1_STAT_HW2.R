mu=1
a=0

beta=.5
M=3.1

variance=3
sigma=sqrt(variance)
area=1-pnorm(a,mu,sigma)
C=1/area
K=C/sqrt(2*pi*variance)

count=0
kk=1
f=matrix(NA,1000)
while (kk<1001)
{
z=rexp(1,beta)
fz=K*exp(-(z-mu)^2/(2*variance))
Mgz=M*beta*exp(-beta*z)
r=fz/Mgz
if(r>runif(1,min=0,max=1))
{f[kk]=z; kk=kk+1}
count=count+1
}
hist(f)
mean(f)
var(f)
