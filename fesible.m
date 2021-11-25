clc;
clear;

dL=1;
dR=0;
non=0.8;
gamma=0.2;
r=1.2;
M=3;
rs=r/M;
v=0.25;
a1=non*(1-v);
a2=non;
delta=0.2;
delta0=0.15;
mu=0.5;
deltau=0.5;
E=abs(a1-a2);


[OmegaVal,OmegaVal0,Omegaval1,Omegaval2,Omegaval3,Omegaval4,Omegaval5,Omegaval6,Omegaval7,Omegaval8,Omegaval9,Omegaval10]=LMIs_th2(M,dL,dR,delta,rs,a1,a2,deltau,delta0,gamma,mu,E)
