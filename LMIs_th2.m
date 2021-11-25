function [OmegaVal,OmegaVal0,Omegaval1,Omegaval2,Omegaval3,Omegaval4,Omegaval5,Omegaval6,Omegaval7,Omegaval8,Omegaval9,Omegaval10]=LMIs_th2(M,dL,dR,delta,rs,a1,a2,deltau,delta0,gamma,mu,E)

%参数
edelta=exp(-2*delta*rs);

%variables
lambda=sdpvar(1);
lambda1=sdpvar(1);
P1=sdpvar(1,1);
P2=sdpvar(1,1);
P3=sdpvar(1,1);
P4=sdpvar(1,1);
P5=sdpvar(1,1);
P=sdpvar(1,1);
R1=sdpvar(1,1);
R2=sdpvar(1,1);
W1=sdpvar(1,1);
W2=sdpvar(1,1);


%constraints
Phi=blkvar;%分块矩阵变量类型
Phi(1,1)=-(max([dL,dR])*pi^2/(4-3*dL*dR))*P-(edelta/rs)*P4+2*delta*P1+P3+2*E*P2;
Phi(1,2)=P1+P2*E-P2;
Phi(1,3)=(edelta/rs)*P4-W1;
Phi(1,4)=W1;
Phi(2,2)=-2*P2+rs*P4;
Phi(2,3)=-W1;
Phi(2,4)=W1;
Phi(3,3)=-edelta*P3-edelta*P4/rs;
Phi(4,4)=-pi^2*lambda/deltau^2;
Phi=sdpvar(Phi);

Phi1=blkvar;%分块矩阵变量类型
Phi1(1,1)=-(max([dL,dR])*pi^2/(4-3*dL*dR))*P-(edelta/rs)*P4+2*delta*P1+P3+2*E*P2;
Phi1(1,2)=P1+P2*E-P2;
Phi1(1,3)=(edelta/rs)*P4-W2;
Phi1(1,4)=W2;
Phi1(2,2)=-2*P2+rs*P4;
Phi1(2,3)=-W2;
Phi1(2,4)=W2;
Phi1(3,3)=-edelta*P3-edelta*P4/rs;
Phi1(4,4)=-pi^2*lambda/deltau^2;
Phi1=sdpvar(Phi1);

Psi=blkvar;
Psi(1,1)=-2*gamma*(1-delta)*P2+P5+max([dL,dR])*P;
Psi=sdpvar(Psi);

Gamma=blkvar;
Gamma(1,1)=-edelta*P5+lambda;
Gamma=sdpvar(Gamma);

Psi2=blkvar;
Psi2(1,1)=-2*gamma*P2+lambda1;
Psi2=sdpvar(Psi2);

Theata=blkvar;
Theata(1,1)=2*a1*P2-R1'-R1+2*delta0*P2+.5*mu*(W1'+W1);
Theata(1,2)=-R1;
Theata(2,2)=-(pi^2/deltau^2);
Theata=sdpvar(Theata);

Theata1=blkvar;
Theata1(1,1)=2*a2*P2-R2'-R2+2*delta0*P2+.5*mu*(W2'+W2);
Theata1(1,2)=-R2;
Theata1(2,2)=-(pi^2/deltau^2);
Theata1=sdpvar(Theata1);

%solution
LMIs=[lambda>0,lambda1>0,P1>=0, P2>=0,P3>=0, P4>=0, P5>=0, P>=0, W1>0,W2>0,R1>0,R2>0,Phi<=0,Phi1<=0, Psi<=0, Gamma<=0, Psi2<=0, Theata<=0,Theata1<=0]; 
options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,[],options); 

OmegaVal=[]; 
OmegaVal0=[];
Omegaval1=[]; 
Omegaval2=[]; 
Omegaval3=[];
Omegaval4=[]; 
Omegaval5=[];
Omegaval6=[];
Omegaval7=[];
Omegaval8=[];
Omegaval9=[];
Omegaval10=[];


if sol.problem == 0
    primal=check(LMIs); 
if min(primal)>=0 && all(primal([1,2,3,4,5,6,7,8,9,10,11,12])>0)
       OmegaVal=double(lambda);
       OmegaVal0=double(lambda1);
         Omegaval1=double(P1);
       Omegaval2=double(P2);
       Omegaval3=double(P3);
        Omegaval4=double(P4);
        Omegaval5=double(P5);
        Omegaval6=double(P);
        Omegaval7=double(W1);
        Omegaval8=double(W2);
        Omegaval9=double(R1);
        Omegaval10=double(R2);
end
else
    yalmiperror(sol.problem); 
end 