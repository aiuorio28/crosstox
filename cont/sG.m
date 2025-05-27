function r=sG(p,u) 
% compute pde-part of residual
R=u(1:p.np); % extract the first component
T=u(p.np+1:2*p.np); % extract the second component
par=u(p.nu+1:end); % extract parameters
% par=[g,gamma,Rc,d,s,c,k,dR,dT,sigma]';
g=par(1); gamma=par(2); Rc=par(3); d=par(4); s=par(5); c=par(6); k=par(7);
dR=par(8); dT=par(9); psigma=par(10);

Tc=Rc*c*(d+s)/k;
sigma=psigma*dR;
tT=min(T/Tc,1);

f1=(g-gamma*tT).*R.*(1-R/Rc)-(d+s*tT).*R;
f2=c*(d+s*tT).*R-k*T;

%-----------------
K=p.mat.K;
r=[dR*K*R-sigma*K*(tT.*R);...
   dT*K*T]-p.mat.M*[f1;f2];
%-----------------
end