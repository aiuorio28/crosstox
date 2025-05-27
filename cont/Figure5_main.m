clc
keep pphome 
%% Figure 6 (bifurcation diagram - parameter s)
% parameters kept fixed
g=10; % growth rate
d=1; % mortality rate                              ! g>d !
c=0.5; % conversion factor
k=1; % toxixity depletion
Rc=6; % sort of carrying capacity for R
dR=3/0.9; % diffusion parameter (biomass)
dT=0.05; % diffusion parameter (toxicity)

%% parameters varied
gamma=0.1; % growth inhibition
s=0.5; % extramortality (bifurcation parameter in Fig.6)
psigma=0.85;%bifurcation parameter in Fig.5

par=[g,gamma,Rc,d,s,c,k,dR,dT,psigma]';

lx=4;  % domain (interval)
p=[];
p=ToxicityX1Dinit(p,lx,400,par); 
p.Om=2*lx;
p=setfn(p,'1_stationary'); 
p.nc.ilam=10; % continuation wrt p_sigma
%% contiunation of the trivial branch
p.nc.dsmax=1e-3; 
p.sw.bifcheck=2;
p.sol.ds=0.01; % starting stepsize
p.nc.lammax=1; 
p.nc.lammin=0.8; 
p=cont(p,1500); % continuation of the homogeneous branch

% File list
currentFolder=pwd;
cd('1_stationary')
ptFileList=dir('bpt*.mat');
[nn,~]=size(ptFileList);
cd(currentFolder)

%% continuation from bif points
for i=1:nn
    BPT_list=['bpt' num2str(i)];
    Branch_i_u=['bpt' num2str(i) '_u'];
    Branch_i_d=['bpt' num2str(i) '_d'];
    p=swibra('1_stationary',BPT_list,Branch_i_d,1e-2);  p.nc.dsmax=1e-1; p.sw.bifcheck=1; p.plot.pcmp=1; p=cont(p,100);
    p=swibra('1_stationary',BPT_list,Branch_i_u,1e-2);  p.nc.dsmax=1e-2; p.sw.bifcheck=1; p.plot.pcmp=1; p=cont(p,500);
end

%% Postprocessing, plot BifDiagram 
nfig=5;
figure(nfig);
clf(nfig); 
cmp=2; 
box on
hold on
plotbra('1_stationary',nfig,cmp,'cl','k'); 

% in grey
for i=[1, 2:nn]
    Branch_i_u=['bpt' num2str(i) '_u'];
    Branch_i_d=['bpt' num2str(i) '_d'];
    plotbra(Branch_i_d,nfig,cmp,'cl',[0.55 0.57 0.67]);
    %plotbra(Branch_i_u,nfig,cmp,'cl',[0.55 0.57 0.67]);
end
plotbra('bpt2_d',nfig,cmp,'cl',[1,0.55,0]);
%plotbra('bpt2_u',nfig,cmp,'cl',[1,0.55,0]);
xlabel('\sigma/d_R')
ylabel('||R||_{L^1}')
axis([0 0.85 302 347])