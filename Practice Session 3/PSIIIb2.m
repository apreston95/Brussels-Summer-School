% Code that will be solved using CSD version 2.4.0 (developped  by Sijmen Duineveld)

% addpath ../../CoRRAM-M
% addpath ../../CSD_v02.4.0
% addpath ../../CSD_v02.4.0/subfun

clear all;
close all;

%% parameters



% ss hours is 1
par.betaa=.98;
par.gammaa = -1/2;
par.deltaa=.025;
par.omega=.5;
par.alphaa=.7;
% par.nu=1.08;
par.nu=1.09;
par.sigmaa=.1;
par.chi=1;
par.stilde=max(-1.1*(par.deltaa-par.gammaa)/par.deltaa*par.chi,1);

par.rho     = 0.9;
par.sig   = .01;


% compute SS
ellss=1;
css=par.chi;
hss=css/par.deltaa;
lambdass = par.betaa*par.gammaa / ((par.deltaa-par.gammaa)*hss +par.stilde) ...
    *1/(par.betaa*(1-par.deltaa)-1); 

par.omega0=par.alphaa*css/(css-par.gammaa*hss+par.stilde) ...
    + par.alphaa*lambdass*css;


disp('     hss       ellss     css    lambdass    omega0');
disp([hss ellss  css lambdass  par.omega0]);


% SS check
% f1 = par.alphaa*css/(ellss*(css-par.gammaa*hss)+stilde) - par.omega0*ellss^par.omega +lambdass*par.alphaa*css/ellss 
% f2 = -lambdass+par.betaa*(-par.gammaa/(css - par.gammaa*hss+par.stilde) +(1-par.deltaa)*lambdass)
% f3 = hss-(1-par.deltaa)*hss-css
% f4 = css -exp(0)*par.chi*ellss^(par.alphaa*(1+par.nu*exp(-1/2*((ellss-1)/par.sigmaa)^2)))










%% Model
[MOD] = EmergenceDurable_CSD();

% Parameters:
MOD.par_val = [par.betaa,par.gammaa,par.deltaa,par.omega,par.alphaa ...
               par.nu,par.sigmaa,par.chi,par.stilde,par.rho,par.sig,par.omega0];
 
% Steady state:

MOD.SS_vec = [hss,0,css,lambdass,ellss];

par.opt.order = 3;



% [SOL,NUM]  = pert_ana_csd(MOD,par.rho,par.opt.order,par.sig);
[SOL,NUM]  = pert_ana_csd_lim(MOD,par.rho,par.opt.order,par.sig);


%% Simulate transitionnal dynamics


ini_T   = 1;
TT      = 2000;

h = NaN(ini_T+TT+1,1);
chishock = NaN(ini_T+TT+1,1);
c = NaN(ini_T+TT,1);
lambda = NaN(ini_T+TT,1);
ell = NaN(ini_T+TT,1);

%Starting values:

h(1,1)   = hss*1.001;
% hc1 =   39.3202;     
% hc2 =   40.6659;
% h(1,1)   = hc2*1.01;

chishock(1,1)     = 0;


for it = ini_T:ini_T+TT
    [h(it+1,:),commandV,chishock(it+1,:)] = eval_sol_csd(SOL,h(it,:),chishock(it,:),0,par.opt.order);
    c(it,:) = commandV(1);
    lambda(it,:) = commandV(3);
    ell(it,:) = commandV(3);
end


%Remove last row for LK and LA:
h = h(1:end-1,:);
chishock = chishock(1:end-1,:);

% The two elements of the 2-cycle
hc1=min(h(end-1),h(end));
hc2=max(h(end-1),h(end));


%% PLOTS

figure
plot(h(ini_T:300),ell(ini_T:300),'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',10);
hold on
plot(hss,ellss,'o',...
    'linewidth',1,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
plot(h(300),ell(300),'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('dimgray'),...
    'MarkerSize',12);
plot(h(301),ell(301),'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('dimgray'),...
    'MarkerSize',12);
set(gca,'FontName','Times','FontSize',16);
xlabel('$h_t$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell_t$','interpreter','latex','fontsize',20,'FontName','Times');
hold off


figure
% subplot(211),plot(ell(ini_T:TT),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
% set(gca,'FontName','Times','FontSize',16);
% xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
% ylabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');

plot(ell(ini_T:30),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');

print('ell_sim','-dpng');


%% 


figure
plot(h(ini_T:TT-1),h(ini_T+1:TT),'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',10);
hold on
for t=ini_T:TT-1
plot([h(t) h(t+1)],[h(t+1) h(t+1)],'--','linewidth',1,'color',rgb('lightgray'));
plot([h(t) h(t)],[h(t) h(t+1)],'--','linewidth',1,'color',rgb('lightgray'));
end
plot(h(ini_T:TT-1),h(ini_T+1:TT),'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',10);
plot(hss,hss,'o',...
    'linewidth',1,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
plot(h(300),h(301),'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('dimgray'),...
    'MarkerSize',12);
plot(h(301),h(302),'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('dimgray'),...
    'MarkerSize',12);
set(gca,'FontName','Times','FontSize',16);
xlabel('$h_{t}$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h_{t+1}$','interpreter','latex','fontsize',20,'FontName','Times');
aa=axis;
plot(aa(1:2),aa(3:4),'-','linewidth',1,'color',rgb('black'));
axis(aa);
hold off




figure
plot(ell(ini_T:30),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
title('(b2) $\nu$ at the bifurcation','interpreter','latex','fontsize',20,'FontName','Times');
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell_t$','interpreter','latex','fontsize',20,'FontName','Times');
pbaspect([3 1 1]);
print -depsc2 EmLCSimplef2


% return

%% Trace out the (h_t,h_{t+1}) locus  

ini_T   = 1;
TT      = 50;

h = NaN(ini_T+TT+1,1);
chishock = NaN(ini_T+TT+1,1);
c = NaN(ini_T+TT,1);
ell = NaN(ini_T+TT,1);
chishock(1,1)     = 0;

htot=[];
htotp1=[];

%Starting values:

for hstart=linspace(hc1*.995,hc2*1.005,1000)
% for hstart=linspace(hss*.99,hss*1.01,100)

    h(1,1)   = hstart;

    for it = ini_T:ini_T+TT
        [h(it+1,:),commandV,chishock(it+1,:)] = eval_sol_csd(SOL,h(it,:),chishock(it,:),0,par.opt.order);
        c(it,:) = commandV(1);
        lambda(it,:) = commandV(3);
        ell(it,:) = commandV(3);
    end

    %Remove last row for LK and LA:
    h = h(1:end-1,:);

    htot=[htot;h(ini_T:TT-1)];
    htotp1=[htotp1;h(ini_T+1:TT)];

end



figure
plot(htot,htotp1,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerEdgeColor',rgb('lightgray'), 'MarkerSize',4);
hold on
plot(hss,hss,'o',...
    'linewidth',1,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
plot(hc1,hc2,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',12);
plot(hc2,hc1,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',12);
set(gca,'FontName','Times','FontSize',16);
xlabel('$h_{t}$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h_{t+1}$','interpreter','latex','fontsize',20,'FontName','Times');
aa=[min(htot) max(htot) min(htot) max(htot)];
plot(aa(1:2),aa(1:2),'--','linewidth',1,'color',rgb('black'));
axis([aa(1:2) aa(1:2)]);
hold off
print('Phase_Plot_Durables','-dpng');


%%
function [MOD] = EmergenceDurable_CSD()

%% BLOCK 1: Define parameters and variables
%Parameters:
syms betaa gammaa deltaa omega alphaa nu sigmaa chi stilde rho sig omega0;
%Variables:
syms c_t  c_n h_t h_n lambda_t lambda_n ell_t ell_n chishock_t chishock_n;
%% BLOCK 2: MODEL EQUATIONS (excl. stochastic process)

f1 = alphaa*c_t/ell_t* 1/((c_t-gammaa*h_t)+stilde) - omega0*ell_t^omega +lambda_t*alphaa*c_t/ell_t;

f2 = -lambda_t+betaa*(-gammaa/(c_n - gammaa*h_n+stilde) +(1-deltaa)*lambda_n);

f3 = h_n-(1-deltaa)*h_t-c_t;

f4 = c_t -exp(chishock_t)*chi*ell_t^(alphaa*(1+nu*exp(-1/2*((ell_t-1)/sigmaa)^2)));
%% BLOCK 3: ASSIGNMENTS
%Model equations:
MOD.FS = [f1;f2;f3;f4];
%Endogenous state variables:
MOD.XX = [h_t];
%Exogenous state variables:
MOD.ZZ = [chishock_t];
%Control variables:
MOD.YY = [c_t, lambda_t, ell_t];

MOD.XXn = [h_n];
MOD.ZZn = [chishock_n];
MOD.YYn = [c_n, lambda_n, ell_n];

% Variable names (strings):
MOD.var_bs_nms  = {'h','chishock','c','lambda','ell'};
% Parameter names (strings):
MOD.par_nms     = {'betaa','gammaa','deltaa','omega ','alphaa ',...
    'nu ','sigmaa ','chi ','stilde ','rho ','sig ','omega0'};

end


