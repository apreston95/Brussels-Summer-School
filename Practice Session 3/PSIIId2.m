% This code writes down a model that nests complementarity, 
% habit and durables) and computes transitional dynamics
% Only in the case of Limit Cycles
% Code that will be solved using CSD version 2.4.0 (developped  by Sijmen Duineveld)

addpath ../../CoRRAM-M
addpath ../../CSD_v02.4.0
addpath ../../CSD_v02.4.0/subfun

clear all;
close all;

%% parameters

% ss hours is 1
par.chi     =   1;     % scale parameter for externality in durable sector
par.mu      =   1;
par.betaa    =   .98;    % discount factor
par.alphaa   =   .7;    % Returns to scale to labor in durable sector
% par.omega   =   .5;    % Disutility of labor
par.omega   =   1;    % Disutility of labor
par.gammaa   =   1;
par.stilde  =   1;

% par.nu      =   1.885;
par.nu      =   1.886345;


par.deltad   =   .025;   % depreciation
par.deltah   =   .5;   % depreciation
par.gammad  =   .4;
par.gammah  =   .4;
par.psii     =   .5;

par.sigmaa   =   4;



par.rho     =   .9;
par.sig     =   .01;




% compute SS
ellss=1;
css=par.chi;
dss=par.gammad/par.deltad*css;
hss=par.gammah/par.deltah*css;
sss=par.psii*dss-(1-par.psii)*par.gammaa*hss + css+par.stilde;
lambdadss=1/(1-par.betaa*(1-par.deltad)) *par.betaa*par.psii/sss;
lambdahss=1/(1-par.betaa*(1-par.deltah)) *(- par.betaa*(1-par.psii)*par.gammaa/sss);
par.omega0= par.alphaa*par.chi*(1/sss+par.gammad*lambdadss+par.gammah*lambdahss);

disp('     dss       hss       ellss     sss    lambdadss  lambdahss');
disp([dss hss ellss  sss lambdadss lambdahss ]);
disp('     css     omega0');
disp([css par.omega0]);


% f1 = 1/sss * par.alphaa*css/ellss - par.omega0*ellss^par.omega +par.alphaa*css/ellss*(par.gammad*lambdadss+par.gammah*lambdahss)
% f2 = -lambdadss+par.betaa*(par.psii/sss +(1-par.deltad)*lambdadss)
% f3 = -lambdahss+par.betaa*(-(1-par.psii)*par.gammaa/sss +(1-par.deltah)*lambdahss)
% f4 = hss-(1-par.deltah)*hss-par.gammah*css
% f5 = dss-(1-par.deltad)*dss-par.gammad*css
% f6 = -sss+par.psii*dss+css-(1-par.psii)*par.gammaa*hss+par.stilde
% f7 = css -par.chi*ellss^(par.alphaa*(1+par.nu*exp(-1/2*((ellss-1)/par.sigmaa)^2)))




% exponent on ell in the production function  
xell=[.5:.0001:1.5];
% xexpo=par.alphaa*(1+par.nu*exp(-1/2*((xell-ellss)/par.sigmaa).^2));
xexpo=(par.nu*exp(-1/2*((xell-ellss)/par.sigmaa).^2));

figure
plot(xell,xexpo,'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\nu$ in $\ell^{(1+\nu)\alpha}$','interpreter','latex','fontsize',20,'FontName','Times');







%% Model
[MOD] = EmergenceTwoStocks_CSD();

% Parameters:
MOD.par_val = [par.betaa,par.gammaa,par.omega,par.alphaa ...
               par.nu,par.sigmaa,par.chi,par.stilde,par.rho,par.sig,par.omega0, ...
               par.deltad,par.deltah,par.gammad,par.gammah,par.psii];
 
% Steady state:

MOD.SS_vec = [dss,hss,0,css,ellss,lambdadss,lambdahss,sss];

par.opt.order = 3;



% [SOL,NUM]  = pert_ana_csd(MOD,par.rho,par.opt.order,par.sig);
[SOL,NUM]  = pert_ana_csd_lim(MOD,par.rho,par.opt.order,par.sig);


%% Simulate transitionnal dynamics


ini_T   = 1;
TT      = 100000;

d           = NaN(ini_T+TT+1,1);
h           = NaN(ini_T+TT+1,1);
chishock    = NaN(ini_T+TT+1,1);
ell         = NaN(ini_T+TT,1);
lambdad     = NaN(ini_T+TT,1);
lambdah     = NaN(ini_T+TT,1);
s           = NaN(ini_T+TT,1);


dss          = 16.306832564816908;
hss          = 0.815341628230099;
ellss        = 1.013698082770772;

% Starting values:
% d(1,1)          = dss*1.0001;
d(1,1)          = 16.2;
h(1,1)          = hss;
% d(1,1)          = 16;
% h(1,1)          = .55;




chishock(1,1)   = 0;


for it = ini_T:ini_T+TT
    [stateV,commandV,chishock(it+1,:)] = eval_sol_csd(SOL,[d(it,:);h(it,:)],chishock(it,:),0,par.opt.order);
    d(it+1,:)       = stateV(1);
    h(it+1,:)       = stateV(2);
    c(it,:)         = commandV(1);
    ell(it,:)       = commandV(2);
    lambdad(it,:)   = commandV(3);
    lambdah(it,:)   = commandV(4);
    s(it,:)         = commandV(5);
end


%Remove last row for LK and LA:
d           = d(1:end-1,:);
h           = h(1:end-1,:);
chishock    = chishock(1:end-1,:);


%% PLOTS

figure
plot(h(ini_T:TT),d(ini_T:TT),'-','color',rgb('lightgrey'),'linewidth',2);
hold on
plot(hss,16.2,'o', ...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('dimgray'),...
    'MarkerSize',10);
plot(hss,dss,'o', ...
    'linewidth',1,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',10);
set(gca,'FontName','Times','FontSize',16);
xlabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$d$','interpreter','latex','fontsize',20,'FontName','Times');
hold off

figure('Name','first 100')
plot(h(ini_T:100),d(ini_T:100),'-','color',rgb('darkgrey'),'linewidth',2);
hold on
plot(hss,dss,'o', ...
    'linewidth',1,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',10);
set(gca,'FontName','Times','FontSize',16);
xlabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$d$','interpreter','latex','fontsize',20,'FontName','Times');
hold off


figure('Name','Last 200')
plot(h(TT-200:TT),d(TT-200:TT),'-','color',rgb('lightgrey'),'linewidth',2);
hold on
plot(hss,dss,'o', ...
    'linewidth',1,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',10);
set(gca,'FontName','Times','FontSize',16);
xlabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$d$','interpreter','latex','fontsize',20,'FontName','Times');
hold off

figure
subplot(311),plot(d(TT-99:TT),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$d$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(312),plot(h(TT-99:TT),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(313),plot(ell(TT-99:TT),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');


figure
subplot(211),plot(d(ini_T:TT),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$d$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(212),plot(h(ini_T:TT),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');


% end






function [MOD] = EmergenceTwoStocks_CSD()

%% BLOCK 1: Define parameters and variables
%Parameters:
syms betaa gammaa omega alphaa nu sigmaa chi stilde rho sig omega0 deltad deltah gammad gammah psii;

%Variables:
syms d_t d_n h_t h_n chishock_t chishock_n c_t c_n ell_t ell_n  lambdad_t lambdad_n  lambdah_t lambdah_n s_t s_n;


%% BLOCK 2: MODEL EQUATIONS (excl. stochastic process)

f1 = 1/s_t * alphaa*c_t/ell_t - omega0*ell_t^omega +alphaa*c_t/ell_t*(gammad*lambdad_t+gammah*lambdah_t);

f2 = -lambdad_t+betaa*(psii/s_n +(1-deltad)*lambdad_n);

f3 = -lambdah_t+betaa*(-(1-psii)*gammaa/s_n +(1-deltah)*lambdah_n);

f4 = h_n-(1-deltah)*h_t-gammah*c_t;

f5 = d_n-(1-deltad)*d_t-gammad*c_t;

f6 = -s_t+psii*d_t+c_t-(1-psii)*gammaa*h_t+stilde;

f7 = c_t -exp(chishock_t)*chi*ell_t^(alphaa*(1+nu*exp(-1/2*((ell_t-1)/sigmaa)^2)));



%% BLOCK 3: ASSIGNMENTS
%Model equations:
MOD.FS = [f1;f2;f3;f4;f5;f6;f7];

%Endogenous state variables:
MOD.XX = [d_t,h_t];
%Exogenous state variables:
MOD.ZZ = [chishock_t];
%Control variables:
MOD.YY = [c_t, ell_t,lambdad_t,lambdah_t,s_t];

MOD.XXn = [d_n,h_n];
MOD.ZZn = [chishock_n];
MOD.YYn = [c_n, ell_n,lambdad_n,lambdah_n,s_n];

% Variable names (strings):
MOD.var_bs_nms  = {'d','h','chishock','c','ell','lambdad','lambdah','s'};
% Parameter names (strings):
MOD.par_nms     = {'betaa','gammaa','omega ','alphaa ',...
    'nu ','sigmaa ','chi ','stilde ','rho ','sig ','omega0', ...
    'deltad','deltah','gammad','gammah','psii'};



end


