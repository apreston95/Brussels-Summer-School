% Model with habit,
% produces hysteresis
% for simulations to low and high SS, 
% use with dynare



close all;


var h ell c lambda F;    

parameters ellss00 beta gamma omega0 omega delta alpha chi nu sigma stilde gammah;


ellss00=1;
beta=.98;
gammah = 1/2;
omega=.5;
delta=.9;
gamma=1/2;
alpha=.7;
% nu=0;
nu=3;
% nu=2.6906;
sigma=.1;
chi=1;
stilde= max(-1.1*(delta-gamma*gammah)/delta*chi,1);






% exponent on ell in the production function  
xell=[.5:.01:1.5];
xexpo=alpha*(1+nu*exp(-1/2*((xell-ellss00)/sigma).^2));

figure
plot(xell,xexpo,'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('Exponent','interpreter','latex','fontsize',20,'FontName','Times');
print('Exponent_Plot','-dpng');



% SS

Fs=ellss00^(alpha*nu);
cs=chi*Fs*ellss00^(alpha);
hs=gammah*cs/delta;
lambdas = beta*gamma / ((delta/gammah-gamma)*hs +stilde) *1/(beta*(1-delta)-1); 

omega0=alpha*cs/(cs-gamma*hs+stilde) + alpha*gammah*lambdas*cs  ;



disp('     hs        ells    lambdas     cs      omega0');
disp([hs ellss00 lambdas cs omega0]);


%% Multiple SS

sell=[.94:.00001:1.045];
snu=nu;

inu=0;
for ssnu=snu;
    inu=inu+1;
    ill=0;

    for ssell=sell;
        ill=ill+1;
        sc          =   chi*ssell^(alpha*(1+ssnu*exp(-1/2*((ssell-ellss00)/sigma)^2)));
        sh          =   gammah*sc/delta;
        slambda     =   beta*gamma / ((delta/gammah-gamma)*sh +stilde) ...
            *1/(beta*(1-delta)-1);
        Well(ill,inu)   =   alpha*sc/ssell*1/(sc-gamma*sh+stilde) - omega0*ssell^omega +slambda*gammah*alpha*sc/ssell;
    end

end

[ssval ssindex]=sort(abs(Well));

% sell(ssindex(1:10))

sstates=sort([sell(ssindex(1)) sell(ssindex(2)) sell(ssindex(3))]);

ellss1= sstates(1);
ellss0= sstates(2);
ellss2= sstates(3);



css0          =   chi*ellss0^(alpha*(1+nu*exp(-1/2*((ellss0-ellss00)/sigma)^2)));
hss0          =   gammah*css0/delta;
lambdass0     =   beta*gamma / ((delta/gammah-gamma)*hss0 +stilde) ...
    *1/(beta*(1-delta)-1);
Fss0             = ellss0^(alpha*nu*exp(-1/2*((ellss0-ellss00)/sigma)^2));

css1          =   chi*ellss1^(alpha*(1+nu*exp(-1/2*((ellss1-ellss00)/sigma)^2)));
hss1          =   gammah*css1/delta;
lambdass1     =   beta*gamma / ((delta/gammah-gamma)*hss1 +stilde) ...
    *1/(beta*(1-delta)-1);
Fss1             = ellss1^(alpha*nu*exp(-1/2*((ellss1-ellss00)/sigma)^2));

css2          =   chi*ellss2^(alpha*(1+nu*exp(-1/2*((ellss2-ellss00)/sigma)^2)));
hss2          =   gammah*css2/delta;
lambdass2     =   beta*gamma / ((delta/gammah-gamma)*hss2 +stilde) ...
    *1/(beta*(1-delta)-1);
Fss2             = ellss2^(alpha*nu*exp(-1/2*((ellss2-ellss00)/sigma)^2));

sc          =   chi*sell.^(alpha*(1+nu*exp(-1/2*((sell-ellss00)/sigma).^2)));
sh          =   gammah.*sc/delta;


figure
plot(sell,Well,'-','color',rgb('darkgrey'),'linewidth',4);
hold on
plot(ellss0,0,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',14);
plot(ellss1,0,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',14);
plot(ellss2,0,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',14);
aa=axis;
plot(aa(1:2),[0 0],'--','linewidth',1,'color',rgb('black'));
axis(aa);
plot(ellss0,0,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',14);
plot(ellss1,0,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',14);
plot(ellss2,0,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',14);
set(gca,'FontName','Times','FontSize',16);
xlabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('0 at the SS','interpreter','latex','fontsize',20,'FontName','Times');
hold off

print('Multiple_SS','-dpng');


figure
plot(sh,Well,'-','color',rgb('darkgrey'),'linewidth',4);
hold on
plot(hss0,0,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',14);
plot(hss1,0,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',14);
plot(hss2,0,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',14);
set(gca,'FontName','Times','FontSize',16);
aa=axis;
plot(aa(1:2),[0 0],'--','linewidth',1,'color',rgb('black'));
axis(aa);
plot(hss0,0,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',14);
plot(hss1,0,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',14);
plot(hss2,0,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',14);
xlabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('0 at the SS','interpreter','latex','fontsize',20,'FontName','Times');
hold off





model;
alpha*c/ell*1/(c-gamma*h(-1)+stilde) - omega0*ell^omega +lambda*gammah*alpha*c/ell = 0;
-lambda+beta*(-gamma/(c(+1) - gamma*h+stilde) +(1-delta)*lambda(+1)) = 0;
c = chi*F*ell^alpha;
h-(1-delta)*h(-1)-gammah*c=0;
F=ell^(alpha*nu*exp(-1/2*((ell-ellss00)/sigma)^2));
end; 









%% Deterministic simulations 
% Middle SS

initval;
h       =   hss0;
c       =   css0;
ell    =   ellss0;
F       =   Fss0;
lambda  =   lambdass0;    
end;


steady;
check;


% Low SS

initval;
h       =   hss1;
c       =   css1;
ell    =   ellss1;
F       =   Fss1;
lambda  =   lambdass1;    
end;



steady;
check;
hs=oo_.steady_state(1);
ells=oo_.steady_state(2);
cs=oo_.steady_state(3);
lambdas=oo_.steady_state(4);
Fs=oo_.steady_state(5);



initval;
h       =   hss0*.99999;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

endval;
h       =   hs;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

perfect_foresight_setup(periods=300);
perfect_foresight_solver(stack_solve_algo=0,maxit=20,no_homotopy);

Vell1=oo_.endo_simul(2,2:end)';
Vh1=oo_.endo_simul(1,2:end)';


% High SS

initval;
h       =   hss2;
c       =   css2;
ell    =   ellss2;
F       =   Fss2;
lambda  =   lambdass2;    
end;


steady;
check;


hs=oo_.steady_state(1);
ells=oo_.steady_state(2);
cs=oo_.steady_state(3);
lambdas=oo_.steady_state(4);
Fs=oo_.steady_state(5);



initval;
h       =   hss0*1.00001;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

endval;
h       =   hs;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

perfect_foresight_setup(periods=300);
perfect_foresight_solver(stack_solve_algo=0,maxit=20,no_homotopy);

Vell2=oo_.endo_simul(2,2:end)';
Vh2=oo_.endo_simul(1,2:end)';


%% Figures

figure
subplot(221),plot(Vh1,'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(222),plot(Vh1(1:30),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(223),plot(Vell1,'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(224),plot(Vell1(1:30),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');

print('Sim_SS1','-dpng');



figure
subplot(221),plot(Vh2,'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(222),plot(Vh2(1:30),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(223),plot(Vell2,'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(224),plot(Vell2(1:30),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');

print('Sim_SS2','-dpng');


figure
plot(Vell1(1:20),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold on;
plot(Vell2(1:20),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
plot(Vell2(1:20)./Vell2(1:20),'LineStyle','--','color',rgb('black'),'linewidth',1);
set(gca,'FontName','Times','FontSize',16);
title('(b2) $\nu$ at the bifurcation','interpreter','latex','fontsize',20,'FontName','Times');
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell_t$','interpreter','latex','fontsize',20,'FontName','Times');
axis([1 20 Vell1(20)-.01 Vell2(20)+.01])
pbaspect([3 1 1]);
hold off




figure
plot([hss0*1.00001; Vh1(1:29)],Vell1(1:30),'--o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',10);
hold on
plot([hss0*.99999;Vh2(1:29)],Vell2(1:30),'-o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('darkgray'),...
    'MarkerSize',10);
plot(hss0,ellss0,'o',...
    'linewidth',3,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
plot(hss1,ellss1,'o',...
    'linewidth',3,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
plot(hss2,ellss2,'o',...
    'linewidth',3,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
hold off
set(gca,'FontName','Times','FontSize',16);
xlabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');

Tend=100;

VVh1=[hss0*1.00001;Vh1(1:Tend)];
VVh2=[hss0*0.99999;Vh2(1:Tend)];

figure
plot(VVh1(1:Tend-1),VVh1(2:Tend),'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',10);
hold on
plot([VVh1(1) VVh1(2)],[VVh1(2) VVh1(2)],'--','linewidth',1,'color',rgb('lightgray'));
plot([VVh2(1) VVh2(2)],[VVh2(2) VVh2(2)],'--','linewidth',1,'color',rgb('lightgray'));
for t=2:Tend-1
plot([VVh1(t) VVh1(t+1)],[VVh1(t+1) VVh1(t+1)],'--','linewidth',1,'color',rgb('lightgray'));
plot([VVh1(t) VVh1(t)],[VVh1(t) VVh1(t+1)],'--','linewidth',1,'color',rgb('lightgray'));
plot([VVh2(t) VVh2(t+1)],[VVh2(t+1) VVh2(t+1)],'--','linewidth',1,'color',rgb('lightgray'));
plot([VVh2(t) VVh2(t)],[VVh2(t) VVh2(t+1)],'--','linewidth',1,'color',rgb('lightgray'));
end
plot(VVh1(1:Tend-1),VVh1(2:Tend),'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',10);
plot(VVh2(1:Tend-1),VVh2(2:Tend),'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerSize',10);
plot(hss0,hss0,'o',...
    'linewidth',1,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
plot(hss1,hss1,'o',...
    'linewidth',1,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
plot(hss2,hss2,'o',...
    'linewidth',1,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
set(gca,'FontName','Times','FontSize',16);
xlabel('$h_{t}$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h_{t+1}$','interpreter','latex','fontsize',20,'FontName','Times');
aa=axis;
plot(aa(1:2),aa(3:4),'-','linewidth',1,'color',rgb('black'));
axis(aa);
hold off


% return

% Tracing out the h,h(+1) locus


% Staring from below SS0


htot=[];
htotp1=[];

for hstart=linspace(hss1*.95,hss0*.9999,500)
% for hstart=linspace(hss*.99,hss*1.01,100)


initval;
h       =   hss1;
c       =   css1;
ell    =   ellss1;
F       =   Fss1;
lambda  =   lambdass1;    
end;

steady;
check;
hs=oo_.steady_state(1);
ells=oo_.steady_state(2);
cs=oo_.steady_state(3);
lambdas=oo_.steady_state(4);
Fs=oo_.steady_state(5);


initval;
h       =   hstart;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

endval;
h       =   hs;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

perfect_foresight_setup(periods=100);
perfect_foresight_solver(stack_solve_algo=0,maxit=20,no_homotopy);

Vell1=oo_.endo_simul(2,2:end)';
Vh1=oo_.endo_simul(1,2:end)';

htot=[htot;oo_.endo_simul(1,1:end-1)'];
htotp1=[htotp1;oo_.endo_simul(1,2:end)'];

end





% Staring from above SS0

for hstart=linspace(hss0*1.0001,hss2*1.05,500)
% for hstart=linspace(hss*.99,hss*1.01,100)


initval;
h       =   hss2;
c       =   css2;
ell    =   ellss2;
F       =   Fss2;
lambda  =   lambdass2;    
end;

steady;
check;
hs=oo_.steady_state(1);
ells=oo_.steady_state(2);
cs=oo_.steady_state(3);
lambdas=oo_.steady_state(4);
Fs=oo_.steady_state(5);


initval;
h       =   hstart;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

endval;
h       =   hs;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

perfect_foresight_setup(periods=100);
perfect_foresight_solver(stack_solve_algo=0,maxit=20,no_homotopy);

Vell1=oo_.endo_simul(2,2:end)';
Vh1=oo_.endo_simul(1,2:end)';

htot=[htot;oo_.endo_simul(1,1:end-1)'];
htotp1=[htotp1;oo_.endo_simul(1,2:end)'];

end



figure
plot(htot,htotp1,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('lightgray'),...
    'MarkerEdgeColor',rgb('lightgray'), 'MarkerSize',4);
hold on
plot(hss0,hss0,'o',...
    'linewidth',1,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
plot(hss1,hss1,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
plot(hss2,hss2,'o',...
    'linewidth',1,'color',rgb('dimgray'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
set(gca,'FontName','Times','FontSize',16);
xlabel('$h_{t}$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h_{t+1}$','interpreter','latex','fontsize',20,'FontName','Times');
aa=[min(htot) max(htot) min(htot) max(htot)];
plot(aa(1:2),aa(1:2),'--','linewidth',1,'color',rgb('black'));
axis([aa(1:2) aa(1:2)]);
hold off

print('Phase_Plot','-dpng');








%% Bifurcation plot


xrange=linspace(2.69,2.6906,200);

    inu=0;

for xnu=xrange % for nu
    inu=inu+1;
    nu=xnu;        
    Vnu(inu)=nu;

initval;
h       =   hss0;
c       =   css0;
ell    =   ellss0;
F       =   Fss0;
lambda  =   lambdass0;    
end;

fprintf(['nu = %s'],num2str(nu));
steady;
check;

Veig(:,inu)=oo_.dr.eigval(:,1);

end



figure
plot(Vnu,(Veig(1,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold on
plot(Vnu,(Veig(2,:)),'LineStyle','-','color',rgb('dimgrey'),'linewidth',4);
% plot(Vnu,abs(Veig(3,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
% plot(Vnu,abs(Veig(4,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
plot(Vnu,Vnu./Vnu,'LineStyle','--','color',rgb('black'),'linewidth',1);
% axis([-inf inf 0 4]);
set(gca,'FontName','Times','FontSize',16);
xlabel('$\nu$','interpreter','latex','fontsize',16,'FontName','Times');
ylabel('$\lambda$','interpreter','latex','fontsize',16,'FontName','Times');
hold off

print('Bifurcation_Plot','-dpng');





%% Deterministic simulations with nu = 0;

nu=0;

% Middle SS

initval;
h       =   hss0;
c       =   css0;
ell    =   ellss0;
F       =   Fss0;
lambda  =   lambdass0;    
end;


steady;
check;




initval;
h       =   hss0*1.1;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

endval;
h       =   hs;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

perfect_foresight_setup(periods=300);
perfect_foresight_solver(stack_solve_algo=0,maxit=20,no_homotopy);

Vell1nu0=oo_.endo_simul(2,2:end)';
Vh1nu0=oo_.endo_simul(1,2:end)';



initval;
h       =   hss0*.9;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

endval;
h       =   hs;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

perfect_foresight_setup(periods=300);
perfect_foresight_solver(stack_solve_algo=0,maxit=20,no_homotopy);

Vell2nu0=oo_.endo_simul(2,2:end)';
Vh2nu0=oo_.endo_simul(1,2:end)';


figure
plot(Vell1nu0(1:10),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold on;
plot(Vell2nu0(1:10),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
plot(Vell1nu0(1:10)./Vell1nu0(1:10),'LineStyle','--','color',rgb('black'),'linewidth',1);
set(gca,'FontName','Times','FontSize',16);
title('(b1) $\nu=0$','interpreter','latex','fontsize',20,'FontName','Times');
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell_t$','interpreter','latex','fontsize',20,'FontName','Times');
% axis([1 10 .99 Inf])
pbaspect([3 1 1]);
hold off



