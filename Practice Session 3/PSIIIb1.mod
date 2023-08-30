% Franck Portier March 2023


close all;


var h ell c lambda F;    

parameters ellss beta gamma omega0 omega delta alpha chi nu sigma stilde gammah;



ells=1;
ellss=ells;
beta=.98;
gamma = -1/2;
gammah=1;
omega=.5;
delta=.025;
alpha=.7;
nu=1.147405;
% nu=1.1536595;
sigma=.1;
chi=1;

stilde= max(-1.1*(delta-gamma*gammah)/delta*chi,1);


% exponent on ell in the production function  
xell=[.5:.01:1.5];
xexpo=alpha*(1+nu*exp(-1/2*((xell-ellss)/sigma).^2));

figure
plot(xell,xexpo,'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('Exponent','interpreter','latex','fontsize',20,'FontName','Times');



% SS

Fs=ells^(alpha*nu);
cs=chi*Fs*ells^(alpha);
hs=gammah*cs/delta;
lambdas = beta*gamma / ((delta/gammah-gamma)*hs +stilde) *1/(beta*(1-delta)-1); 

omega0=alpha*cs/(cs-gamma*hs+stilde) + alpha*gammah*lambdas*cs  ;


disp('     hs        ells    lambdas     cs      omega0');
disp([hs ells lambdas cs omega0]);


model;
alpha*c/ell*1/(c-gamma*h(-1)+stilde) - omega0*ell^omega +lambda*gammah*alpha*c/ell = 0;
-lambda+beta*(-gamma/(c(+1) - gamma*h+stilde) +(1-delta)*lambda(+1)) = 0;
c = chi*F*ell^alpha;
h-(1-delta)*h(-1)-gammah*c=0;
F=ell^(alpha*nu*exp(-1/2*((ell-ellss)/sigma)^2));
end; 


initval;
h       =   hs;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;



steady;
check;


% Deterministic simulation in a case where nu is large so that there are 
% complex eigenvalues, but the steady state is stable

hs=oo_.steady_state(1);
ells=oo_.steady_state(2);
cs=oo_.steady_state(3);
lambdas=oo_.steady_state(4);
Fs=oo_.steady_state(5);

ellss=ells;



initval;
h       =   hs*1.001;
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

perfect_foresight_setup(periods=3000);
perfect_foresight_solver(stack_solve_algo=0,maxit=20,no_homotopy);

Vell=oo_.endo_simul(2,2:end)';
Vh=oo_.endo_simul(1,2:end)';



figure
subplot(221),plot(Vh,'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(222),plot(Vh(1:30),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(223),plot(Vell,'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(224),plot(Vell(1:30),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');







%% Bifurcation plot


xrange=linspace(1.145,1.14729,200);

    inu=0;

for xnu=xrange % for nu
    inu=inu+1;
    nu=xnu;        
    Vnu(inu)=nu;

initval;
h       =   hs;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;

fprintf(['nu = %s'],num2str(nu));
steady;
check;

Veig(:,inu)=oo_.dr.eigval(:,1);

end



figure
plot(Vnu,(Veig(1,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold on
plot(Vnu,[Veig(2,1:158) NaN(1,42)],'LineStyle','-','color',rgb('dimgrey'),'linewidth',4);
plot(Vnu,[NaN(1,158) Veig(2,159:200)],'LineStyle','-','color',rgb('dimgrey'),'linewidth',4);
% plot(Vnu,abs(Veig(3,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
% plot(Vnu,abs(Veig(4,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
plot(Vnu,Vnu./Vnu,'LineStyle','--','color',rgb('black'),'linewidth',1);
plot(Vnu,-Vnu./Vnu,'LineStyle','--','color',rgb('black'),'linewidth',1);
axis([-inf inf -5 5]);
set(gca,'FontName','Times','FontSize',16);
xlabel('$\nu$','interpreter','latex','fontsize',16,'FontName','Times');
ylabel('$\lambda$','interpreter','latex','fontsize',16,'FontName','Times');
hold off




%% Deterministic simulations with nu = 0;

nu=0;

% Middle SS

initval;
h       =   hs;
c       =   cs;
ell    =   ells;
F       =   Fs;
lambda  =   lambdas;    
end;


steady;
check;




initval;
h       =   hs*1.1;
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

Vellnu0=oo_.endo_simul(2,2:end)';
Vhnu0=oo_.endo_simul(1,2:end)';

figure
plot(Vellnu0(1:60),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold on;
plot(Vellnu0(1:60)./Vellnu0(1:60),'LineStyle','--','color',rgb('black'),'linewidth',1);
set(gca,'FontName','Times','FontSize',16);
title('(b1) $\nu=0$','interpreter','latex','fontsize',20,'FontName','Times');
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell_t$','interpreter','latex','fontsize',20,'FontName','Times');
axis([1 60 .93 1.01])
pbaspect([3 1 1]);
hold off


