% Model that produces sunspots
% use with dynare
% Franck Portier April 2023


close all;


var ell  h;    

parameters beta omega0 omega alpha mu0 mu  nu delta;


beta        =   0.95;
omega       =   0.5;
alpha       =   0.7;
delta       =   0.5;
mu0         =   0.1;
mu          =   1;

% nu          =   0;
% nu          =   6.2;
nu          =   6.5;

% SS for nu=0, to set omega0

ells        =   1;
hs          =   ells^alpha/delta;
lambdas     =   -beta*mu0/(beta*(1-delta)-1) * hs^mu;
omega0      =   alpha*(lambdas+1);

disp('    omega0')
disp(omega0)



model;
- (omega0/alpha*ell^omega -1/ell)*ell^(1-alpha) + beta*mu0*ell(+1)^(alpha*nu)*h^mu + beta*(1-delta)*((omega0/alpha*ell(+1)^omega -1/ell(+1))*ell(+1)^(1-alpha)) = 0;
h = (1-delta)*h(-1) + ell^alpha;
end; 



% SS

initval;
ell     =   ells;
h       =   hs;
end;


steady;
check;




%% Bifurcation plot


xrange=linspace(4,8,200);

    inu=0;

for xnu=xrange % for nu
    inu=inu+1;
    nu=xnu;        
    Vnu(inu)=nu;

initval;
h       =   hs;
ell    =   ells;
end;

fprintf(['nu = %s'],num2str(nu));
steady;
check;

Veig(:,inu)=oo_.dr.eigval(:,1);

end



figure
plot(Vnu,(Veig(1,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold on
plot(Vnu,Veig(2,:),'LineStyle','-','color',rgb('dimgrey'),'linewidth',4);
% plot(Vnu,abs(Veig(3,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
% plot(Vnu,abs(Veig(4,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
plot(Vnu,Vnu./Vnu,'LineStyle','--','color',rgb('black'),'linewidth',1);
plot(Vnu,-Vnu./Vnu,'LineStyle','--','color',rgb('black'),'linewidth',1);
axis([-inf inf 0 2]);
set(gca,'FontName','Times','FontSize',16);
xlabel('$\nu$','interpreter','latex','fontsize',16,'FontName','Times');
ylabel('$\lambda$','interpreter','latex','fontsize',16,'FontName','Times');
hold off

print('Eigvalues','-dpng');





%% Deterministic simulations with nu = 0;

nu=0;

% Middle SS

initval;
h       =   hs;
ell    =   ells;
end;


steady;
check;




initval;
h       =   hs*1.1;
ell    =   ells;
end;

endval;
h       =   hs;
ell    =   ells;
end;

perfect_foresight_setup(periods=300);
perfect_foresight_solver(stack_solve_algo=0,maxit=20,no_homotopy);

Vellnu0=oo_.endo_simul(1,2:end)';

figure
plot(Vellnu0(1:20),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold on;
plot(Vellnu0(1:20)./Vellnu0(1:20),'LineStyle','--','color',rgb('black'),'linewidth',1);
set(gca,'FontName','Times','FontSize',16);
title('(b1) $\nu=0$','interpreter','latex','fontsize',20,'FontName','Times');
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell_t$','interpreter','latex','fontsize',20,'FontName','Times');
% axis([1 60 .93 1.01])
pbaspect([3 1 1]);
hold off



