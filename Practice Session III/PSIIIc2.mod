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

nu          =   6.4;

% SS for nu=0, to set omega0

ells        =   1;
hs          =   ells^alpha/delta;
lambdas     =   -beta*mu0/(beta*(1-delta)-1) * hs^mu;
omega0      =   alpha*(lambdas+1);

disp('    omega0')
disp(omega0)



% As I am only here looking at indterminate solutions, I lag the first equation by 1 
% so that ell "looks" predetermined, and I choose arbitrarily the value of ell(1)

model;
- (omega0/alpha*ell(-1)^omega -1/ell(-1))*ell(-1)^(1-alpha) + beta*mu0*ell^(alpha*nu)*h^mu + beta*(1-delta)*((omega0/alpha*ell^omega -1/ell)*ell^(1-alpha)) = 0;
h = (1-delta)*h(-1) + ell(-1)^alpha;
end; 



% SS

initval;
ell     =   ells;
h       =   hs;
end;


steady;
check;

hs=oo_.steady_state(2);
ells=oo_.steady_state(1);








% Deterministic simulations when there are sunspots

nu          =   6.4;
hstart       =   hs*.98;



%%%%
initval;
h       =   hstart;
ell    =   ells;
end;

endval;
h       =   hs;
ell    =   ells;
end;

perfect_foresight_setup(periods=1000);
perfect_foresight_solver(stack_solve_algo=0,maxit=20);



Vell1=[oo_.endo_simul(1,1:end)'];
Vh1=[oo_.endo_simul(2,1:end)'];


%%%%
initval;
h       =   hstart;
ell    =   ells*1.01;
end;

endval;
h       =   hs;
ell    =   ells;
end;

perfect_foresight_setup(periods=1000);
perfect_foresight_solver(stack_solve_algo=0,maxit=20);

Vell2=[oo_.endo_simul(1,1:end)'];
Vh2=[oo_.endo_simul(2,1:end)'];

%%%%
initval;
h       =   hstart;
ell    =   ells*.99;
end;

endval;
h       =   hs;
ell    =   ells;
end;

perfect_foresight_setup(periods=1000);
perfect_foresight_solver(stack_solve_algo=0,maxit=50,no_homotopy);

Vell3=[oo_.endo_simul(1,1:end)'];
Vh3=[oo_.endo_simul(2,1:end)'];




Tend=150;

figure
subplot(211),plot(Vh1(1:Tend),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold on;
plot(Vh2(1:Tend),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
plot(Vh3(1:Tend),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold off
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h$','interpreter','latex','fontsize',20,'FontName','Times');

subplot(212),plot(Vell1(1:Tend),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold on;
plot(Vell2(1:Tend),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
plot(Vell3(1:Tend),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold off
set(gca,'FontName','Times','FontSize',16);
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');

print('Det_Sim','-dpng');



figure
plot(Vell1(1:Tend),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold on;
plot(Vell2(1:Tend),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
plot(Vell3(1:Tend),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
plot(Vell2(1:Tend)./Vell2(1:Tend),'LineStyle','--','color',rgb('black'),'linewidth',1);
set(gca,'FontName','Times','FontSize',16);
title('(b2) $\nu$ at the bifurcation','interpreter','latex','fontsize',20,'FontName','Times');
xlabel('Periods','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell_t$','interpreter','latex','fontsize',20,'FontName','Times');
% axis([1 20 Vell1(20)-.01 Vell2(20)+.01])
pbaspect([3 1 1]);
hold off



% return

% Tracing out the h,h(+1) locus



hstart=hs*.99;


initval;
h       =   hstart;
ell    =   ells;
end;

endval;
h       =   hs;
ell    =   ells;
end;

perfect_foresight_setup(periods=1000);
perfect_foresight_solver(stack_solve_algo=0,maxit=20,no_homotopy);


htota=[oo_.endo_simul(2,1:end-1)'];
htotp1a=[oo_.endo_simul(2,2:end)'];




hstart=hs*1.01;

initval;
h       =   hstart;
ell    =   ells;
end;

endval;
h       =   hs;
ell    =   ells;
end;

perfect_foresight_setup(periods=1000);
perfect_foresight_solver(stack_solve_algo=0,maxit=20,no_homotopy);


htotb=[oo_.endo_simul(2,1:end-1)'];
htotp1b=[oo_.endo_simul(2,2:end)'];







figure
plot(htota,htotp1a,'-',...
    'linewidth',6,'color',rgb('dimgray'),'MarkerFaceColor',rgb('dimgray'),...
    'MarkerEdgeColor',rgb('dimgray'), 'MarkerSize',4,...
    'LineStyle','-','color',rgb('darkgrey'));
hold on
plot(htotb,htotp1b,'-',...
    'linewidth',6,'color',rgb('dimgray'),'MarkerFaceColor',rgb('dimgray'),...
    'MarkerEdgeColor',rgb('dimgray'), 'MarkerSize',4,...
    'LineStyle','-','color',rgb('darkgrey'));
plot(hs,hs,'o',...
    'linewidth',1,'color',rgb('black'),'MarkerFaceColor',rgb('black'),...
    'MarkerSize',12);
set(gca,'FontName','Times','FontSize',16);
xlabel('$h_{t}$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$h_{t+1}$','interpreter','latex','fontsize',20,'FontName','Times');
aa=[min([htota;htotb]) max([htota;htotb]) min([htota;htotb]) max([htota;htotb])];
plot(aa(1:2),aa(1:2),'--','linewidth',1,'color',rgb('black'));
axis([aa(1:2) aa(1:2)]);
axis([1.995 2.005 1.995 2.005]);
hold off

print('Phase_Plot_BF','-dpng');


