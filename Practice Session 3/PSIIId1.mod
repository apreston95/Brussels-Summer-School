% This code writes down a model that nests complementarity, 
% habit and durables) and
% solves the eigenvalues for various level of externality
% here I am putting c and not h as an externality for habit 
% Franck Portier October 2020


close all;

var d h c ell lambdad lambdah s;    

parameters nu alpha psi chi omega omega0 beta deltad deltah gammad gammah gamma stilde mu sigma;


mmod = 1;   %: 1: LC, 2: Hysteresis, 3: Sunspots



if mmod==1;
    % LC
    % If loop over nu,
    xrange=[1.8:.0001:1.91];

    chi     =   1;     % scale parameter for externality in durable sector
    mu      =   1;
    beta    =   .98;    % discount factor
    alpha   =   .7;    % Returns to scale to labor in durable sector
    omega   =   1;    % Disutility of labor
    gamma   =   .5;
    stilde  =   1;
    % nu      =   0;
    % nu      =   1.886345;
    nu      =   1.889;
    deltad   =   .025;   % depreciation
    deltah   =   .5;   % depreciation
    gammad  =   .4;
    gammah  =   .4;
    psi     =   .5;

    sigma   =   4;
end

if mmod==2;
    % Hysteresis
    % If loop over nu,
    xrange=[0:.001:1];

    chi     =   1;     % scale parameter for externality in durable sector
    beta    =   .98;    % discount factor
    alpha   =   .7;    % Returns to scale to labor in durable sector
    omega   =   .5;    % Disutility of labor
    gamma   =   1;
    stilde  =   1;
    % nu      =   0;
    % nu      =   1.886345;
    nu      =   .6;
    mu      =   1;
    deltad   =   .025;   % depreciation
    deltah   =   .05;   % depreciation
    gammad  =   .05;
    gammah  =   .08;
    psi     =   .1;

    sigma   =   4;
end;


if mmod==3;
    % Sunspots
    % If loop over nu,
    xrange=[1.06:.0001:1.19];
    
    chi     =   1;     % scale parameter for externality in durable sector
    beta    =   .98;    % discount factor
    alpha   =   .7;    % Returns to scale to labor in durable sector
    omega   =   .5;    % Disutility of labor
    gamma   =   1;
    stilde  =   1;
    % nu      =   0;
    % nu      =   1.886345;
    nu      =   .6;
    mu      =   .5;
    deltad   =   .025;   % depreciation
    deltah   =   .5;   % depreciation
    gammad  =   .4;
    gammah  =   .4;
    psi     =   .5;

    sigma   =   4;

    
end;

% exponent on ell in the production function  
ellss=1;
xell=[.5:.01:1.5];
xexpo=alpha*(1+nu*exp(-1/2*((xell-ellss)/sigma).^2));
xexter=xell.^xexpo;


figure
plot(xell,xexpo,'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\nu$ in $\ell^{(1+\nu)\alpha}$','interpreter','latex','fontsize',20,'FontName','Times');

figure
plot(xell,xexter,'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
set(gca,'FontName','Times','FontSize',16);
xlabel('$\ell$','interpreter','latex','fontsize',20,'FontName','Times');
ylabel('$\ell^{(1+\nu)\alpha}$','interpreter','latex','fontsize',20,'FontName','Times');
hold off;




model;
% # s= psi*d(-1)+c-(1-psi)*gamma*h(-1)+stilde;
1/s * alpha*c/ell - omega0*ell^omega +alpha*c/ell*(gammad*lambdad+gammah*lambdah)=0;
-lambdad+beta*(psi/s(+1) +(1-deltad)*lambdad(+1))=0;
-lambdah+beta*(-(1-psi)*gamma*mu*h^(mu-1)*c(+1)^(1-mu)/s(+1) +(1-deltah)*lambdah(+1))=0;
h-(1-deltah)*h(-1)-gammah*c=0;
d-(1-deltad)*d(-1)-gammad*c=0;
-s+psi*d(-1)+c-(1-psi)*gamma*h(-1)^mu*c^(1-mu)+stilde=0;
c - chi*ell^(alpha*(1+nu*exp(-1/2*((ell-1)/sigma)^2)))=0;
end; 


%%

% compute SS
ellss=1;
css=chi;
dss=gammad/deltad*css;
hss=gammah/deltah*css;
sss=psi*dss-(1-psi)*gamma*hss^mu*css^(1-mu) + css+stilde;
lambdadss=1/(1-beta*(1-deltad)) *beta*psi/sss;
lambdahss=1/(1-beta*(1-deltah)) *(- beta*(1-psi)*gamma*mu*hss^(mu-1)*css^(1-mu)/sss);
omega0= alpha*chi*(1/sss+gammad*lambdadss+gammah*lambdahss);


disp('     dss       hss       ellss     sss     lambdadss  lambdahss  css');
disp([dss hss ellss  sss lambdadss lambdahss css]);
disp('    omega0');
disp([omega0]);




initval;
d       =   dss;
h       =   hss;
c       =   css;
ell    =   ellss;
lambdad  =   lambdahss;    
lambdah  =   lambdahss;    
s       =   sss;
end;



steady;
check;

% var d h c ell lambdah lambdad  s;    


ds=oo_.steady_state(1);
hs=oo_.steady_state(2);
cs=oo_.steady_state(3);
ells=oo_.steady_state(4);
lambdads=oo_.steady_state(5);
lambdahs=oo_.steady_state(6);
ss=oo_.steady_state(7);

ellss=ells;

disp('     dss       hss       ellss     sss     lambdadss  lambdahss  css');
disp([ds hs ells  ss lambdads lambdahs cs]);



%% Loop over nu


inu=0;

for xnu=xrange % for nu
    inu=inu+1;
    nu=xnu;        
    Vnu(inu)=nu;
    fprintf(['nu = %s'],num2str(nu));


    steady;
    check;

    Vells(inu)=oo_.steady_state(4);

    Veig(:,inu)=oo_.dr.eigval(:,1);

end

th = 0:pi/50:2*pi;
xunit = cos(th);
yunit = sin(th);


VVeig=Veig;

for i=1:size(Veig,2)
[sv si]=sort(abs(VVeig(:,i)));
Veig(:,i)=VVeig(si,i);
end

Reig=real(Veig);
Ieig=imag(Veig);




figure
plot(xunit, yunit,'linewidth',5,'color',rgb('lightgray'));
hold on
plot(Reig(4,:), Ieig(4,:),'.','linewidth',2,'markersize',8,'color',rgb('SlateGray'));
plot(Reig(3,:), Ieig(3,:),'.','linewidth',1,'markersize',8,'color',rgb('black'));
plot(Reig(2,:), Ieig(2,:),'.','linewidth',2,'markersize',8,'color',rgb('SlateGray'));
plot(Reig(1,:), Ieig(1,:),'.','linewidth',1,'markersize',8,'color',rgb('black'));
set(gca,'FontName','Times','FontSize',16);
xlabel('Re($\lambda$)','interpreter','latex','fontsize',16,'FontName','Times');
ylabel('Im($\lambda$)','interpreter','latex','fontsize',16,'FontName','Times');
set(gca,'DataAspectRatio',[1,1,1])
xl=xlim;
yl=ylim;
line([0 0],[-2 2],'LineStyle','--','color',rgb('lightgray'),'linewidth',2);
line([-2 2],[0 0],'LineStyle','--','color',rgb('lightgray'),'linewidth',2);
plot(Reig(4,:), Ieig(4,:),'.','linewidth',1,'markersize',18,'color',rgb('SlateGray'));
plot(Reig(3,:), Ieig(3,:),'.','linewidth',1,'markersize',18,'color',rgb('black'));
plot(Reig(2,:), Ieig(2,:),'.','linewidth',1,'markersize',18,'color',rgb('SlateGray'));
plot(Reig(1,:), Ieig(1,:),'.','linewidth',1,'markersize',18,'color',rgb('black'));
plot(Reig(1:4,1), Ieig(1:4,1),'p','linewidth',1,'markersize',12,'color','black','MarkerFaceColor',rgb('black'));
plot(Reig(1:4,end), Ieig(1:4,end),'o','linewidth',1,'markersize',12,'color',rgb('SlateGray'),'MarkerFaceColor',rgb('SlateGray'));
axis(2*[-1 1 -1 1]);

if mmod==2



figure
plot(xunit, yunit,'linewidth',5,'color',rgb('lightgray'));
hold on
plot(Reig(4,:), Ieig(4,:),'.','linewidth',2,'markersize',8,'color',rgb('SlateGray'));
plot(Reig(3,:), Ieig(3,:),'.','linewidth',1,'markersize',8,'color',rgb('black'));
plot(Reig(2,:), Ieig(2,:),'.','linewidth',2,'markersize',8,'color',rgb('SlateGray'));
plot(Reig(1,:), Ieig(1,:),'.','linewidth',1,'markersize',8,'color',rgb('black'));
set(gca,'FontName','Times','FontSize',16);
xlabel('Re($\lambda$)','interpreter','latex','fontsize',16,'FontName','Times');
ylabel('Im($\lambda$)','interpreter','latex','fontsize',16,'FontName','Times');
set(gca,'DataAspectRatio',[1,1,1])
xl=xlim;
yl=ylim;
line([0 0],[-2 2],'LineStyle','--','color',rgb('lightgray'),'linewidth',2);
line([-2 2],[0 0],'LineStyle','--','color',rgb('lightgray'),'linewidth',2);
plot(Reig(4,:), Ieig(4,:),'.','linewidth',1,'markersize',18,'color',rgb('SlateGray'));
plot(Reig(3,:), Ieig(3,:),'.','linewidth',1,'markersize',18,'color',rgb('black'));
plot(Reig(2,:), Ieig(2,:),'.','linewidth',1,'markersize',18,'color',rgb('SlateGray'));
plot(Reig(1,:), Ieig(1,:),'.','linewidth',1,'markersize',18,'color',rgb('black'));
plot(Reig(1:4,1), Ieig(1:4,1),'p','linewidth',1,'markersize',12,'color','black','MarkerFaceColor',rgb('black'));
plot(Reig(1:4,end), Ieig(1:4,end),'o','linewidth',1,'markersize',12,'color',rgb('SlateGray'),'MarkerFaceColor',rgb('SlateGray'));
axis([.96 1.07 -.1 .1]);

end


figure
plot(Vnu,abs(Veig(1,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold on
plot(Vnu,abs(Veig(2,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
plot(Vnu,abs(Veig(3,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
plot(Vnu,abs(Veig(4,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
plot(Vnu,Vnu./Vnu,'LineStyle','--','color',rgb('black'),'linewidth',1);
% axis([-inf inf 0 4]);
if mmod==1
axis([Vnu(1) Vnu(end) -inf inf]);
end
if mmod==2
axis([Vnu(1) Vnu(end) .97 1.05]);
end
if mmod==3
axis([Vnu(1) Vnu(end) .6 2.4]);
end
set(gca,'FontName','Times','FontSize',16);
xlabel('$\nu$','interpreter','latex','fontsize',16,'FontName','Times');
ylabel('$|\lambda|$','interpreter','latex','fontsize',16,'FontName','Times');
hold off



