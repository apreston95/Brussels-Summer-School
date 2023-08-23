% EABCN 2023 - The Macroeconomics of Complemetarities Estimate spectrum
% Practice Session II, item 1 (PSII1)
% PSII1b: Simulates the economy with a limit cycle
% Piecewise linear F function


clear all;
close all;
rng('default');
rng(1);
options=optimset('Display','off','MaxFunEvals',10000,'MaxIter',1000);

% Choose the simulation
% if mmd==1  : Transitional dynamics
% if mmd==2  : iid shocks
% if mmd==3  : Persistent shock

mmd=3;


N=270;
N0=50;



% model

if mmd==1  % Transitional dynamics
    sigeps=0.0;
    rho=0;
end

if mmd==2  % iid shocks
    sigeps=0.005;
    rho=0;
end

if mmd==3   % Persistent shock
    sigeps=0.01;
    rho=0.96;
end


alpha1=.03;
delta=.1;
alpha2=.05;



% steady state
Iss=1.;
Xss=Iss/delta;
alpha0= (1+alpha1/delta-alpha2)*Iss;


I1=.5;
I2=1.5;

b0=.01;
b2=.1;
b1=.96;


% Dynamics eignevalues at the SS:
M= [(alpha2-alpha1)/(1-b1) -alpha1*(1-delta)/(1-b1);...
    1 1-delta];
disp('eig(M)');
disp(eig(M));


vI=.01:.01:1.9;
for i=1:length(vI);
    vIi(i)=BigF(vI(i),Iss,I1,I2,b0,b1,b2) + alpha0 -alpha1*Xss + alpha2*Iss;
end
vIi0=vIi;



h=figure;
plot(vI,vIi,'linewidth',3,'color',rgb('dimgray'));
hold on
plot(vI,vI,'--','linewidth',3,'color',rgb('black'));
hold off
xlabel('$I$','FontName','Times','FontSize',20,'interpreter','latex')
ylabel({'$I_i$'},'FontName','Times','FontSize',20,'interpreter','latex')
set(gca,'FontName','Times','FontSize',18)


% Non stochastic case

% initial conditions
I(1)=Iss*1.01;
X(1)=NaN;
X(2)=Xss*1;

epsshock=randn(N,1)*0;
shock(1)=0;



for i=2:N
    shock(i)=rho*shock(i-1)+epsshock(i);
    ct(i)=shock(i)+alpha0-alpha1*X(i)+alpha2*I(i-1);
    I(i)=fsolve('PSII1bx',Iss,options,ct(i),I1,I2,Iss,b0,b1,b2);


    for j=1:length(vI);
        vIi(j)=BigF(vI(j),Iss,I1,I2,b0,b1,b2) + alpha0 -alpha1*X(i) + alpha2*I(i-1);
    end
    X(i+1)=(1-delta)*X(i) + I(i);
end

% spectrum of I - No shocks

% Zero Padding
% padding=0;
padding=512;
% padding=1024;

% Smoothing window
%             = 0 : no smoothing is performed
%             = 1 : Blackman-Tukey window
%             = 2 : Parzen window
%             = 3 : Tukey-Hanning window
SpectrWindow=3;

% Estimate confidence bands
ConfBand=0;
% ConfBand=1;

SpecNS=Spectrum(log(I(2:N)'),padding,ConfBand,SpectrWindow);

% Growth rate autocorrelation

dI=diff(log(I(2:N)'));
[acfNS,xxNS,yyNS]=autocorr(dI,10,[],2);



h=figure;
plot(xxNS(2:end),acfNS(2:end),'linewidth',3,'color',rgb('dimgray'))
hold on;
plot(xxNS(2:end),ones(1,10)*yyNS(1),'linewidth',3,'color',rgb('lightgray'))
plot(xxNS(2:end),ones(1,10)*yyNS(2),'linewidth',3,'color',rgb('lightgray'))
xlabel('Lags','FontName','Times','FontSize',20,'interpreter','latex')
set(gca,'FontName','Times','FontSize',18)
title('Autocorrelation of $I$ growth rate ','FontName','Times','FontSize',20,'interpreter','latex')



%% simulation


NS=N+100;


    disp(j)
    clear I X shock ct R dI


    I(1)=Iss*1.05;
    X(1)=NaN;
    X(2)=Xss*1;


    epsshock=randn(NS,1)*sigeps;
    shock(1)=0;

    for i=2:NS
        shock(i)=rho*shock(i-1)+epsshock(i);
        ct(i)=shock(i)+alpha0-alpha1*X(i)+alpha2*I(i-1);
        I(i)=fsolve('PSII1bx',Iss,options,ct(i),I1,I2,Iss,b0,b1,b2);


        X(i+1)=(1-delta)*X(i) + I(i);
    end


    Is=I(2+99:N+100)';
    Xs=X(2+99:N+100)';



    h=figure;
    plot([Xs(2:N0)'],'linewidth',3,'color',rgb('black'))
    hold on
    plot([ Is(2:N0)'],'linewidth',3,'color',rgb('lightgray'))
    hold off
    xlabel('period','FontName','Times','FontSize',20,'interpreter','latex')
    legend('X','I');
    set(gca,'FontName','Times','FontSize',18)
    
    print('TS_Plot_Persistent','-dpng');


    h=figure;
    plot(Xs(2:N),Is(2:N),'x','markersize',12,'linewidth',3,'color',rgb('dimgray'))
    xlabel('$X$','FontName','Times','FontSize',20,'interpreter','latex')
    ylabel({'$I$'},'FontName','Times','FontSize',20,'interpreter','latex')
    set(gca,'FontName','Times','FontSize',18)
    
    print('Scatter_Plot_Persistent','-dpng');


    % spectrum of I


    Spec=Spectrum(log(Is(2:N)),padding,ConfBand,SpectrWindow);
    
    print('Spectrum_Persistent','-dpng');


    % Growth rate autocorrelation

    dI=diff(log(Is));
    [acf,xx,yy]=autocorr(dI,10,[],2);



    h=figure;
    plot(xx(2:end),acf(2:end),'linewidth',3,'color',rgb('dimgray'))
    hold on;
    plot(xx(2:end),ones(1,10)*yy(1),'linewidth',3,'color',rgb('lightgray'))
    plot(xx(2:end),ones(1,10)*yy(2),'linewidth',3,'color',rgb('lightgray'))
    xlabel('Lags','FontName','Times','FontSize',20,'interpreter','latex')
    set(gca,'FontName','Times','FontSize',18)
    title('Autocorrelation of $I$ growth rate ','FontName','Times','FontSize',20,'interpreter','latex')
    
    print('ACF_Persistent','-dpng');





