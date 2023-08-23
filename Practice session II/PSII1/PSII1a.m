% EABCN 2023 - The Macroeconomics of Complemetarities Estimate spectrum
% Practice Session II, item 1 (PSII1)
% PSII1a: Study the evolution of eigenvalues as a function of strategic
% interactions
% Piecewise linear F function

clear all;
close all;

%  Type of bifurcation
 bifurc = 1;   % Hopf
% bifurc = 2;   % Flip
% bifurc = 3;   % Fold


if bifurc == 1
    alpha1=.03;
    delta=.1;
    alpha2=.05;
end

if bifurc == 2
    alpha1=.03;
    delta=.1;
    alpha2=.0006;
end

if bifurc == 3
    alpha1=.025;
    delta=.2;
    alpha2=.6312;
end




% steady state
Iss=1;
Xss=Iss/delta;
alpha0= (1+alpha1/delta-alpha2)*Iss;

% Shape of the F function
% it is required that I1<Iss<I2
% it has slope b0 below I1, b1 between I1 and I2 (that inclides the SS) and
% b2 after  
% F(I) = a0+ b0 I if I<I1
%        a1+ b1 I if I1<I<I2
%        a2+ b2 I if I>I2
% Given b0, b1, b2 we find a0,a1 and a2 by imposing:
% F(Iss)=0
% Continuity at I1 and I2  

I1=.5;
I2=1.5;
b0=.01;
b2=.1;
b1=.6;  % b1 is the slope of the best response function at the SS







% Function BigF gives F(I)  
vI=.01:.01:2;


for i=1:length(vI);
    vIi(i)=BigF(vI(i),Iss,I1,I2,b0,b1,b2) + alpha0 -alpha1*Xss + alpha2*Iss;
end
vIi0=vIi;

for i=1:length(vI);
    vFi(i)=BigF(vI(i),Iss,I1,I2,b0,b1,b2);
end

h=figure;
plot(vI,vFi,'linewidth',3,'color',rgb('gray'));
hold off
xlabel('$I$','FontName','Times','FontSize',20,'interpreter','latex')
ylabel({'$F(I)$'},'FontName','Times','FontSize',20,'interpreter','latex')
xline(1,'--')
set(gca,'FontName','Times','FontSize',18)

print('F_Plot','-dpng');

h=figure;
plot(vI,vIi,'linewidth',3,'color',rgb('dimgray'));
hold on
plot(vI,vI,'--','linewidth',3,'color',rgb('black'));
hold off
xlabel('$I$','FontName','Times','FontSize',20,'interpreter','latex')
ylabel({'$I_i$'},'FontName','Times','FontSize',20,'interpreter','latex')
set(gca,'FontName','Times','FontSize',18)

print('BR_Plot','-dpng');



%% Loop over b1
if bifurc == 1
    Vb1=0:.001:.999;
end
if bifurc == 2
    Vb1=0.9:.0001:.999;
end
if bifurc == 3
    Vb1=0.45:.000001:.55;
end

ib1=0;
for b1=Vb1
    ib1=ib1+1;
    % Eignevalues at the SS:
    M= [(alpha2-alpha1)/(1-b1) -alpha1*(1-delta)/(1-b1);...
        1 1-delta];


    Veig(:,ib1)=eig(M);

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
plot(Reig(2,:), Ieig(2,:),'.','linewidth',1,'markersize',18,'color',rgb('SlateGray'));
plot(Reig(1,:), Ieig(1,:),'.','linewidth',1,'markersize',18,'color',rgb('black'));
plot(Reig(1:2,1), Ieig(1:2,1),'p','linewidth',1,'markersize',12,'color','black','MarkerFaceColor',rgb('black'));
plot(Reig(1:2,end), Ieig(1:2,end),'o','linewidth',1,'markersize',12,'color',rgb('SlateGray'),'MarkerFaceColor',rgb('SlateGray'));
axis(2*[-1 1 -1 1]);

print('Eig_Plot','-dpng');


figure
plot(Vb1,abs(Veig(1,:)),'LineStyle','-','color',rgb('darkgrey'),'linewidth',4);
hold on
plot(Vb1,abs(Veig(2,:)),'LineStyle','-','color',rgb('black'),'linewidth',4);
plot(Vb1,Vb1./Vb1,'LineStyle','--','color',rgb('black'),'linewidth',1);
axis([Vb1(1) Vb1(end) 0 2]);
set(gca,'FontName','Times','FontSize',16);
xlabel('$b_1$','interpreter','latex','fontsize',16,'FontName','Times');
ylabel('$|\lambda|$','interpreter','latex','fontsize',16,'FontName','Times');
hold off

print('SpecRadius_Plot','-dpng');





