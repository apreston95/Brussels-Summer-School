% EABCN 2023 - The Macroeconomics of Complemetarities Estimate spectrum
% Practice Session I, item 2 (PSI2)
% PSI2: Heuristic estimation of spectral density


clear all
close all;
clc;

pprint=1;



%% variable and sample

dateQ=[1800:.25:2050]';
SampQ=1:length(dateQ);

TimeAdjust(1)=0;
TimeAdjust(4)=.25;
TimeAdjust(7)=.5;
TimeAdjust(10)=.75;


% Hours Worked
Series1=NaN(length(dateQ),1);
GSeries = getFredData('HOANBS', '1947-01-01', '2023-04-01');
DateSeries= datevec(GSeries.Data(:,1));
Sampseries = find(dateQ==DateSeries(1,1)+TimeAdjust(DateSeries(1,2))):find(dateQ==DateSeries(end,1)+TimeAdjust(DateSeries(end,2)));
Series1(Sampseries) = GSeries.Data(:,2);
% Population
Series2=NaN(length(dateQ),1);
GSeries = getFredData('CNP16OV', '1948-01-01', '2023-07-01','lin','q','avg');
DateSeries= datevec(GSeries.Data(:,1));
Sampseries = find(dateQ==DateSeries(1,1)+TimeAdjust(DateSeries(1,2))):find(dateQ==DateSeries(end,1)+TimeAdjust(DateSeries(end,2)));
Series2(Sampseries) = GSeries.Data(:,2);
% Series specific adjustment
Samp=SampQ(find(dateQ==1948)):SampQ(find(dateQ==2023));
% Hours per capita
Series=log(Series1./Series2);
NameSeries= 'Hours Worked per capita';

var=Series(Samp);
dstart=1948;
dend=2023;





% making the series of even size
if rem(length(var),2)==0
else
    dstart=dstart+.25;
end

ddate=[dstart:.25:dend]';
samp=find(ddate==dstart):find(ddate==dend);
T=length(samp);


%% constructing the frequencies

for j=1:T/2
    omega(j)=2*pi*j/T;
end
period=2*pi./omega;

for j=1:T/2
    for t=1:T
        X(t,2*j-1)=cos(omega(j)*t);
        X(t,2*j)=sin(omega(j)*t);
    end
end

X(:,end)=ones(T,1);

beta=inv(X'*X)*X'*var(samp);

figure;
plot(ddate,X(:,1),'linewidth',4,'color',rgb('darkgrey'))
% xlabel('Time','FontName','Times','FontSize',12)
% ylabel('log','FontName','Times','FontSize',12)
set(gca,'FontName','Times','FontSize',12)
xy=axis;
axis([1948 2022 xy(3) xy(4)]);
pbaspect([2 1 1])
if pprint==1;print -depsc2 cos1;end


figure;
plot(ddate,X(:,1:2),'linewidth',4,'color',rgb('darkgrey'))
% xlabel('Time','FontName','Times','FontSize',12)
% ylabel('log','FontName','Times','FontSize',12)
set(gca,'FontName','Times','FontSize',12)
xy=axis;
axis([1948 2022 xy(3) xy(4)]);
pbaspect([2 1 1])
if pprint==1;print -depsc2 cossin1;end

figure;
plot(ddate,X(:,3),'linewidth',4,'color',rgb('darkgrey'))
% xlabel('Time','FontName','Times','FontSize',12)
% ylabel('log','FontName','Times','FontSize',12)
set(gca,'FontName','Times','FontSize',12)
xy=axis;
axis([1948 2022 xy(3) xy(4)]);
pbaspect([2 1 1])
if pprint==1;print -depsc2 cos2;end

figure;
plot(ddate,X(:,3:4),'linewidth',4,'color',rgb('darkgrey'))
% xlabel('Time','FontName','Times','FontSize',12)
% ylabel('log','FontName','Times','FontSize',12)
set(gca,'FontName','Times','FontSize',12)
xy=axis;
axis([1948 2022 xy(3) xy(4)]);
pbaspect([2 1 1])
if pprint==1;print -depsc2 cossin2;end

figure;
plot(ddate,X(:,5),'linewidth',4,'color',rgb('darkgrey'))
% xlabel('Time','FontName','Times','FontSize',12)
% ylabel('log','FontName','Times','FontSize',12)
set(gca,'FontName','Times','FontSize',12)
xy=axis;
axis([1948 2022 xy(3) xy(4)]);
pbaspect([2 1 1])
if pprint==1;print -depsc2 cos3;end

figure;
plot(ddate,X(:,7),'linewidth',4,'color',rgb('darkgrey'))
% xlabel('Time','FontName','Times','FontSize',12)
% ylabel('log','FontName','Times','FontSize',12)
set(gca,'FontName','Times','FontSize',12)
xy=axis;
axis([1948 2022 xy(3) xy(4)]);
pbaspect([2 1 1])
if pprint==1;print -depsc2 cos4;end

figure;
plot(ddate,X(:,9),'linewidth',4,'color',rgb('darkgrey'))
% xlabel('Time','FontName','Times','FontSize',12)
% ylabel('log','FontName','Times','FontSize',12)
set(gca,'FontName','Times','FontSize',12)
xy=axis;
axis([1948 2022 xy(3) xy(4)]);
pbaspect([2 1 1])
if pprint==1;print -depsc2 cos5;end

figure;
plot(ddate,X(:,11),'linewidth',4,'color',rgb('darkgrey'))
% xlabel('Time','FontName','Times','FontSize',12)
% ylabel('log','FontName','Times','FontSize',12)
set(gca,'FontName','Times','FontSize',12)
xy=axis;
axis([1948 2022 xy(3) xy(4)]);
pbaspect([2 1 1])
if pprint==1;print -depsc2 cos6;end


figure;
plot(ddate,X(:,295),'linewidth',.11,'color',rgb('darkgrey'))
% xlabel('Time','FontName','Times','FontSize',12)
% ylabel('log','FontName','Times','FontSize',12)
set(gca,'FontName','Times','FontSize',12)
xy=axis;
axis([1948 2022 xy(3) xy(4)]);
pbaspect([2 1 1])
if pprint==1;print -depsc2 cos7;end



%% Approximate spectral density

svar=std(var(samp));
for j=1:T/2-1
    sp(j)= (std(beta(2*j-1)*X(:,2*j-1)+beta(2*j)*X(:,2*j))/svar)^2;
end
j=T/2;
sp(j)= (std(beta(2*j-1)*X(:,2*j-1))/svar)^2;

sp=100*sp;

figure;
plot(ddate,var(samp),'linewidth',4,'color',rgb('darkgrey'))
% xlabel('Time','FontName','Times','FontSize',12)
ylabel('log','FontName','Times','FontSize',12)
set(gca,'FontName','Times','FontSize',12)
xy=axis;
axis([1948 2022 xy(3) xy(4)]);
pbaspect([2 1 1])
if pprint==1;print -depsc2 UShourspc;end
% ll=legend('$\pi$','$\ell$');
% set(ll,'FontName','Times','FontSize',32,'interpreter', 'latex')
% set(ll,'color','none');

% figure
% plot(period,sp,'linewidth',4,'color',rgb('darkgrey'))
% xlabel('Periodicity in quarters','FontName','Times','FontSize',18)
% ylabel('Share of total variance','FontName','Times','FontSize',18)
% set(gca,'FontName','Times','FontSize',18)
% 
% figure
% plot(period(4:end),sp(4:end),'linewidth',4,'color',rgb('darkgrey'))
% xlabel('Periodicity in quarters','FontName','Times','FontSize',18)
% ylabel('Share of total variance','FontName','Times','FontSize',18)
% set(gca,'FontName','Times','FontSize',18)
% 
% figure
% plot(period(2:end),sp(2:end),'linewidth',4,'color',rgb('darkgrey'))
% xlabel('Periodicity in quarters','FontName','Times','FontSize',18)
% ylabel('Share of total variance','FontName','Times','FontSize',18)
% set(gca,'FontName','Times','FontSize',18)


figure
bar(period(1:end),sp(1:end),'linewidth',8,'edgecolor',rgb('darkgrey'),'facecolor',rgb('darkgrey'),'barwidth',1)
xlabel('Periodicity in quarters','FontName','Times','FontSize',14,'interpreter','latex')
ylabel('Share of total variance (\%)','FontName','Times','FontSize',14,'interpreter','latex')
set(gca,'FontName','Times','FontSize',14)
axis([2 300 -Inf Inf])
pbaspect([2 1 1])
if pprint==1;print -dpng SpectrumFull;end

figure
bar(period(4:end),sp(4:end),'linewidth',8,'edgecolor',rgb('darkgrey'),'facecolor',rgb('darkgrey'),'barwidth',1)
xlabel('Periodicity in quarters','FontName','Times','FontSize',14,'interpreter','latex')
ylabel('Share of total variance (\%)','FontName','Times','FontSize',14,'interpreter','latex')
set(gca,'FontName','Times','FontSize',14);
axis([2 75 -Inf Inf])
pbaspect([2 1 1])
if pprint==1;print -dpng SpectrumShort;end


%% Fit

fprintf('\n  Per. \t  Index\n')
for i=1:size(period,2)
fprintf('\t %3.1f',period(i));
fprintf('\t %i\n',i);
end
fprintf('\n');


[a b]=sort(-sp);


fprintf('\n Share \t Per. \t Freq. \t Index')
for i=1:size(a,2)
fprintf('\n %2.1f',-a(i));
fprintf('\t %3.1f',period(b((i))));
fprintf('\t %1.3f',omega(b((i))));
fprintf('\t %i',(b((i))));
end
fprintf('\n');

fprintf('\n Share \t Per. \t Freq.')
for i=1:size(a,2)
fprintf('\n % 2.2f',-a(i));
fprintf('\t & %3.1f',period(b((i))));
fprintf('\t & %1.3f\\\\',omega(b((i))));
end
fprintf('\n');

fprintf('\n Share \t Per. ')
for i=1:20
fprintf('\n % 2.2f',-a(i));
fprintf('\t & %3.1f\\\\',period(b((i))));
end
fprintf('\n');

fprintf('\n');

fprintf('\n Share \t Per. \t Share \t Per. ')
for i=1:10
fprintf('\n % 2.2f',-a(i));
fprintf('\t & %3.1f',period(b((i))));
fprintf('\t & % 2.2f',-a(i+10));
fprintf('\t & %3.1f\\\\',period(b((i+10))));
end
fprintf('\n');


% disp('     share   period      freq     index');
% disp([-a' period(b)' omega(b)' b']);
% 


%%
for nfirst=[1:9 20 30 40]

varfit= beta(T)*X(:,T);
for j=1:nfirst
    varfit=varfit + beta(2*b(j)-1)*X(:,2*b(j)-1)+beta(2*b(j))*X(:,2*b(j));
end

figure;
plot(ddate,var(samp),'linewidth',4,'color',rgb('black'))
hold on
plot(ddate,varfit,'linewidth',8,'color',rgb('darkgrey'))
hold off
ylabel('log','FontName','Times','FontSize',12)
% xlabel('Time','FontName','Times','FontSize',12)
set(gca,'FontName','Times','FontSize',12)
xy=axis;
axis([1948 2022 xy(3) xy(4)]);
pbaspect([2 1 1])
if nfirst==1
 Raxis=axis;   
end
axis(Raxis);
if nfirst==1
 Raxis=axis;   
end
axis(Raxis);
if nfirst==1
eval(['title(''With the First Main Frequency'');']);
else    
eval(['title(''With the ' num2str(nfirst) ' Main Frequencies'');']);
end
if pprint==1;eval(['print -dpng fnfirst' num2str(nfirst) ';']);end;
end



%%
% 
% % deviations from the 1st contributor
% 
% varfit1= beta(T)*X(:,T);
% for j=1
%     varfit1=varfit1 + beta(2*b(j)-1)*X(:,2*b(j)-1)+beta(2*b(j))*X(:,2*b(j));
% end
% 
% figure;
% plot(ddate,var,'linewidth',4,'color',rgb('black'))
% hold on
% plot(ddate,varfit1,'linewidth',8,'color',rgb('darkgrey'))
% hold off
% eval(['title('' first contibutor'');']);
% % xlabel('Time','FontName','Times','FontSize',12)
% ylabel('log','FontName','Times','FontSize',12)
% set(gca,'FontName','Times','FontSize',12)
% 
% 
% var2=(var-varfit1)*100;
% 
% varfit2= zeros(T,1);
% for j=3:5
%     varfit2=varfit2 + (beta(2*b(j)-1)*X(:,2*b(j)-1)+beta(2*b(j))*X(:,2*b(j)))*100;
% end
% 
% figure;
% plot(ddate,var2,'linewidth',4,'color',rgb('black'))
% hold on
% plot(ddate,varfit2,'linewidth',8,'color',rgb('darkgrey'))
% hold off
% title(' XXX');
% % xlabel('Time','FontName','Times','FontSize',12)
% ylabel('% deviations from first contributor','FontName','Times','FontSize',12)
% set(gca,'FontName','Times','FontSize',12)

%% Band pass filter
% 
% % 6-32Q : j=9-48
% varfit632= zeros(T,1);
% for j=9:48
%     varfit632=varfit632 + (beta(2*j-1)*X(:,2*j-1)+beta(2*j)*X(:,2*j))*100;
% end
% 
% figure
% plot_nber(varfit632,dstart,dend);
% xy=axis;
% axis([1948 2022 xy(3) xy(4)]);
% pbaspect([2 1 1])
% if pprint==1;print -depsc2 BP632;end
% 
% varhp=(var(samp)-hpfilter(var(samp),1600))*100;
% 
% % figure
% % plot([varhp varfit632])
% figure
% plot_nber_2series(varhp,varfit632,dstart,dend,'HP','BP(6,32)','NorthEast');
% xy=axis;
% axis([1948 2022 xy(3) xy(4)]);
% pbaspect([2 1 1])
% if pprint==1;print -depsc2 BP632HP;end
% 
% 
% % 20-50Q : j=6-12
% varfit2050= zeros(T,1);
% for j=[6:14]
% % for j=[6]
%     varfit2050=varfit2050 + (beta(2*j-1)*X(:,2*j-1)+beta(2*j)*X(:,2*j))*100;
% end
% 
% figure
% plot_nber(varfit2050,dstart,dend);
% xy=axis;
% axis([1948 2022 xy(3) xy(4)]);
% pbaspect([2 1 1])
% if pprint==1;print -depsc2 BP2050;end
% 
% figure
% plot_nber_2series(varhp,varfit2050,dstart,dend,'HP','BP(20,50)','NorthEast');
% xy=axis;
% axis([1948 2022 xy(3) xy(4)]);
% pbaspect([2 1 1])
% if pprint==1;print -depsc2 BP2050HP;end
% 
% 

%% AR(1)

Y = var(samp);
Y = Y - mean(Y);
Y1 = Y(2:end) ;
T = length(Y1);
X = [ones(T, 1) Y(1:end-1) ];
b = regress(Y1,X)
rho = b(2);
eps = Y1 - X*b;
sigma = std(eps);

% Simulate an AR(1) with the same rho and sigma
T = 300;

Y_Sim(1) = 0;

for i = 2:T
    Y_Sim(i) = rho*Y_Sim(i-1) + sigma*normrnd(0,1);
end

Y_Sim = Y_Sim';

% Now estimate the heuristic spectral density for the AR(1)

for j=1:T/2
    omega(j)=2*pi*j/T;
end
period=2*pi./omega;

for j=1:T/2
    for t=1:T
        X(t,2*j-1)=cos(omega(j)*t);
        X(t,2*j)=sin(omega(j)*t);
    end
end

X(:,end)=ones(T,1);

beta=inv(X'*X)*X'*Y_Sim;

svar=std(Y_Sim);
for j=1:T/2-1
    sp(j)= (std(beta(2*j-1)*X(:,2*j-1)+beta(2*j)*X(:,2*j))/svar)^2;
end
j=T/2;
sp(j)= (std(beta(2*j-1)*X(:,2*j-1))/svar)^2;

sp=100*sp;

figure
bar(period(1:end),sp(1:end),'linewidth',8,'edgecolor',rgb('darkgrey'),'facecolor',rgb('darkgrey'),'barwidth',1)
xlabel('Periodicity in quarters','FontName','Times','FontSize',14,'interpreter','latex')
ylabel('Share of total variance (\%)','FontName','Times','FontSize',14,'interpreter','latex')
set(gca,'FontName','Times','FontSize',14)
axis([2 300 -Inf Inf])
pbaspect([2 1 1])
if pprint==1;print -dpng SpectrumAR1;end







