function  Spec=Spectrum(Series,padding,ConfBand,SpectrWindow)



Series=Series-mean(Series);

% Spectral density using SPECTRAN

WindowLag=[];
if ConfBand == 0
SpectrAlpha=[];
else
SpectrAlpha=[.34];
end

Lsamp=length(Series);


if padding~=0
    Tpad=padding-Lsamp;
    Nsamp=[1:Lsamp Lsamp+1:Lsamp+Tpad];
    Series(Lsamp+1:Lsamp+Tpad)=0;
end

Spec = periodg(Series,SpectrWindow,WindowLag,SpectrAlpha);


% Range of freqiencies to darken on figures
pa=32;
pb=50;

frmn = 1/80;            % lowest frequency to plot
frmx = 1/4;             % highest frequency to plot

xtickvec = [4, 6, 24, pa, 40, pb, 60, 80, 160, 200];   % vector of xtick mark locations
ntick = numel(xtickvec);                                % number of xtick marks
xticklab = cell(1,ntick);                               % cell to hold xtick labels
for j = 1:ntick                                         % for each tick location
    xticklab{j} = num2str(xtickvec(j),'%i');                % label = string value of xtick mark
end
xtickvec = log(xtickvec);                               % get logs of tick marks

freqs = Spec.frq/(2*pi);                % vector of (ordinary) frequencies
xl = [frmn,frmx];                           % min and max of x-axis in frequencies
xax = log(min(1./freqs,realmax));           % get log of periodicities associated with frequencies
xl = log(1./min(xl([2,1]),realmax));        % min and max of x-axis in log of periodicity


cmap = 'bone';                              % color map name for plot
eval(['cmapvals = ' cmap '(20);']);         % color map numeric values

fsz = 16;               % font size for axis/legend
lnwd = 2;               % line width




h=figure;
plot(xax,Spec.f,'-','LineWidth',4,'color',rgb('black'));
hold on;
if ConfBand==1
    eval(['plot(xax,Spec.fconf,''--'',''LineWidth'',1,''color'',rgb(''black''));']);  
end
xlim(xl)
yy=ylim;
area(([xtickvec(4) xtickvec(2) ]),.995*[yy(2) yy(2) ], 'facecolor',rgb('whitesmoke'), 'edgecolor',rgb('whitesmoke'));
area(([xtickvec(6) xtickvec(4)]) ,.995*[yy(2) yy(2) ], 'facecolor',rgb('gainsboro'), 'edgecolor',rgb('gainsboro'));
yy=ylim;
for i=2:6
    plot([xtickvec(i) xtickvec(i)],[0 yy(2)],'--','color',rgb('dimgray'),'LineWidth',.5);
end
xlabel('Periodicity','FontName','Times','FontSize',fsz)
set(gca,'XTick',xtickvec,'XTickLabel',xticklab,'FontName','Times','FontSize',fsz,'xgrid','on')
plot(xax,Spec.f,'-','LineWidth',4,'color',rgb('black'));
if ConfBand==1
    plot(xax,Spec.fconf,'--','LineWidth',1,'color',rgb('black'));
end
hold off


