% EABCN 2023 - The Macroeconomics of Complemetarities Estimate spectrum
% Practice Session I, item 3 (PSI3)
% PSI3: Estimate spectrum

% Start with unemployment, no padding, no smoothing
% Then play with padding and smoothing
% Reasonable choice: padding=512, smoothing=Tukey-Hanning window



clear all;
close all;

addpath spectran24_11_12_13

dateQ=[1800:.25:2050]';
SampQ=1:length(dateQ);

TimeAdjust(1)=0;
TimeAdjust(4)=.25;
TimeAdjust(7)=.5;
TimeAdjust(10)=.75;



% % Unemployment
NameSeries='Unemployment';
Series=NaN(length(dateQ),1);
GSeries = getFredData('UNRATE', '1948-01-01', '2023-07-01','lin','q','avg');
DateSeries= datevec(GSeries.Data(:,1));
Sampseries = find(dateQ==DateSeries(1,1)+TimeAdjust(DateSeries(1,2))):find(dateQ==DateSeries(end,1)+TimeAdjust(DateSeries(end,2)));
Series(Sampseries) = GSeries.Data(:,2);
% Series specific adjustment
Samp=Sampseries(1):Sampseries(end-1);


% % Hours Worked
% Series1=NaN(length(dateQ),1);
% GSeries = getFredData('HOANBS', '1947-01-01', '2023-04-01');
% DateSeries= datevec(GSeries.Data(:,1));
% Sampseries = find(dateQ==DateSeries(1,1)+TimeAdjust(DateSeries(1,2))):find(dateQ==DateSeries(end,1)+TimeAdjust(DateSeries(end,2)));
% Series1(Sampseries) = GSeries.Data(:,2);
% % Population
% Series2=NaN(length(dateQ),1);
% GSeries = getFredData('CNP16OV', '1948-01-01', '2023-07-01','lin','q','avg');
% DateSeries= datevec(GSeries.Data(:,1));
% Sampseries = find(dateQ==DateSeries(1,1)+TimeAdjust(DateSeries(1,2))):find(dateQ==DateSeries(end,1)+TimeAdjust(DateSeries(end,2)));
% Series2(Sampseries) = GSeries.Data(:,2);
% % Series specific adjustment
% Samp=SampQ(find(dateQ==1948)):SampQ(find(dateQ==2023));
% % Hours per capita
% Series=log(Series1./Series2);
% NameSeries= 'Hours Worked per capita';


% Zero Padding
% padding=0;
padding=2^9;
% padding=1024;

% Smoothing window
%             = 0 : no smoothing is performed
%             = 1 : Blackman-Tukey window
%             = 2 : Parzen window
%             = 3 : Tukey-Hanning window
SpectrWindow=3;

% Estimate confidence bands
%ConfBand=0;
 ConfBand=0;

figure;
plot(dateQ(Samp),Series(Samp),'color',rgb('dimgray'),'LineWidth',3);
title(NameSeries);

Spec=Spectrum(Series(Samp),padding,ConfBand,SpectrWindow);
title(NameSeries);

%% Estimating Spectrum of AR(1)

Y_Sim(1) = 0;

for i = 2:300
    Y_Sim(i) = 0.96*Y_Sim(i-1) + 0.05*normrnd(0,1);
end

padding=2^9;


Spec=Spectrum(Y_Sim',padding,ConfBand,SpectrWindow);
title('AR(1)');



