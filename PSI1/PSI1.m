% EABCN 2023 - The Macroeconomics of Complemetarities Estimate spectrum
% Practice Session I, item 1 (PSI1)
% PSI1: Computing conditional probability of being in a reession

% as: austria,        start from 1969.
% au: australia,      start from 1948.
% br: brazil,         start from 1981.
% ca: canada,         start from 1948.
% fr: france,         start from 1953.
% ge: germany,        start from 1948.
% in: india,          start from 1956.
% it: italy,          start from 1956.
% jp: japan,          start from 1953.
% kr: korea,          start from 1962.
% me: mexico,         start from 1979.
% nz: new zealand,    start from 1962.
% ru: russia,         start from 1994.
% sa: south africa,   start from 1962.
% sd: sweden,         start from 1963.
% sw: switzerland,    start from 1956.
% tw: taiwan,         start from 1962.
% uk: united kingdom, start from 1951.
% us: united states,  start from 1948.
% g4: germany,france,uk,italy
% g7:  US, Japan, Germany, UK, France, Italy, Canada

clear all;
close all;
clc;


pprint = 1;


Vxname{1}  = 'recession_us';
Vxname{2}  = 'recession_uk';
Vxname{3}  = 'recession_fr';
Vxname{4}  = 'recession_ca';
Vxname{5}  = 'recession_it';
Vxname{6}  = 'recession_jp';
Vxname{7}  = 'recession_ge';
Vxname{8}  = 'recession_us_long';
Vxname{9}  = 'recession_ca_long';

Vcountry{1}='USA';
Vcountry{2}='UK';
Vcountry{3}='France';
Vcountry{4}='Canada';
Vcountry{5}='Italy';
Vcountry{6}='Japan';
Vcountry{7}='Germany';
Vcountry{8}='USA (Long Sample)';
Vcountry{9}='Canada (Long Sample)';

dddate{1}=[1948:.25:2017]';
dddate{2}=[1951:.25:2017]';
dddate{3}=[1953:.25:2017]';
dddate{4}=[1948:.25:2017]';
dddate{5}=[1956:.25:2017]';
dddate{6}=[1953:.25:2017]';
dddate{7}=[1948:.25:2017]';
dddate{8}=[1854.75:.25:2017.75]';


% Below I am coding the quarterly series of 1 (recession) and 0 (expansion) using
% NBER Datation.pdf and ECRI Datation.pdf

% Recoding US recessions (fpj)
% Using monthly datation, adjustment to quarters made by NBER
% Recession is from the quaeter after the peak to the quarter of the trough


dddate{8}=[1854.75:.25:2017.75]';

recession_us_long=zeros(size(dddate{8}));

recession_us_long(find(dddate{8}==1854.75):find(dddate{8}==1854.75))=1;
recession_us_long(find(dddate{8}==1857.5):find(dddate{8}==1858.75))=1;
recession_us_long(find(dddate{8}==1860.75):find(dddate{8}==1861.25))=1;
recession_us_long(find(dddate{8}==1865.25):find(dddate{8}==1867))=1;
recession_us_long(find(dddate{8}==1869.5):find(dddate{8}==1870.75))=1;
recession_us_long(find(dddate{8}==1873.75):find(dddate{8}==1879))=1;
recession_us_long(find(dddate{8}==1882.25):find(dddate{8}==1885.25))=1;
recession_us_long(find(dddate{8}==1887.5):find(dddate{8}==1888))=1;
recession_us_long(find(dddate{8}==1890.75):find(dddate{8}==1891.25))=1;
recession_us_long(find(dddate{8}==1893.25):find(dddate{8}==1894.25))=1;
recession_us_long(find(dddate{8}==1896):find(dddate{8}==1897.25))=1;
recession_us_long(find(dddate{8}==1899.75):find(dddate{8}==1900.75))=1;
recession_us_long(find(dddate{8}==1903):find(dddate{8}==1904.5))=1;
recession_us_long(find(dddate{8}==1907.5):find(dddate{8}==1908.25))=1;
recession_us_long(find(dddate{8}==1910.25):find(dddate{8}==1912.75))=1;
recession_us_long(find(dddate{8}==1913.25):find(dddate{8}==1914.75))=1;
recession_us_long(find(dddate{8}==1918.75):find(dddate{8}==1919))=1;
recession_us_long(find(dddate{8}==1920.25):find(dddate{8}==1921.5))=1;
recession_us_long(find(dddate{8}==1923.5):find(dddate{8}==1924.5))=1;
recession_us_long(find(dddate{8}==1926.75):find(dddate{8}==1927.75))=1;
recession_us_long(find(dddate{8}==1929.75):find(dddate{8}==1933))=1;
recession_us_long(find(dddate{8}==1937.5):find(dddate{8}==1938.25))=1;
recession_us_long(find(dddate{8}==1945.25):find(dddate{8}==1945.75))=1;


recession_us_long(find(dddate{8}==1949):find(dddate{8}==1949.75))=1;
recession_us_long(find(dddate{8}==1953.5):find(dddate{8}==1954.25))=1;
recession_us_long(find(dddate{8}==1957.75):find(dddate{8}==1958.25))=1;
recession_us_long(find(dddate{8}==1960.5):find(dddate{8}==1961))=1;
recession_us_long(find(dddate{8}==1970):find(dddate{8}==1970.75))=1;
recession_us_long(find(dddate{8}==1974):find(dddate{8}==1975))=1;
recession_us_long(find(dddate{8}==1980.25):find(dddate{8}==1980.5))=1;
recession_us_long(find(dddate{8}==1981.75):find(dddate{8}==1982.75))=1;
recession_us_long(find(dddate{8}==1990.75):find(dddate{8}==1991))=1;
recession_us_long(find(dddate{8}==2001.25):find(dddate{8}==2001.75))=1;
recession_us_long(find(dddate{8}==2008):find(dddate{8}==2009.25))=1;

recession_us = recession_us_long(find(dddate{8}==1948):find(dddate{8}==2017));


% Recoding UK recessions (fpj)
% Using monthly datation
% if the peak is in the first month of the quarter, that quarter is in the recession
% if the trough is not in the first month of the quarter, that quarter is not in the recession

recession_uk=zeros(size(dddate{2}));
recession_uk(find(dddate{2}==1951):find(dddate{2}==1952.25))=1;
recession_uk(find(dddate{2}==1974.75):find(dddate{2}==1975.25))=1;
recession_uk(find(dddate{2}==1979.5):find(dddate{2}==1981.25))=1;
recession_uk(find(dddate{2}==1990.5):find(dddate{2}==1992))=1;
recession_uk(find(dddate{2}==2008.5):find(dddate{2}==2009.75))=1;


% Recoding France recessions (fpj)
% Using monthly datation
% if the peak is in the first month of the quarter, that quarter is in the recession
% if the trough is not in the first month of the quarter, that quarter is not in the recession

recession_fr=zeros(size(dddate{3}));
recession_fr(find(dddate{3}==1958):find(dddate{3}==1959))=1;
recession_fr(find(dddate{3}==1974.25):find(dddate{3}==1975.25))=1;
recession_fr(find(dddate{3}==1979.75):find(dddate{3}==1980.25))=1;
recession_fr(find(dddate{3}==1982.25):find(dddate{3}==1984.75))=1;
recession_fr(find(dddate{3}==1992.25):find(dddate{3}==1993.5))=1;
recession_fr(find(dddate{3}==2002.25):find(dddate{3}==2003.25))=1;
recession_fr(find(dddate{3}==2008.25):find(dddate{3}==2009))=1;
recession_fr(find(dddate{3}==2011.25):find(dddate{3}==2012.75))=1;

% Recoding Italy recessions (fpj)
% Using monthly datation
% if the peak is in the first month of the quarter, that quarter is in the recession
% if the trough is not in the first month of the quarter, that quarter is not in the recession

recession_it=zeros(size(dddate{5}));
recession_it(find(dddate{5}==1964):find(dddate{5}==1965.25))=1;
recession_it(find(dddate{5}==1970.75):find(dddate{5}==1971.5))=1;
recession_it(find(dddate{5}==1974.25):find(dddate{5}==1975))=1;
recession_it(find(dddate{5}==1980.25):find(dddate{5}==1983.25))=1;
recession_it(find(dddate{5}==1992.25):find(dddate{5}==1993.5))=1;
recession_it(find(dddate{5}==2008.5):find(dddate{5}==2009))=1;
recession_it(find(dddate{5}==2011.25):find(dddate{5}==2014.5))=1;

% Recoding Japan recessions (fpj)
% Using monthly datation
% if the peak is in the first month of the quarter, that quarter is in the recession
% if the trough is not in the first month of the quarter, that quarter is not in the recession

recession_jp=zeros(size(dddate{6}));
recession_jp(find(dddate{6}==1953):find(dddate{6}==1954.75))=1;
recession_jp(find(dddate{6}==1973.75):find(dddate{6}==1975))=1;
recession_jp(find(dddate{6}==1992.25):find(dddate{6}==1994))=1;
recession_jp(find(dddate{6}==1997.25):find(dddate{6}==1999.25))=1;
recession_jp(find(dddate{6}==2000.75):find(dddate{6}==2003))=1;
recession_jp(find(dddate{6}==2008.25):find(dddate{6}==2009))=1;
recession_jp(find(dddate{6}==2010.75):find(dddate{6}==2011))=1;
recession_jp(find(dddate{6}==2012.5):find(dddate{6}==2012.75))=1;
recession_jp(find(dddate{6}==2014.25):find(dddate{6}==2014.5))=1;

% Recoding German recessions (fpj)
% Using monthly datation
% if the peak is in the first month of the quarter, that quarter is in the recession
% if the trough is not in the first month of the quarter, that quarter is not in the recession

recession_ge=zeros(size(dddate{7}));
recession_ge(find(dddate{7}==1966.25):find(dddate{7}==1967.25))=1;
recession_ge(find(dddate{7}==1973.75):find(dddate{7}==1975.25))=1;
recession_ge(find(dddate{7}==1980):find(dddate{7}==1982.5))=1;
recession_ge(find(dddate{7}==1991):find(dddate{7}==1994))=1;
recession_ge(find(dddate{7}==1997.25):find(dddate{7}==1999.25))=1;
recession_ge(find(dddate{7}==2001):find(dddate{7}==2003.5))=1;
recession_ge(find(dddate{7}==2008.25):find(dddate{7}==2009))=1;


% Recoding Canada recessions (fpj)
% Using monthly datation, adjustment to quarters made by CD Howe Institute

dddate{9}=[1926:.25:2017]';
recession_ca_long=zeros(size(dddate{9}));
recession_ca_long(find(dddate{9}==1929.25):find(dddate{9}==1933))=1;
recession_ca_long(find(dddate{9}==1937.5):find(dddate{9}==1938.25))=1;
recession_ca_long(find(dddate{9}==1947.25):find(dddate{9}==1948))=1;
recession_ca_long(find(dddate{9}==1951):find(dddate{9}==1951.75))=1;
recession_ca_long(find(dddate{9}==1953.25):find(dddate{9}==1954.25))=1;
recession_ca_long(find(dddate{9}==1957):find(dddate{9}==1958))=1;
recession_ca_long(find(dddate{9}==1960):find(dddate{9}==1961))=1;
recession_ca_long(find(dddate{9}==1974.75):find(dddate{9}==1975))=1;
recession_ca_long(find(dddate{9}==1979.75):find(dddate{9}==1980.25))=1;
recession_ca_long(find(dddate{9}==1981.25):find(dddate{9}==1982.75))=1;
recession_ca_long(find(dddate{9}==1990):find(dddate{9}==1992.25))=1;
recession_ca_long(find(dddate{9}==2008.5):find(dddate{9}==2009.25))=1;


recession_ca = recession_ca_long(find(dddate{9}==1948):find(dddate{9}==2017));





%% Computing conditional correlations

for in=[1:9]


    eval(['rpr = fnrpr(' deblank(Vxname{in}) ',12,90,[1 2 3 4 5],''q'');']);

    figure
    legendInfo = cell(3,1);
    hold on
    j=3;   plot(12:90, rpr(:,j),'-','linewidth',4,'color',rgb('lightgray'))
    legendInfo{j-2} = strcat('$x =\,$ ' , num2str(j), ' quarters');
    j=4;    plot(12:90, rpr(:,j),'-.','linewidth',4,'color',rgb('dimgray'))
    legendInfo{j-2} = strcat('$x =\,$ ' , num2str(j), ' quarters');
    j=5;    plot(12:90, rpr(:,j),'linewidth',4,'color',rgb('dimgray'))
    legendInfo{j-2} = strcat('$x =\,$ ' , num2str(j), ' quarters');
    xlabel('$k$ (Quarters Since Recession)','interpreter','Latex','FontSize',16);
    ylabel('Probability of Recession','interpreter','Latex','FontSize',16);
    ll=legend(legendInfo, 'location','best');
    set(ll,'FontName','Times','FontSize',15,'interpreter', 'latex')
    set(ll,'color','none');
    title(Vcountry{in}); 
    hold off
    if pprint==1,
        eval(['print -dpng ' Vxname{in}]);
    end

% 
%         figure
%         hold on
%         j=5;    plot(12:90, rpr(:,j),'linewidth',4,'color',rgb('dimgray'))
%         xlabel('$k$ (Quarters Since Recession)','interpreter','Latex','FontSize',16);
%         ylabel('Probability of Recession','interpreter','Latex','FontSize',16);
%         set(ll,'FontName','Times','FontSize',15,'interpreter', 'latex')
%         set(ll,'color','none');
%         hold off
%         if pprint==1,
%             pause
%             eval(['print -depsc2 ' Vxname{in}]);
%         end
%         title(Vcountry{in});



end

