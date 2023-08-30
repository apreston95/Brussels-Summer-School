% PLOT_NBER_2SERIES plots 2 series on the same graph, with shaded bars for
% NBER US recessions
%
%( Franck Portier  dec 2012, adapted from Fabrice Collard plot_nber)
%
% plot_nber_2series(x,y,Tstart,Tend,legx,legy,LLloc)
%
% x and y are the two series to plot
%
% Tstart and Tend are first and last period
% 1970    == 1970Q1
% 1970.25 == 1970Q2
% 1970.5  == 1970Q3
% 1970.75 == 1970Q4
%
% legx and legy are the names of the werties for the legend (should be
% strings)
%
% LLloc is the location of the legend (should be a string in the following
% list:
% North               Inside plot box near top
% South               Inside bottom
% East                Inside right
% West                Inside left
% NorthEast           Inside top right (default for 2-D plots)
% NorthWest           Inside top left
% SouthEast           Inside bottom right
% SouthWest           Inside bottom left
% NorthOutside        Outside plot box near top
% SouthOutside        Outside bottom
% EastOutside         Outside right
% WestOutside         Outside left
% NorthEastOutside    Outside top right (default for 3-D plots)
% NorthWestOutside    Outside top left
% SouthEastOutside    Outside bottom right
% SouthWestOutside    Outside bottom left
% Best                Least conflict with data in plot
% BestOutside         Least unused space outside plot
 


function h=plot_nber_2series(x,y,Tstart,Tend,legx,legy,LLloc)


Peak    = [ 1857.25 1860.50 1865.00 1869.25 1873.50 1882.00 1887.25 1890.50 1893.00 1895.75 1899.50 1902.75 1907.25 1910.00 ...
            1913.00 1918.50 1920.00 1923.25 1926.50 1929.50 1937.25 1945.00 1948.75 1953.25 1957.50 1960.25 1969.75 1973.75 ...
            1980.00 1981.50 1990.50 2001.00 2007.75 2019.75
        ];

Trough  = [ 1858.75 1861.50 1867.00 1870.75 1879.00 1885.25 1888.00 1891.25 1894.25 1897.25 1900.75 1904.50 1908.25 1912.75 ...
            1914.75 1919.00 1921.50 1924.50 1927.75 1933.00 1938.25 1945.75 1949.75 1954.25 1958.25 1961.00 1970.75 1975.00 ...
            1980.50 1982.75 1991.00 2001.75 2009.25 2020.25
        ];

i1      = find(Peak>=Tstart);
i1      = i1(1);
i2      = find(Trough<=Tend);
i2      = i2(end);
Peak    = Peak(i1:i2);
Trough  = Trough(i1:i2);
time    = Tstart:0.25:Tend;
COL     = [1 1 1]*0.8;
%plot(time,x,'-b',time,y,'--r','linewidth',2);
h1=plot(time,x,'-k');
hold on;
h=plot(time,y,'--','Color',rgb('DimGray'));
set(h,'linewidth',2);
set(h1,'linewidth',2);
%legend(legx,legy);
LL=legend(legx,legy);
set(LL,'AutoUpdate','off');
set(LL,'fontname','times','FontSize',16,'interpreter', 'latex','location',LLloc);
set(LL,'color','none');
ylim    = get(gca,'ylim');
for i=1:size(Peak,2);
    patch([Peak(i) Trough(i) Trough(i) Peak(i)],[ylim(1) ylim(1) ylim(2) ylim(2)],[0 0 0 0],COL,'edgecolor','none');
end
hold on;
%h       = plot(time,x,'-b',time,y,'--r');
h1       = plot(time,x,'-k');
hold on;
h       = plot(time,y,'--','Color',rgb('DimGray'));
%LL=legend(legx,legy);
set(h,'linewidth',3);
set(h1,'linewidth',3);
axis([Tstart Tend ylim(1) ylim(2)]);
set(gca,'Fontname','Times','Fontsize',16,'layer','top');
box on;


