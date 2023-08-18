function rpr = fnrpr(rcind,kmn,kmx,x,dtyp)

% Inputs:
%   rcind: vector of recession indicators (monthly or quarterly)
%   kmn: minimum value of k (in quarters)
%   kmx: maximum value of k (in quarters)
%   x: window around t+k (in quarters); can be a vector for multiple windows
%   dtyp: rcind type; 'q' = quarterly, 'm' = monthly
%
% Outputs:
%   rpr: recession probabilities; rpr(i,j) corresponds to k = kmn+i-1 and x = x(j)


if strcmp(dtyp,'m')
    fct = 3;
elseif strcmp(dtyp,'q')
    fct = 1;
end

T = numel(rcind);
nx = numel(x);

rpr = zeros(kmx-kmn+1,nx);

for jx = 1:nx       % for each window size
    
    fx = fct*x(jx);         % convert quarters to months if necessary
    wn = ones(2*fx+1,1);    % vector of one equal to total window size (2x+1)
    
    % t-th element of rcwn is true if date t has a recession +/- x quarters
    % around it
    rcwn = logical(conv(rcind,wn,'same'));
    
    for k = kmn:kmx         % for each value of k
        fk = fct*k;                 % convert quarters to months if necessary
        
        % take only dates t of rcind for which the entire t+k+/-x window is
        % in the data, and then convert to logical array
        rcindxk = logical(rcind(1:T-fk-fx));
        
        % take corresponding elements of rcwn, i.e, t-th element
        % corresponds to rcwn at t+k
        rcwnxk = rcwn(fk+1:end-fx);
        
        nrec = nnz(rcindxk);            % number of recessions
        nreck = nnz(rcindxk(rcwnxk));   % number of recessions in window around t+k
        
        rpr(k-kmn+1,jx) = nreck/nrec;   
        
    end
end

