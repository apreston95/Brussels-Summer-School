% M2Qf.m
%
% Purpose: Transform flow monthly series into a quarterly one
%

function yQ = M2Qf(yM)

N = length(yM);
k = 0;
for i = 3:3:N
	k = k + 1;
	yQ(k) = yM(i) + yM(i-1) + yM(i-2);
end
