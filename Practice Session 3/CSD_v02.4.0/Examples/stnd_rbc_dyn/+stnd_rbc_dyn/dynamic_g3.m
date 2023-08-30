function g3 = dynamic_g3(T, y, x, params, steady_state, it_, T_flag)
% function g3 = dynamic_g3(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g3
%

if T_flag
    T = stnd_rbc_dyn.dynamic_g3_tt(T, y, x, params, steady_state, it_);
end
g3_i = zeros(25,1);
g3_j = zeros(25,1);
g3_v = zeros(25,1);

g3_i(1)=1;
g3_i(2)=1;
g3_i(3)=1;
g3_i(4)=1;
g3_i(5)=1;
g3_i(6)=1;
g3_i(7)=1;
g3_i(8)=1;
g3_i(9)=1;
g3_i(10)=1;
g3_i(11)=1;
g3_i(12)=1;
g3_i(13)=1;
g3_i(14)=1;
g3_i(15)=1;
g3_i(16)=1;
g3_i(17)=1;
g3_i(18)=1;
g3_i(19)=1;
g3_i(20)=1;
g3_i(21)=1;
g3_i(22)=2;
g3_i(23)=2;
g3_i(24)=2;
g3_i(25)=2;
g3_j(1)=267;
g3_j(2)=932;
g3_j(3)=928;
g3_j(4)=933;
g3_j(5)=934;
g3_j(6)=884;
g3_j(7)=889;
g3_j(8)=890;
g3_j(9)=944;
g3_j(10)=945;
g3_j(11)=956;
g3_j(12)=400;
g3_j(13)=405;
g3_j(14)=406;
g3_j(15)=460;
g3_j(16)=461;
g3_j(17)=472;
g3_j(18)=1065;
g3_j(19)=1066;
g3_j(20)=1077;
g3_j(21)=1198;
g3_j(22)=267;
g3_j(23)=1;
g3_j(24)=400;
g3_j(25)=666;
g3_v(1)=(-((-params(5))*(-params(5))*(-params(5))*exp((-params(5))*y(3))));
g3_v(2)=(1+params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))-params(3))*params(2)*(-params(5))*(-params(5))*(-params(5))*exp((-params(5))*y(8));
g3_v(3)=params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*params(2)*(-params(5))*(-params(5))*exp((-params(5))*y(8));
g3_v(4)=params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*params(2)*(-params(5))*(-params(5))*exp((-params(5))*y(8));
g3_v(5)=params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1))*params(2)*(-params(5))*(-params(5))*exp((-params(5))*y(8));
g3_v(6)=T(3)*params(1)*(params(1)-1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
g3_v(7)=T(3)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
g3_v(8)=T(3)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g3_v(9)=params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*T(3);
g3_v(10)=T(3)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g3_v(11)=T(3)*params(1)*(-(params(1)-1))*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g3_v(12)=T(1)*params(1)*(params(1)-1)*(params(1)-1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
g3_v(13)=T(1)*params(1)*(params(1)-1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
g3_v(14)=T(1)*params(1)*(params(1)-1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g3_v(15)=T(1)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
g3_v(16)=T(1)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g3_v(17)=T(1)*params(1)*(params(1)-1)*(-(params(1)-1))*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g3_v(18)=T(1)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
g3_v(19)=T(1)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g3_v(20)=T(1)*params(1)*(-(params(1)-1))*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g3_v(21)=T(1)*params(1)*(-(params(1)-1))*(-(params(1)-1))*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g3_v(22)=exp(y(3));
g3_v(23)=(-((1-params(3))*exp(y(1))));
g3_v(24)=exp(y(4));
g3_v(25)=(-exp(y(6)));
g3 = sparse(g3_i,g3_j,g3_v,5,1331);
end
