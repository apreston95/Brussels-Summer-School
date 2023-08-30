function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
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
%   g2
%

if T_flag
    T = stnd_rbc_dyn.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
g2_i = zeros(21,1);
g2_j = zeros(21,1);
g2_v = zeros(21,1);

g2_i(1)=1;
g2_i(2)=1;
g2_i(3)=1;
g2_i(4)=1;
g2_i(5)=1;
g2_i(6)=1;
g2_i(7)=1;
g2_i(8)=1;
g2_i(9)=1;
g2_i(10)=1;
g2_i(11)=1;
g2_i(12)=1;
g2_i(13)=1;
g2_i(14)=1;
g2_i(15)=1;
g2_i(16)=1;
g2_i(17)=1;
g2_i(18)=2;
g2_i(19)=2;
g2_i(20)=2;
g2_i(21)=2;
g2_j(1)=25;
g2_j(2)=85;
g2_j(3)=81;
g2_j(4)=41;
g2_j(5)=86;
g2_j(6)=96;
g2_j(7)=87;
g2_j(8)=107;
g2_j(9)=37;
g2_j(10)=42;
g2_j(11)=92;
g2_j(12)=43;
g2_j(13)=103;
g2_j(14)=97;
g2_j(15)=98;
g2_j(16)=108;
g2_j(17)=109;
g2_j(18)=25;
g2_j(19)=1;
g2_j(20)=37;
g2_j(21)=61;
g2_v(1)=(-((-params(5))*(-params(5))*exp((-params(5))*y(3))));
g2_v(2)=(1+params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))-params(3))*params(2)*(-params(5))*(-params(5))*exp((-params(5))*y(8));
g2_v(3)=T(3)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
g2_v(4)=g2_v(3);
g2_v(5)=params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*T(3);
g2_v(6)=g2_v(5);
g2_v(7)=T(3)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g2_v(8)=g2_v(7);
g2_v(9)=T(1)*params(1)*(params(1)-1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
g2_v(10)=T(1)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
g2_v(11)=g2_v(10);
g2_v(12)=T(1)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g2_v(13)=g2_v(12);
g2_v(14)=T(1)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
g2_v(15)=T(1)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g2_v(16)=g2_v(15);
g2_v(17)=T(1)*params(1)*(-(params(1)-1))*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g2_v(18)=exp(y(3));
g2_v(19)=(-((1-params(3))*exp(y(1))));
g2_v(20)=exp(y(4));
g2_v(21)=(-exp(y(6)));
g2 = sparse(g2_i,g2_j,g2_v,5,121);
end
