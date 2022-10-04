function [s, r] = StochasticTrajectoriesTiev2(s0, q, hb, s_lower, s_upper, n_max, hb_mode, roundS)
%INPUTS:    s0  - number of initially detectable antigens
%           q   - recognition probability
%           hb  - affine penalty
%           s_lower - s lower limit to stop process
%           s_upper - s upper limit to stop process
%           n_max - upper limit on number of iterations before stopping
%           process
%           hb_mode - if 1, then hb=hb
%                     if 0, then hb=hb*rn
%           roundS  - if 1, then state variable S is rounded to nearest
%           integer at end of each period
%                     if 0, then state variable may take rational values.

%OUTPUTS: S - Number of current recognizable antigens
%         R - Number of recognied antigens         


c=-log(1-q);
s=nan(1,n_max+1);
s(1)=s0;
r=nan(1,n_max+1);
n=1;
while n<n_max+1 && s(n)<s_upper && s(n)>s_lower
%I. Simulate r
r(n)=binornd(round(s(n)),q);

%II. Create transition
if hb_mode==1
    s(n+1)=s(n)+(1-c)/c*r(n) + hb;
elseif hb_mode==0
    s(n+1)=s(n)+(1-c)/c*r(n) + hb*r(n);
else
    disp('problem with hb type!')
    break
end

if roundS==1
    s(n+1)=round(s(n+1));
end
n=n+1;

end








end
