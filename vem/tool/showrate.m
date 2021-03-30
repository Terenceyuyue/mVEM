function r = showrate(h,Err,opt1,opt2,strErr)
% 
% e.g. opt1 = 'r-*',  opt2 = 'k.-'
%      opt1: line proverty for error curve 
%      opt2: line proverty for convergence curve
%
% Copyright (C) Terence Yu.

if nargin == 4
    strErr = '||u-u_h||';
end

Err(Err == 0) = 1e-16; % Prevent the case err = 0, log(err) = -Inf.
p = polyfit(log(h(1:end)),log(Err(1:end)),1);
r = p(1);
s = 0.75*Err(1)/h(1)^r;

loglog(1./h,Err,opt1,'linewidth',2);
hold on
loglog(1./h,s*h.^r,opt2,'linewidth',1);

xlabel('log(1/h)');

h_legend = legend(strErr,['O (h^{' num2str(r,'%0.2f') '})'],'location','best');
set(h_legend,'FontSize',10);