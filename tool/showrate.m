function r = showrate(h,err,opt1,opt2,str)
% 
% e.g. opt1 = 'r-*',  opt2 = 'k.-'
% Copyright (C) Long Chen, modified by Terence Yu.

if nargin == 4
    str = '||u-u_h||';
end

err(err == 0) = 1e-16; % Prevent the case err = 0, log(err) = -Inf.
p = polyfit(log(h(1:end)),log(err(1:end)),1);
r = p(1);
s = 0.75*err(1)/h(1)^r;

loglog(1./h,err,opt1,'linewidth',2);
hold on
loglog(1./h,s*h.^r,opt2,'linewidth',1);

xlabel('log(1/h)');

h_legend = legend(str,['O (h^{' num2str(r,'%0.2f') '})'],'location','best');
set(h_legend,'FontSize',10);