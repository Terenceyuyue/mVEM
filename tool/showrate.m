function r = showrate(h,Err,varargin)
% varargin = { opt1, opt2, strErr };
% e.g. opt1 = 'r-*',  opt2 = 'k.-', strErr = '||u-u_h||'
%      opt1: line proverty for error curve 
%      opt2: line proverty for convergence curve
%
% Copyright (C) Long Chen, Modified by Terence Yu.

switch length(varargin)
    case 0  % showrate(h,Err)
        opt1 = 'r-*';  opt2 = 'k.-'; strErr = '||u-u_h||'; 
    case 1  % showrate(h,Err,strErr)
        opt1 = 'r-*';  opt2 = 'k.-'; strErr = varargin{1};   
    case 2 % showrate(h,Err,opt1,opt2)
        opt1 = varargin{1}; opt2 = varargin{2}; strErr = '||u-u_h||';
    case 3 % showrate(h,Err,opt1,opt2,strErr)
        opt1 = varargin{1}; opt2 = varargin{2}; strErr = varargin{3};
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