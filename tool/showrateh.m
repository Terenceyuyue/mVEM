function showrateh(h,varargin)
%showrateh displays convergence rates of an err sequence in the form of
% e.g. 
%  showrateh(h,ErruL2, 'r-*', '||u-u_h||')
%  showrateh(h,ErruL2,'r-*','||u-u_h||',  ErrpL2, 'b-s','||p-p_h||')
% 
% Copyright (C) Terence Yu.

optRate = {'k.-', 'k--', 'k-', 'k:.'};

%% Get Err, opt, strErr
n = (nargin-1)/3; % number of Err
Err = cell(1,n); opt = cell(1,n); strErr = cell(1,n);
s = 1;
for i = 1:n
    Err{i} = varargin{s};
    opt{i} = varargin{s+1};
    strErr{i} = varargin{s+2};
    s = s+3;
end

%% Plot convergence rates
str = cell(1,2*n);
for i = 1:n
    r = showrate(h,Err{i},opt{i},optRate{i});
    str{2*i-1} = strErr{i};
    str{2*i} = ['O (h^{' num2str(r,'%0.2f') '})'];
end

h_legend = legend(str,'location','best');
set(h_legend,'FontSize',10);