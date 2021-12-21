function showrateErr(h,varargin)
%showrateErr displays convergence rates of an err sequence with given names
% such as ErrL2, ErrH1, ErrH2, ErrI
%  e.g. showrateErr(h,ErrL2,ErrH1,ErrH2,ErrI)
%
% Copyright (C) Terence Yu.

%% Determine number of Err
n = 0;
for i = 2:nargin
    varName = lower(inputname(i));
    if mycontains(varName, 'err'),  n = n + 1; end
end

%% Plot convergence rates
str = cell(1,2*n);  
for i = 1:n
    var = upper(inputname(i+1));
    if mycontains(var, 'L2') 
        stri = '||u-u_h||';  
        r = showrate(h,varargin{i},'r-*','k.-');
    end
    if mycontains(var, 'H1')
        stri = '|u-h_h|_1';  
        r = showrate(h,varargin{i},'b-s','k--');
    end
    if mycontains(var, 'H2')
        stri = '|u-u_h|_2'; 
        r = showrate(h,varargin{i},'b-o','r.-');
    end
    if mycontains(var, 'I')
        stri = '||u_I-u_h||_E'; 
        r = showrate(h,varargin{i},'k-d','k:.');        
    end
    str{2*i-1} = stri;
    str{2*i} = ['O (h^{' num2str(r,'%0.2f') '})'];
end

h_legend = legend(str,'location','best');
set(h_legend,'FontSize',10);