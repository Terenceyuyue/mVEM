function varargout = pnorm(varargin)
%PNORM P-Norm of SDPVAR variable with convexity knowledge
%
% PNORM is recommended if your goal is to obtain
% a convex model, since the function PNORM is implemented
% as a so called nonlinear operator. (For p/q ==1,2,inf you should use the
% overloaded norm)
%
% t = pnorm(x,p/q), p/q >= 1
%
% Note, the pnorm is implemented using cpower, which adds
% a large number of variables and constraints

% Author Johan Löfberg
% $Id: pnorm.m,v 1.1 2010-03-26 11:35:23 joloef Exp $

switch class(varargin{1})

    case 'double'
        varargout{1} = sum(varargin{1}.^varargin{2}).^(1/varargin{2});

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        X = varargin{1};
        [n,m] = size(X);
        if isreal(X) & min(n,m)==1
            if varargin{2}>=1
                varargout{1} = yalmip('define',mfilename,varargin{:});
            else
                error('PNORM only applicable for p>=1');
            end
        else
            error('PNORM can only be applied to real vectors.');
        end

    case 'char' % YALMIP send 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2}; % Second arg is the extended operator variable
            X = varargin{3}; % Third arg and above are the args user used when defining t.
            p = varargin{4};

            [p,q] = rat(p);
            absX = sdpvar(length(X),1);
            y = sdpvar(length(X),1);
            F = [-absX < X < absX];

            for i = 1:length(y)
                F = [F,pospower(absX(i),y(i),p,q)];
            end
            F = [F,pospower(t,sum(y),q,p)];

            varargout{1} = F;
            varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive','model','graph');
            varargout{3} = X;
        end
    otherwise
end

function F = pospower(x,t,p,q)
if p>q
    l = ceil(log2(abs(p)));
    r = 2^l-p;
    y = [ones(r,1)*x;ones(q,1)*t;ones(2^l-r-q,1)];
    F = detset(x,y);
else
    l = ceil(log2(abs(q)));
    y = [ones(p,1)*x;ones(2^l-q,1)*t;ones(q-p,1)];
    F = detset(t,y);
end
