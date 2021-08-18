function varargout = log2(varargin)
%log2 (overloaded)

% Author Johan Löfberg
% $Id: log2.m,v 1.7 2007-08-02 18:16:26 joloef Exp $
switch class(varargin{1})

    case 'double' 
        error('Overloaded SDPVAR/NORM CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'      

        X = varargin{3};
        F = F + set(X > eps);

        operator = struct('convexity','concave','monotonicity','increasing','definiteness','none','model','callback');
        operator.convexhull = @convexhull;
        operator.bounds = @bounds;
        operator.derivative = @derivative;

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/LOG2 called with CHAR argument?');
end

function df = derivative(x)
df = (1./x)/log(2);

function [L,U] = bounds(xL,xU)
if xL < 0
    % The variable is not bounded enough yet
    L = -inf;
elseif xL==0
    L = -inf;
else
    L = log2(xL);
end
if xU < 0
    % This is an infeasible problem
    L = inf;
    U = -inf;
else
    U = log2(xU);
end

function [Ax, Ay, b] = convexhull(xL,xU)
fL = log2(xL);
fU = log2(xU);
dfL = (1/(xL))/log(2);
dfU = (1/(xU))/log(2);
[Ax,Ay,b] = convexhullConcave(xL,xU,fL,fU,dfL,dfU);