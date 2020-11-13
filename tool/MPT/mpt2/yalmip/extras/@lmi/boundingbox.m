function varargout = boundingbox(F,ops)
%BOUNDINGBOX Computes bounding box of a constraint object
%
% If two outputs are requested, the numerical bounds are returned
% [L,U] = boundingbox(F)
%
% If three outputs are requrested, the symbolic model is appended
% [L,U,B] = boundingbox(F)
%
% If only one output is requested, only the symbolic model is returned
% B = boundingbox(F)



% Author Johan Löfberg
% $Id: boundingbox.m,v 1.1 2004-12-08 00:07:15 johanl Exp $

x = recover(depends(F));

if nargin < 2
    ops = sdpsettings('verbose',0);
end

for i = 1:length(x);
    sol = solvesdp(F,x(i),ops);
    L(i,1) = double(x(i));
    sol = solvesdp(F,-x(i),ops);
    U(i,1) = double(x(i));
end

switch nargout
    case 0
        [L < x < U]
    case 1
        varargout{1} = [L < x < U];
    case 2
        varargout{1} = L;
        varargout{2} = U;
    case 3
        varargout{1} = L;
        varargout{2} = U;
        varargout{3} = [L < x < U];
    otherwise
end