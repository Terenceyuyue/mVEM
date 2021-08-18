function varargout = boundingbox(F)
%BOUNDINGBOX Computes bounding box of a constraint object
%
% If two outputs are requested, the numerical box is returned
% [L,U] = boundingbox(F)
%
% If three outputs are requrested, the symbolic model is appended
% [L,U,C] = boundingbox(F)
%
% If only one output is requested, only the symbolic model is reyrned
% C = boundingbox(F)

% $Id: boundingbox.m,v 1.2 2008-02-14 14:53:36 joloef Exp $

switch nargout
    case 0
        [f] = boundingbox(set(F))
    case 1
        [f] = boundingbox(set(F));
        varargout{1} = f;
    case 2
        [L,U] = boundingbox(set(F));
        varargout{1} = L;
        varargout{2} = U;
    case 3
        [L,U,f] = boundingbox(set(F));
        varargout{1} = L;
        varargout{2} = U;
        varargout{3} = f;
    otherwise
end