function varargout = or(varargin)
%XOR (overloaded)

% Author Johan Löfberg 
% $Id: xor.m,v 1.2 2007-05-02 12:33:54 joloef Exp $   

% Models OR using a nonlinear operator definition
varargout{1} = set(yalmip('define','lmixor',varargin{:}) == 1);
