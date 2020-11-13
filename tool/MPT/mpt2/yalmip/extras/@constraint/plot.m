function varargout = plot(varargin)
%PLOT  Plots the (projection of the) feasible region of a set of constraints
%
% p = plot(C,x,c,n,options)
%
% C:  Constraint object
% x:  Plot variables, optional. [At most three variables]
% c:  color, optional. [double] ([r g b] format) or char from 'rymcgbk'
% n:  #vertices, optional [double ]
% options: options structure from sdpsettings, optional
% Author Johan Löfberg
% $Id: plot.m,v 1.1 2007-02-28 16:20:33 joloef Exp $

varargin{1} = lmi(varargin{1});
plot(varargin{:});
