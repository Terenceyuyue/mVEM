function x = uncertain(x)
%UNCERTAIN Declares a variable as uncertain
%
%   F = UNCERTAIN(w) is used to describe the set of uncertain variables
%   in an uncertain program
%
%   INPUT
%    w : SDPVAR object
%
%   OUTPUT
%    F : Constraint object
%
%   EXAMPLE
%    sdpvar x w
%    F = [x + w <= 1, -0.5 <= w <= 0.5, uncertain(w)];
%    solvesdp(F,-x) 
%
%   See also SOLVESDP, ROBUSTIFY

% Author Johan Löfberg
% $Id: uncertain.m,v 1.3 2006-08-18 15:01:04 joloef Exp $

x.typeflag = 15;
x = lmi(x);