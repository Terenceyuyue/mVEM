function F = set(varargin)
%set               Defines a constraint (the feasible set)
%   
%    F = SET               Creates an empty SET-object
%
%  Constraints can be generated using string notation
%    F = SET('X>Y')        Constrains X-Y to be positive semi-definite if X-Y is Hermitian,
%                          interpreted as element-wise constraint otherwise
%    F = SET('X==Y')       Element-wise equality constraint  
%    F = SET('||X||<Y')    Create second order cone constraint (X and Y column vectors)
%
%  One can also use overloaded >, < and ==
%    F = SET(X > Y)          Constrains X-Y to be positive semi-definite if X-Y is Hermitian,
%                          interpreted as element-wise constraint otherwise
%    F = SET(X==Y)         Element-wise equality constraint
%    F = SET(CONE(X,Y))    Second order cone constraint (X and Y column vectors)
%
%  Variables can be constrained to be integer or binary
%    F = SET(INTEGER(X))
%    F = SET(BINARY(X))
%
%  Multiple constraints are obtained with overloaded plus
%    F = set(X > 0) + set(CONE(X(:),1)) + SET(X(1,1) == 1/2)
%
%  Double-sided constraint (and extensions) can easily be defined
%  The following two comands give equivalent problems
%    F = set(X > 0 > Y > Z < 5 < W)
%    F = set(X > 0) + set(0 > Y) + set(Y > Z) + set(Z < 5) + set(5 < W)
%
%
%  General info
%    A constraint can be tagged with a name or description 
%    F = SET(X > Y,'tag')  Gives the constraint a tag (used in display/checkset)  
%
%    The right-hand side and left-hand side can be interchanged. Supports {>,<,==}.
%
%    All inequalities are interpreted as non-strict.
%
%    For notational purposes though, both >= and > are supported (as well as < and <=)
%
%    Any valid expression built using DOUBLEs & SDPVARs can be used on both sides.
%
%    The advantage of using the string notation approach is that more information will be  
%    shown when the SET is displayed (and in checkset)
%
%    See also   DUAL, SOLVESDP, INTEGER, BINARY

switch nargin
case 0
    F = lmi;
case 1
    F = lmi(varargin{1});
case 2
    F = lmi(varargin{1},varargin{2});
case 3
    F = lmi(varargin{1},varargin{1},varargin{3});
case 4
    F = lmi(varargin{1},[],[],1);
    
otherwise
end
    