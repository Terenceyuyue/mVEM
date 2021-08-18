function F = constraint(X,quantifier,Y)
% Internal class for constraint list

% Author Johan Löfberg
% $Id: constraint.m,v 1.8 2009-04-29 07:48:12 joloef Exp $

superiorto('sdpvar');
superiorto('double');

if isa(X,'blkvar')
    X = sdpvar(X);
end
if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

% Evaluate and save a SET object
switch quantifier
    case '>'
        Z = X - Y;
       % C = set(X > 0);
    case '>='
        Z = X - Y;
       % C = set(X >= 0);
    case '<'
        Z = Y - X;
       % C = set(Z > 0);
    case '<='
        Z = Y - X;
       % C = set(Z >= 0);
    case '=='
        Z = Y - X;
       % C = set(Z == 0);
        %case {'>','>='}
        %    Z = X - Y;
        %case {'<','<=','=='}
        %    Z = Y - X;
    otherwise
        error('Quantifier not supported')
end

if isequal(Z,0)
    warning('Constraint evaluated to trivial true.')
    F = set([]);
    return
end

switch quantifier
case {'>','<'}
    F.strict(1) = 1;
case {'>=','<=','=='}
    F.strict(1) = 0;
otherwise
    error('Quantifier not supported')
end

F.List={X,quantifier,Y};
F.Evaluated{1} = Z;
F.ConstraintID = yalmip('ConstraintID');
F.tag{1} = '';
F = class(F,'constraint');
	