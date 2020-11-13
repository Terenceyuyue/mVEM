function expr = plus(expr1, expr2)
%PLUS Sumation operator for MPTAFFEXPR objects

% Copyright is with the following author(s):
%
%(C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%         kvasnica@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

superiorto('mptvar');

if isa(expr2, 'mptaffexpr') & ~isa(expr1, 'mptaffexpr')
    dummy = expr2;
    expr2 = expr1;
    expr1 = dummy;
end

expr = expr1;
if isa(expr1, 'mptaffexpr') & isa(expr2, 'mptaffexpr')
    
    gX1 = expr1.gX;
    gU1 = expr1.gU;
    gC1 = expr1.gC;
    gX2 = expr2.gX;
    gU2 = expr2.gU;
    gC2 = expr2.gC;
    
    nx1 = size(gX1, 2);
    nx2 = size(gX2, 2);
    nu1 = size(gU1, 2);
    nu2 = size(gU2, 2);

    if nx1>0 & nx2>0,
        if nx1 ~= nx2,
            error('Number of states in an affine expression must remain constant!');
        end
    end
    if nu1>0 & nu2>0,
        if nu1 ~= nu2,
            error('Number of inputs in an affine expression must remain constant!');
        end
    end
    
    gX = zeros(1, max(nx1, nx2));
    gU = zeros(1, max(nu1, nu2));
    gC = 0;
    
    for ii = 1:nx1,
        gX(ii) = gX1(ii);
    end
    for ii = 1:nu1,
        gU(ii) = gU1(ii);
    end
    for ii = 1:nx2,
        gX(ii) = gX(ii) + gX2(ii);
    end
    for ii = 1:nu2,
        gU(ii) = gU(ii) + gU2(ii);
    end
        
    gC = gC1 + gC2;
    
    expr.gX = gX;
    expr.gU = gU;
    expr.gC = gC;
    
elseif isa(expr1, 'mptaffexpr') & isa(expr2, 'double')
    expr.gC = expr.gC - expr2;
    return
    
elseif isa(expr1, 'mptaffexpr') & isa(expr2, 'mptvar')
    
    var = expr2;
    
    gX = expr.gX;
    gU = expr.gU;
    gC = expr.gC;
    
    ind = get(var, 'varindex');
    if ind==0
        if length(var)>1
            error(sprintf('Variable "%s" must be referenced through an index!', char(var)));
        end
        ind = 1;
    end
    
    if isstate(var)
        if size(gX, 2) < ind,
            gX = [gX zeros(1, ind - size(gX,2))];
        end
        gX(ind) = gX(ind) + get(var, 'multiplier');
    elseif isinput(var)
        if size(gU, 2) < ind,
            gU = [gU zeros(1, ind - size(gU,2))];
        end
        gU(ind) = gU(ind) + get(var, 'multiplier');
    else
        error(sprintf('Unhandled MPTVAR type "%s"!', get(var, vartype)));
    end
    
    expr.gX = gX;
    expr.gU = gU;
    expr.gC = gC;
   
    return
    
else
    error(sprintf('Cannot sum objects of classes "%s" and "%s"!', class(expr1), class(expr2)));
end
