function expr = plus(var1, var2)
%PLUS Sumation operator for MPTVAR objects

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


if ~(isa(var1, 'mptvar') | isa(var1, 'double'))
    error(sprintf('Variables of class "%s" not supported!', class(var1)));
end
if ~(isa(var2, 'mptvar') | isa(var2, 'double'))
    error(sprintf('Variables of class "%s" not supported!', class(var2)));
end

if isa(var1, 'mptvar') & isa(var2, 'mptvar')
    nvar1 = length(var1);
    nvar2 = length(var2);
    if nvar1 ~= 1,
        error(sprintf('"%s" must be of length 1!', char(var1)));
    end
    if nvar2 ~= 1,
        error(sprintf('"%s" must be of length 1!', char(var2)));
    end
    if nvar1 ~= nvar2
        error('Vectors must have identical number of compontents!');
    end

    norig1 = length(var1.fromvar);
    norig2 = length(var2.fromvar);

    if var1.varindex==0 
        if length(var1)>1
            error(sprintf('Variable "%s" must be referenced through an index!', char(var1)));
        else
            var1.varindex = 1;
            norig1 = 1;
        end
    end
    if var2.varindex==0,
        if length(var2)>1,
            error(sprintf('Variable "%s" must be referenced through an index!', char(var2)));
        else
            var2.varindex = 1;
            norig2 = 1;
        end
    end

    nx = 0;
    nu = 0;
    if isstate(var1) & isstate(var2)
        if norig1 ~= norig2,
            error(sprintf('Variables "%s" and "%s" must be of same length!', ...
                char(var1.fromvar), char(var2.fromvar)));
        end
        nx = norig1;
    elseif isstate(var1) & isinput(var2)
        nx = norig1;
        nu = norig2;
    elseif isstate(var2) & isinput(var1)
        nx = norig2;
        nu = norig1;
    else
        nu = norig1;
    end
    
    gX = zeros(1, nx);
    gU = zeros(1, nu);
    
    if isstate(var1)
        gX(var1.varindex) = gX(var1.varindex) + var1.multiplier;
    elseif isinput(var1)
        gU(var1.varindex) = gU(var1.varindex) + var1.multiplier;
    else
        error(sprintf('Unhandled variable type "%s"!', var1.vartype));
    end
    if isstate(var2)
        gX(var2.varindex) = gX(var2.varindex) + var2.multiplier;
    elseif isinput(var2)
        gU(var2.varindex) = gU(var2.varindex) + var2.multiplier;
    else
        error(sprintf('Unhandled variable type "%s"!', var2.vartype));
    end
    gC = 0;
    expr = mptaffexpr(nx, nu, gX, gU, gC);
    return
    
end

if isa(var1, 'double') & isa(var2, 'mptvar')
    dummy = var2;
    var2 = var1;
    var1 = dummy;
end

if isa(var1, 'mptvar')
    if ~isa(var2, 'double')
        error(sprintf('Variables of class "%s" not supported!', class(var2)));
    else
        % var1 is MPTVAR, var2 is DOUBLE
        nvar1 = length(var1);
        nvar2 = length(var2);
        if nvar1 ~= 1,
            error(sprintf('"%s" must be of length 1!', char(var1)));
        end
        if nvar2 ~= 1,
            error(sprintf('"%s" must be of length 1!', char(var2)));
        end
        if nvar1 ~= nvar2
            error('Vectors must have identical number of compontents!');
        end
        norig1 = length(var1.fromvar);
        norig2 = length(var2);
        if var1.varindex==0 
            if length(var1)>1
                error(sprintf('Variable "%s" must be referenced through an index!', char(var1)));
            else
                var1.varindex = 1;
                norig1 = 1;
            end
        end
        nx = 0;
        nu = 0;
        if isstate(var1)
            nx = norig1;
        elseif isinput(var1)
            nu = norig1;
        else
            error(sprintf('Unhandled variable type "%s"!', var1.vartype));
        end
        
        gX = zeros(1, nx);
        gU = zeros(1, nu);
        
        if isstate(var1)
            gX(var1.varindex) = gX(var1.varindex) + var1.multiplier;
        elseif isinput(var1)
            gU(var1.varindex) = gU(var1.varindex) + var1.multiplier;
        else
            error(sprintf('Unhandled variable type "%s"!', var1.vartype));
        end

        gC = -var2;
        expr = mptaffexpr(nx, nu, gX, gU, gC);
        return
        
            
        var = var1;
        var.multiplier = var2;
    end
else
    error('Unhandled case!');
end    
