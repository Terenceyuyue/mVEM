function str = char(expr)
%CHAR Convert a given MPTAFFEXPR object to a string representation

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

gX = expr.gX;
gU = expr.gU;
gC = expr.gC;

nx = size(gX,2);
nu = size(gU,2);
str = '';

ineqstr = '<=';

[dummy, firstnonzerou] = find(gU~=0);
isgXempty = isempty(gX) | all(gX==0);
isgUempty = isempty(gU) | all(gU==0);
if ~isgXempty,
    nonzero = gX(find(gX~=0));
    if ~isempty(nonzero),
        if nonzero(1)<0,
            gX = -gX;
            gC = -gC;
            gU = -gU;
            ineqstr = '>=';
        end
    end
elseif ~isgUempty,
    nonzero = gU(find(gU~=0));
    if ~isempty(nonzero)
        if gU(firstnonzerou(1))<0,
            gU = -gU;
            gC = -gC;
            ineqstr = '>=';
        end
    end
end

[dummy, firstnonzerox] = find(gX~=0);
for ii=1:nx
    coef = gX(ii);
    if coef==0
        continue
    elseif ii==firstnonzerox(1),
        if coef==1,
            str = [str sprintf('x(%d)', ii)];
        elseif coef==-1,
            str = [str sprintf('-x(%d)', ii)];
        elseif coef < 0,
            str = [str sprintf('-%f*x(%d)', abs(coef), ii)];
        else
            str = [str sprintf('%f*x(%d)', coef, ii)];
        end
        
    elseif coef==1
        if ii>1,
            str = [str sprintf(' + x(%d)', ii)];
        else
            str = [str sprintf('x(%d)', ii)];
        end
    elseif coef==-1
        if ii==1,
            str = [str sprintf('-x(%d)', ii)];
        else
            str = [str sprintf(' - x(%d)', ii)];
        end
    elseif coef < 0,
        if ii==1,
            str = [str sprintf('%f*x(%d)', coef, ii)];
        else
            str = [str sprintf(' - %f*x(%d)', abs(coef), ii)];
        end
    else
        if ii==1,
            str = [str sprintf('%f*x(%d)', coef, ii)];
        else
            str = [str sprintf(' + %f*x(%d)', coef, ii)];
        end
    end        
    
end

allXzero = all(gX==0);
[dummy, firstnonzerou] = find(gU~=0);

for ii=1:nu,
    coef = gU(ii);
    if coef==0
        continue
    elseif ii==firstnonzerou(1) & allXzero,
        if coef==1,
            str = [str sprintf('u(%d)', ii)];
        elseif coef==-1,
            str = [str sprintf('-u(%d)', ii)];
        elseif coef < 0,
            str = [str sprintf('-%f*u(%d)', abs(coef), ii)];
        else
            str = [str sprintf('%f*u(%d)', coef, ii)];
        end
    elseif coef==1
        if ii>1,
            str = [str sprintf(' + u(%d)', ii)];
        elseif nx>0,
            str = [str sprintf(' + u(%d)', ii)];
        else
            str = [str sprintf('u(%d)', ii)];
        end
    elseif coef==-1
        if ii==1,
            str = [str sprintf('-u(%d)', ii)];
        else
            str = [str sprintf(' - u(%d)', ii)];
        end
    elseif coef < 0,
        if ii==1,
            if nx==0,
                str = [str sprintf('%f*u(%d)', coef, ii)];
            else
                str = [str sprintf(' - %f*u(%d)', abs(coef), ii)];
            end
        else
            str = [str sprintf(' - %f*u(%d)', abs(coef), ii)];
        end
    else
        if ii==1,
            if nx==0,
                str = [str sprintf('%f*u(%d)', coef, ii)];
            else
                str = [str sprintf(' + %f*u(%d)', coef, ii)];
            end
        else
            str = [str sprintf(' + %f*u(%d)', coef, ii)];
        end
    end        
    
end    

str = [str ' ' ineqstr ' ' num2str(gC)];
