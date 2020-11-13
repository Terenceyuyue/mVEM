function display(var)
%DISPLAY Display function for MPTVAR objects

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

if var.varlength>1,
    fprintf('\n%s vector (%dx1)\n\n', var.vartype, var.varlength);
else
    if var.varindex == 0 | 1
        if var.multiplier==0,
            fprintf('\n0\n\n');
        elseif var.multiplier==1
            fprintf('\n%s\n\n', var.varname);
        else
            fprintf('\n%f*%s\n\n', var.multiplier, var.varname);
        end
    else
        if var.multiplier==0,
            fprintf('\n0\n\n');
        elseif var.multiplier==1
            fprintf('\n%s(%d)\n\n', var.varname, var.varindex);
        else
            fprintf('\n%f*%s(%d)\n\n', var.multiplier, var.varname, var.varindex);
        end
    end
end
