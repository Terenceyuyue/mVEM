function ctrl = copyvalidfields(ctrlStruct, allowedfields)
%COPYVALIDFIELDS copies only fields listed in 'allowedfields' from ctrlStruct to ctrl

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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

ctrl = ctrlStruct;

f = fields(ctrlStruct);
lf = length(f);
field_allowed = zeros(lf,1); 
for ii=1:lf;
    for jj=1:length(allowedfields),
        if strcmp(f{ii}, allowedfields{jj}),
            field_allowed(ii) = 1;
            break
        end
    end
end
if any(field_allowed==0)
    ctrl.details.additional = [];
end
for ii=1:lf,
    if field_allowed(ii)==0
        notallowed = getfield(ctrlStruct, f{ii});
        ctrl.details.additional = setfield(ctrl.details.additional, f{ii}, notallowed);
        ctrl = rmfield(ctrl, f{ii});
    end
end
