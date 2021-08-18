function [out,dyn] = mpt_setSysStruct(dynamics, field, value)
%MPT_SETSYSSTRUCT MPT Studio helper function

% $Revision: 1.1 $ $Date: 2005/02/23 12:43:34 $
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

global mpt___sysStruct

dyn = 0;
if nargout==1 & nargin==0,
    out = mpt___sysStruct;
    return
end

if nargin==1 & isstruct(dynamics),
    mpt___sysStruct = dynamics;
    return
end

if nargin==1 & ischar(dynamics)
    out = [];
    if strcmpi(dynamics, 'validate')
        out = '';
        [res, msg, dyn] = sub_validateSysStruct(mpt___sysStruct);
        out = msg;
        return
        
    elseif strcmpi(dynamics, 'nu')
        if isfield(mpt___sysStruct, 'B'),
            if iscell(mpt___sysStruct.B),
                for ii=1:length(mpt___sysStruct.B),
                    if ~isempty(mpt___sysStruct.B{ii}),
                        out = size(mpt___sysStruct.B{1}, 2);
                        return
                    end
                end
            elseif ~isempty(mpt___sysStruct.B)
                out = size(mpt___sysStruct.B, 2);
                return
            end
            if isfield(mpt___sysStruct, 'umax')
                out = length(mpt___sysStruct.umax);
                return
            elseif isfield(mpt___sysStruct, 'umin')
                out = length(mpt___sysStruct.umin);
            end
        end
        
    elseif strcmpi(dynamics, 'nx')
        if isfield(mpt___sysStruct, 'A'),
            if iscell(mpt___sysStruct.A),
                for ii=1:length(mpt___sysStruct.A),
                    if ~isempty(mpt___sysStruct.A{ii}),
                        out = size(mpt___sysStruct.B{1}, 1);
                        return
                    end
                end
            elseif ~isempty(mpt___sysStruct.A),
                out = size(mpt___sysStruct.A, 1);
                return
            end
            if isfield(mpt___sysStruct, 'xmax')
                out = length(mpt___sysStruct.xmax);
                return
            elseif isfield(mpt___sysStruct, 'xmin')
                out = length(mpt___sysStruct.xmin);
            end
            
        end
        
    elseif strcmpi(dynamics, 'ny')
        if isfield(mpt___sysStruct, 'C'),
            if iscell(mpt___sysStruct.C),
                for ii=1:length(mpt___sysStruct.C),
                    if ~isempty(mpt___sysStruct.C{ii}),
                        out = size(mpt___sysStruct.C{1}, 1);
                        return
                    end
                end

            elseif ~isempty(mpt___sysStruct.C)
                out = size(mpt___sysStruct.C, 1);
                return
            end
            if isfield(mpt___sysStruct, 'ymax')
                out = length(mpt___sysStruct.ymax);
                return
            elseif isfield(mpt___sysStruct, 'ymin')
                out = length(mpt___sysStruct.ymin);
            end

        end

    elseif strcmpi(dynamics, 'ndyn')
        out = 1;
        if isfield(mpt___sysStruct, 'A')
            if iscell(mpt___sysStruct.A)
                out = length(mpt___sysStruct.A);
            end
        end
        return
    end
    out = [];
    return
end

if isempty(value)
    return
end
ffield = field;
convert_to_cell = 0;
field_is_cell = 0;
if strcmp(field,'A') | strcmp(field,'B') | strcmp(field,'f') | ...
        strcmp(field,'C') | strcmp(field,'D') | strcmp(field,'g') | ...
        strcmp(field,'guardX') | strcmp(field,'guardU') | strcmp(field,'guardC') | ...
        strcmp(field, 'InputName') | strcmp(field, 'OutputName') | strcmp(field, 'StateName') | ...
        strcmp(field, 'Uset'),
    
    if isfield(mpt___sysStruct, field),
        val = getfield(mpt___sysStruct, field);
        field_is_cell = iscell(val);
    end
    if field_is_cell,
        field_str = sprintf('%s{%d}', field, dynamics);
    elseif dynamics~=1,
        field_str = sprintf('%s{%d}', field, dynamics);
        convert_to_cell = 1;
    else
        field_str = sprintf('%s', field);
    end
else
    field_str = field;
end

try
    if ~isfield(mpt___sysStruct, ffield),
        curvalue = [];
    else
        curvalue = getfield(mpt___sysStruct, ffield);
    end
    if ~ischar(value),
        value = mat2str(value);
    end
    if isempty(value),
        value = [];
    end
    if ischar(value),
        if strcmpi(field, 'StateName') | strcmpi(field, 'InputName') | strcmpi(field, 'OutputName'),
            evalvalue = value;
        else
            eval(['evalvalue = ' value ';']);
        end
    else
        eval(['evalvalue = ' value ';']);
    end
    newvalue = curvalue;
    if ~iscell(curvalue) & convert_to_cell,
        newvalue = {curvalue};
        newvalue{dynamics} = evalvalue;
    end
    if ~convert_to_cell,
        if iscell(curvalue),
            newvalue{dynamics} = evalvalue;
        else
            newvalue = evalvalue;
        end
    end
    
    mpt___sysStruct = setfield(mpt___sysStruct, ffield, newvalue);
    
catch
    errordlg(lasterr);
end



%================================================================
function [res, msg, dyn] = sub_validateSysStruct(sysStruct)

msg = '';
dyn = 0;
try
    Options.verbose = 0;
    Options.guierrors = 1;
    T = evalc('sysStruct = mpt_verifySysStruct(sysStruct, Options);');
    if iscell(sysStruct.A),
        for ii=1:length(sysStruct.A),
            if isempty(sysStruct.A),
                res = 0;
                dyn = ii;
                msg = '"A" must not be empty!';
            elseif isempty(sysStruct.B),
                res = 0;
                dyn = ii;
                msg = '"B" must not be empty!';
            elseif isempty(sysStruct.C),
                res = 0;
                dyn = ii;
                msg = '"C" must not be empty!';
            elseif isempty(sysStruct.D),
                res = 0;
                dyn = ii;
                msg = '"D" must not be empty!';
            elseif isempty(sysStruct.guardX) | isempty(sysStruct.guardU) | isempty(sysStruct.guardC),
                res = 0;
                dyn = ii;
                msg = 'Active region must not be empty!';
            end
        end
    end
    res = 1;
catch
    msg = lasterr;
    res = 0;
end