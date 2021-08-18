function mpt_exportST(ctrl, fname)
%MPT_EXPORTST Exports a search tree to C code
%
% mpt_exportST(ctrl, fname)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Creates a standalone C-file which contains definition of the search tree and
% code which identifies control action associated to a given state.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl   - MPT explicit controller
% fname  - Name of the header file to be generated ("searchtree.c" by default)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%

% Copyright is with the following author(s):
%
% (C) 2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

error(nargchk(1,2,nargin)); 

global mptOptions
if ~isstruct(mptOptions)
    mpt_error;
end

if ~isa(ctrl, 'mptctrl')
    error('Input must be an MPTCTRL controller object!');
end
if ~cancompile(ctrl)
    error('This controller cannot be exported to C code.');
end
if nargin<2,
    fname = 'searchtree.c';
end
if ~isfield(ctrl.details, 'searchTree'),
    error('Search tree is not present in the controller. Use "mpt_searchTree(ctrl)" to compute it.');
end

template_file = which('mpt_searchTree.c');
if isempty(template_file),
    error('Cannot find file "mpt_searchTree.c". Check your path setup.');
end

ST = ctrl.details.searchTree'; ST = ST(:);
Fi = ctrl.Fi;
Gi = ctrl.Gi;
nx = ctrl.details.dims.nx;
nu = ctrl.details.dims.nu;
nr = length(ctrl.Pn);

% only keep the first element of the open-loop sequence of inputs
for ii=1:length(Fi),
    Fi{ii} = Fi{ii}(1:nu,:);
    Gi{ii} = Gi{ii}(1:nu);
end

out = '';
out = strvcat(out, sprintf('#define MPT_NU %d', nu));
out = strvcat(out, sprintf('#define MPT_NX %d', nx));

template = fileread(template_file);

ctr = 0;
out = strvcat(out, 'static float MPT_ST[] = {');

o = '';
for ii = 1:length(ST),
    o = [o sprintf('%e,\t', ST(ii))];
    if mod(ii, 5)==0 | ii==length(ST),
        out = strvcat(out, o);
        o = '';
    end
    
end
out = strvcat(out, '};');

out = strvcat(out, 'static float MPT_F[] = {');
for ii = 1:nr,
    F = Fi{ii};
    for jj = 1:nu,
        f = F(jj, :);
        o = '';
        for kk = 1:length(f),
            o = [o sprintf('%e,\t', f(kk))];
        end
        out = strvcat(out, o);
    end
end
out = strvcat(out, '};');

out = strvcat(out, 'static float MPT_G[] = {');
for ii = 1:nr,
    F = Gi{ii};
    for jj = 1:nu,
        f = F(jj, :);
        o = '';
        for kk = 1:length(f),
            o = [o sprintf('%e,\t', f(kk))];
        end
        out = strvcat(out, o);
    end
end
out = strvcat(out, '};');

% convert the concatenated array into one string, add line breaks
out_nl = [];
for ii = 1:size(out, 1),
    out_nl = [out_nl sprintf('%s\n', deblank(out(ii, :)))];
end

% now replace the placeholder in mpt_searchTree.c
template = strrep(template, '/* placeholder, do not edit or remove!!! */', out_nl);

outfid = fopen(fname, 'w');
if outfid < 0,
    error(sprintf('Cannot open file "%s" for writing!', fname));
end
fprintf(outfid, '%s', template);
fclose(outfid);

fprintf('Output written to "%s".\n', fname);
