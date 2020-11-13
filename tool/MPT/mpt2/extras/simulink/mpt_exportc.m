function mpt_exportc(ctrl, fname)
%MPT_EXPORTC Exports an explicit controller to C code
%
% mpt_exportc(ctrl, fname)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Creates a header file which contains information about a given explicit
% controller. The header file must then be compiled with mpt_getInput.c, which
% contains the region identification procedure. Note that if you want to compile
% multiple controllers, you should use different names of the header files and
% you should also modify mpt_getInput.c to take the appropriate header.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl   - MPT explicit controller
% fname  - Name of the header file to be generated ("mpt_getInput.h" by default)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%

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
    fname = 'mpt_getInput.h';
end

Pn = ctrl.Pn;
[Hn,Kn] = double(Pn);
if ~iscell(Hn),
    Hn = {Hn};
    Kn = {Kn};
end
nctotal = 0;
for ii=1:length(Pn),
    nctotal = nctotal + nconstr(Pn(ii));
end
nr = length(Pn);
Fi = ctrl.Fi;
Gi = ctrl.Gi;
nx = ctrl.details.dims.nx;
nu = ctrl.details.dims.nu;
ny = ctrl.details.dims.ny;
nref = 0;
nxt = nx;

for ii=1:length(Fi),
    Fi{ii} = Fi{ii}(1:nu,:);
    Gi{ii} = Gi{ii}(1:nu);
end
deltau = isfield(ctrl.sysStruct, 'dumode');
tracking = ctrl.probStruct.tracking;
if tracking==1,
    nxt = ctrl.sysStruct.dims.nx;
elseif deltau,
    nxt = nx - nu;
end
if tracking>0,
    nxt = ctrl.sysStruct.dims.nx;
    if isfield(ctrl.probStruct, 'Qy'),
        nref = ctrl.sysStruct.dims.ny;
    else
        nref = ctrl.sysStruct.dims.nx;
    end
end

if isfield(ctrl.sysStruct, 'Ts'),
    Ts = ctrl.sysStruct.Ts;
else
    Ts = 1;
end

fid = fopen(fname, 'w');
if fid<0,
    error('Cannot open file for writing!');
end
fprintf(fid, '#define mpt_getInput_h\n\n');
fprintf(fid, '#define MPT_NR %d\n', nr);
fprintf(fid, '#define MPT_NX %d\n', nx);
fprintf(fid, '#define MPT_NU %d\n', nu);
fprintf(fid, '#define MPT_NY %d\n', ny);
fprintf(fid, '#define MPT_NXT %d\n', nxt);
fprintf(fid, '#define MPT_NREF %d\n', nref);
fprintf(fid, '#define MPT_TS %f\n', Ts);
fprintf(fid, '#define MPT_DUMODE %d\n', deltau);
fprintf(fid, '#define MPT_TRACKING %d\n', tracking);
fprintf(fid, '#define MPT_ABSTOL %e\n', mptOptions.abs_tol);

ctr = 0;
fprintf(fid, '\nstatic float MPT_H[] = {\n');
for ii = 1:nr,
    Hi = Hn{ii};
    nc = size(Hi, 1);
    for jj = 1:nc,
        h = Hi(jj, :);
        for kk = 1:length(h),
            ctr = ctr + 1;
            if ctr<nctotal*nx,
                fprintf(fid, '%e,\t', h(kk));
            else
                fprintf(fid, '%e ', h(kk));
            end
            if mod(ctr, 5)==0,
                fprintf(fid, '\n');
            end
        end
    end
end
fprintf(fid, '};\n\n');

ctr = 0;
fprintf(fid, 'static float MPT_K[] = {\n');
for ii = 1:nr,
    Ki = Kn{ii};
    nc = size(Ki, 1);
    for jj = 1:nc,
        ctr = ctr + 1;
        if ctr<nctotal,
            fprintf(fid, '%e,\t', Ki(jj));
        else
            fprintf(fid, '%e ', Ki(jj));
        end
        if mod(ctr, 5)==0,
            fprintf(fid, '\n');
        end
    end
end
fprintf(fid, '};\n\n');

ctr = 0;
fprintf(fid, 'static int MPT_NC[] = {\n');
for ii = 1:nr,
    if ii < nr,
        fprintf(fid, '%d,\t', nconstr(Pn(ii)));
    else
        fprintf(fid, '%d ', nconstr(Pn(ii)));
    end
    if mod(ii, 5)==0,
        fprintf(fid, '\n');
    end
end
fprintf(fid, '};\n\n');

nctotal = nx*nu*nr;
ctr = 0;
fprintf(fid, 'static float MPT_F[] = {\n');
for ii = 1:nr,
    F = Fi{ii};
    for jj = 1:nu,
        f = F(jj, :);
        for kk = 1:nx,
            ctr = ctr + 1;
            if ctr<nctotal,
                fprintf(fid, '%e,\t', f(kk));
            else
                fprintf(fid, '%e ', f(kk));
            end
            if mod(ctr, 5)==0,
                fprintf(fid, '\n');
            end
        end
    end
end
fprintf(fid, '};\n\n');

ctr = 0;
fprintf(fid, 'static float MPT_G[] = {\n');
for ii = 1:nr,
    G = Gi{ii};
    for jj = 1:nu,
        ctr = ctr + 1;
        if ctr<nctotal,
            fprintf(fid, '%e,\t', G(jj));
        else
            fprintf(fid, '%e ', G(jj));
        end
        if mod(ctr, 5)==0,
            fprintf(fid, '\n');
        end
    end
end
fprintf(fid, '};\n');

fclose(fid);
