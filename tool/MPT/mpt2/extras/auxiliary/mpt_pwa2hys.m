function mpt_pwa2hys(sysStruct, fname)
%MPT_PWA2HYS Converts a PWA system described in sysStruct into a HYSDEL model
%
% mpt_pwa2hys(sysStruct, filename)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Converts a PWA system described into a HYSDEL model.
%
% NOTE! you need to compile the output file with HYSDEL manually.
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%  sysStruct   - MPT system structure
%  filename    - name of HYSDEL file to create
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%

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

% still to fix:
%  * if sysStruct.(C,D,g) is constant for every dynamics, we can get rid of
%    auxliary ZY variables and rewrite the OUTPUT section to use state variables

fprintf('Warning: This function is experimental and known to be buggy in some cases, use with caution.\n');

sysStruct = mpt_verifySysStruct(sysStruct);

if ~iscell(sysStruct.A),
    error('System must be Piecewise-Affine!');
end

if ~isfield(sysStruct, 'xmin')
    error('"sysStruct.xmin" must be defined!');
end
if ~isfield(sysStruct, 'xmax')
    error('"sysStruct.xmax" must be defined!');
end
    

% obtain information about the system
[nx,nu,ny,ndyn,nbool,ubool] = mpt_sysStructInfo(sysStruct);
unonbool = setdiff(1:nu, ubool);

modelname = 'sysStruct2mld';

%=================================================================
% create symbolic names for states, inputs and outputs
Xstr = cell(1,nx);
for ii=1:nx,
    Xstr{ii} = sprintf('x%d', ii);
end

Ustr = cell(1,nu);
for ii=1:nu,
    Ustr{ii} = sprintf('u%d', ii);
end

Ystr = cell(1,ny);
for ii=1:ny,
    Ystr{ii} = sprintf('y%d', ii);
end



%=================================================================
% get the list of unique delta variables
Guards = [];
for idyn = 1:ndyn,
    Gx = sysStruct.guardX{idyn};
    Gu = sysStruct.guardU{idyn};
    Gc = sysStruct.guardC{idyn};
    nc = size(Gx,1);
    Guards = [Guards; Gx Gu Gc];
end
Guards = unique(Guards, 'rows');
nG = size(Guards, 1);
oppositeGuards = [];

for ii=1:nG-1,
    for jj=ii+1:nG,
        G1 = Guards(ii,:);
        G2 = Guards(jj,:);
        if all(G1==-G2),
            oppositeGuards = [oppositeGuards; ii jj];
        end
    end
end
for ii=1:nG,
    if isempty(find(oppositeGuards==ii)),
        oppositeGuards = [oppositeGuards; ii 0];
    end
end

%=================================================================
% define MUST section

% look if there are guardlines which occur in every single dynamics. if
% so, these constraints will go to the must section
repetitiveguards = zeros(nG,1);
remove_opposite = [];
MUSTsection = '';
for idyn = 1:ndyn,
    RG=sub_regionGuards(sysStruct, Guards, idyn);
    for ii = 1:length(RG),
        repetitiveguards(RG(ii)) = repetitiveguards(RG(ii)) + 1;
    end
end
must_guards_ind = find(repetitiveguards==ndyn);
if ~isempty(must_guards_ind)
    must_guards = Guards(must_guards_ind,:);
    if ~isempty(ubool),
        % check if these guards are defined only on continuous inputs
        for ii = 1:length(must_guards_ind),
            guard = must_guards(ii,:);
            gX = guard(1:nx);
            gU = guard(nx+1:nx+nu);
            gC = guard(end);
            if any(gU(ubool)~=0)
                error('Problem: non-zero guards on boolean inputs in MUST section! please contact mpt@control.ee.ethz.ch for help.');
            end
        end
    end

    for ii = 1:length(must_guards_ind),
        [row, col] = find(oppositeGuards==must_guards_ind(ii));
        if ~isempty(row)
            if col==1,
                % all ok, remove this guardline from oppositeGuards
                remove_opposite = [remove_opposite; row];
            else
                % this guardline has an "opposite" guardline, check if that one
                % belongs to the must section as well
                oppguard = oppositeGuards(row, 1);
                if isempty(find(must_guards_ind==oppguard))
                    % if not, we have a problem
                    error('Problem: a MUST section guardline has an opposite which does not belong to MUST section! please contact mpt@control.ee.ethz.ch for help.');
                else
                    % all ok, remove this guardline from oppositeGuards
                    remove_opposite = [remove_opposite; row];
                end
            end
        end
    end
    remove_opposite = reshape(oppositeGuards(remove_opposite,:), length(remove_opposite)*2, 1);
    remove_opposite(find(remove_opposite==0)) = [];
    for ii = 1:length(remove_opposite)
        guardindex = remove_opposite(ii);
        guard = Guards(guardindex,:);
        
        gX = guard(1:nx);
        gU = guard(nx+1:nx+nu);
        gC = guard(end);
        Gstr = '';
        for ix = 1:nx,
            Gstr = [Gstr sprintf('(%f * %s) + ', gX(ix), Xstr{ix})];
        end
        for iu = 1:unonbool,
            Gstr = [Gstr sprintf('(%f * %s)', gU(iu), Ustr{iu})];
            if iu<length(unonbool),
                Gstr = [Gstr ' + '];
            end
        end
        Gstr = [Gstr ' <= ' sprintf('%f;', gC)];
        MUSTsection = strvcat(MUSTsection, Gstr);
    end
        
end


%=================================================================
% define AD section

% number of uniqe delta variables
uniqueDeltas = size(oppositeGuards, 1);
uniqueDeltas = 0;
for ii = 1:size(oppositeGuards,1)
    if ~isempty(remove_opposite)
        if ~isempty(find(remove_opposite==ii))
            % this guardline belongs to must section, skip it here
            continue
        end
    end
    uniqueDeltas = uniqueDeltas + 1;
end
    
Dstr = cell(1,uniqueDeltas);
for ii=1:uniqueDeltas,
    Dstr{ii} = sprintf('d%d', ii);
end
    
ADsection = '';
for ii=1:uniqueDeltas,
    
    guard = Guards(oppositeGuards(ii,1),:);
    gX = guard(1:nx);
    gU = guard(nx+1:nx+nu);
    gC = guard(end);
    Gstr = sprintf('%s = ', Dstr{ii});
    
    % first look if there are some deltas which depend on binary inputs. if so,
    % no guards on states are allowed
    if ~isempty(ubool),
        if any(gU(ubool)~=0),
            nonbool = setdiff(1:nu, ubool);
            if any(gX~=0)
                error('Sorry, non-zero state guardlines not allowed in combination with boolean inputs!');
            elseif any(gU(nonbool)~=0)
                error('Sorry, non-zero guardlines on continuous states not allowed in combination with guards on boolean inputs!');
            end
        end
    end

    if ~isempty(ubool) & any(gU(ubool)~=0)
        % detected a guard associated to boolean input, handle it properly. keep
        % it mind that the conditions above exclude a case when there are also
        % additional guards on state and/or continuous inputs
        Gstr = '';
        for iu = ubool,
            if abs(gC)>0.5,
                Gstr = [Gstr sprintf('(%s)', Ustr{iu})];
            else
                Gstr = [Gstr sprintf('(~%s)', Ustr{iu})];
            end
            if iu < max(ubool)
                Gstr = [Gstr ' & '];
            end
        end
        if ~isempty(find(oppositeGuards(:,1)==ii))
            % remove the "opposite" guardline to this guardline
            oppositeGuards(ii,2) = 0;
        end
        Dstr{ii} = Gstr;
        
    else
        for ix = 1:nx,
            Gstr = [Gstr sprintf('(%f * %s) + ', gX(ix), Xstr{ix})];
        end
        for iu = 1:unonbool,
            Gstr = [Gstr sprintf('(%f * %s)', gU(iu), Ustr{iu})];
            if iu<length(unonbool),
                Gstr = [Gstr ' + '];
            end
        end
        Gstr = [Gstr ' <= ' sprintf('%f;', gC)];
        ADsection = strvcat(ADsection, Gstr);
    end

end

%=================================================================
% define DA section - auxiliary variables for state update

ZXstr = {};
zcount = 0;
DAsection = '';

for idyn=1:ndyn,
    RG = sub_regionGuards(sysStruct, Guards, idyn);
    
    % RG now holds indices to 'Guards'
    for ix=1:nx,
        zcount = zcount + 1;
        zxstr = sprintf('zx%d_%d', ix, zcount);
        ZXstr{ix}{idyn} = zxstr;
        dvarname = '(';
        for jj=1:length(RG),
            rg = RG(jj);

            if ~isempty(find(remove_opposite==rg)),
                continue
            end
            
            if isempty(find(oppositeGuards==rg))
                % this situation can happen because we removed "opposite"
                % guardlines which are defined over boolean inputs. in that case
                % we don't need to process this guardline, since the definition
                % of it is already comprised in another DA variable
                continue
            end

            [row, col] = find(oppositeGuards==rg);
            
            if isempty(find(oppositeGuards(:,2)==rg))
                % this delta variable has no 'opposite' term
                dvarname = [dvarname Dstr{row}];
            else
                % this delta variable has an 'opposite' term, i.e. it is defined
                % as a negation of some other delta variable
                dvarname = [dvarname '(~' Dstr{row} ')'];
            end
            dvarname = [dvarname ' & '];
        end
        if length(dvarname)==1
            error(sprintf('No affine update equation for dynamics %d. This should not happen, please contact mpt@control.ee.ethz.ch for help.', idyn));
        end
        if dvarname(end-1:end)=='& ',
            dvarname = dvarname(1:end-3);
        end
        dvarname = [dvarname ')'];
        zupdate = sub_getZupdate(sysStruct, idyn, ix, Xstr, Ustr);
        DAline = sprintf('%s = { IF %s THEN %s ELSE 0 };', ...
            zxstr, dvarname, zupdate);        
        DAsection = strvcat(DAsection, DAline);
    end
end



%=================================================================
% define DA section - auxiliary variables for outputs
clear ix
ZYstr = {};
zcount = 0;

for idyn=1:ndyn,
    RG = sub_regionGuards(sysStruct, Guards, idyn);
    
    % RG now holds indices to 'Guards'
    for iy=1:ny,
        zcount = zcount + 1;
        zystr = sprintf('zy%d_%d', iy, zcount);
        ZYstr{iy}{idyn} = zystr;
        dvarname = '(';
        for jj=1:length(RG),
            rg = RG(jj);
            
            if ~isempty(find(remove_opposite==rg)),
                continue
            end

            if isempty(find(oppositeGuards==rg))
                % this situation can happen because we removed "opposite"
                % guardlines which are defined over boolean inputs. in that case
                % we don't need to process this guardline, since the definition
                % of it is already comprised in another DA variable
                continue
            end
            [row, col] = find(oppositeGuards==rg);
            if isempty(find(oppositeGuards(:,2)==rg))
                % this delta variable has no 'opposite' term
                dvarname = [dvarname Dstr{row}];
            else
                % this delta variable has an 'opposite' term, i.e. it is defined
                % as a negation of some other delta variable
                dvarname = [dvarname '(~' Dstr{row} ')'];
            end
            dvarname = [dvarname ' & '];
        end
        
        if length(dvarname)==1
            error(sprintf('No affine update equation for dynamics %d. This should not happen, please contact mpt@control.ee.ethz.ch for help.', idyn));
        end
        if dvarname(end-1:end)=='& ',
            dvarname = dvarname(1:end-3);
        end
        
        dvarname = [dvarname ')'];
        zupdate = sub_getYupdate(sysStruct, idyn, iy, Xstr, Ustr);
        DAline = sprintf('%s = { IF %s THEN %s ELSE 0 };', ...
            zystr, dvarname, zupdate);        
        DAsection = strvcat(DAsection, DAline);
    end
end



%=================================================================
% define CONTINUOUS section
CONTsection = '';
ctr = 0;

for ix=1:nx,
    xvar = Xstr{ix};
    contline = sprintf('%s = ', xvar);
    for idyn=1:ndyn,
        zvar = ZXstr{ix}{idyn};
        contline = [contline zvar];
        if idyn<ndyn
            contline = [contline ' + '];
        else
            contline = [contline ';'];
        end
    end
    CONTsection = strvcat(CONTsection, contline);
end



%=================================================================
% define OUTPUT section
OUTsection = '';
ctr = 0;
for iy=1:ny,
    yvar = Ystr{iy};
    outline = sprintf('%s = ', yvar);
    for idyn=1:ndyn,
        zvar = ZYstr{iy}{idyn};
        outline = [outline zvar];
        if idyn<ndyn
            outline = [outline ' + '];
        else
            outline = [outline ';'];
        end
    end
    
    OUTsection = strvcat(OUTsection, outline);
end






HYS = '';
HYS = strvcat(HYS, sprintf('SYSTEM %s {', modelname));
HYS = strvcat(HYS, 'INTERFACE {');
HYS = strvcat(HYS, 'STATE {');

% write states
for ii=1:nx,
    HYS = strvcat(HYS, sprintf('REAL %s [%f, %f];', Xstr{ii}, sysStruct.xmin(ii), sysStruct.xmax(ii)));
end

HYS = strvcat(HYS, '}'); % end of STATE
HYS = strvcat(HYS, 'INPUT {');

% decide which inputs are boolean and which real
utype = repmat('R', nu, 1);
if isfield(sysStruct, 'Uset')
    for ii=1:nu
        if ~iscell(sysStruct.Uset),
            uset = sysStruct.Uset;
        else
            uset = sysStruct.Uset{ii};
        end
        if any(~isinf(uset))
            % either an integer or boolean input
            if uset==[0 1] | uset==[1 0]
                % boolean input
                utype(ii) = 'B';
            else
                error(sprintf('%s.Uset{%d} must be a boolean input!', inputname(1), ii));
            end
        end
    end
end

% write inputs
for ii=1:nu,
    if utype(ii)=='R'
        HYS = strvcat(HYS, sprintf('REAL %s [%f, %f];', Ustr{ii}, sysStruct.umin(ii), sysStruct.umax(ii)));
    else
        HYS = strvcat(HYS, sprintf('BOOL %s;', Ustr{ii}));
    end
end

HYS = strvcat(HYS, '}'); % end of INPUT
HYS = strvcat(HYS, 'OUTPUT {');

% write outputs
for ii=1:ny,
    HYS = strvcat(HYS, sprintf('REAL %s;', Ystr{ii}));
end

HYS = strvcat(HYS, '}'); % end of OUTPUT
HYS = strvcat(HYS, '}'); % end of INTERFACE
HYS = strvcat(HYS, 'IMPLEMENTATION {');
HYS = strvcat(HYS, 'AUX {');

% write auxs
for ii=1:nx,
    for jj=1:ndyn
        HYS = strvcat(HYS, sprintf('REAL %s;', ZXstr{ii}{jj}));
    end
end
for ii=1:ny,
    for jj=1:ndyn
        HYS = strvcat(HYS, sprintf('REAL %s;', ZYstr{ii}{jj}));
    end
end
for ii=1:length(Dstr),
    if isempty(findstr(Dstr{ii}, 'u')),
        % only write boolean variables if they are are not associated to boolean
        % inputs
        HYS = strvcat(HYS, sprintf('BOOL %s;', Dstr{ii}));
    end
end


HYS = strvcat(HYS, '}'); % end of AUX
HYS = strvcat(HYS, 'AD {');

% write ads
HYS = strvcat(HYS, ADsection);
HYS = strvcat(HYS, '}'); % end of AD
HYS = strvcat(HYS, 'DA {');

% write das
HYS = strvcat(HYS, DAsection);

HYS = strvcat(HYS, '}'); % end of DA
HYS = strvcat(HYS, 'CONTINUOUS {');

% write continuous
HYS = strvcat(HYS, CONTsection);

HYS = strvcat(HYS, '}'); % end of CONTINUOUS
HYS = strvcat(HYS, 'OUTPUT {');

% write outputs
HYS = strvcat(HYS, OUTsection);

HYS = strvcat(HYS, '}'); % end of OUTPUT

% write must section

if ~isempty(MUSTsection),
    HYS = strvcat(HYS, 'MUST {');
    HYS = strvcat(HYS, MUSTsection);
    HYS = strvcat(HYS, '}');
end

HYS = strvcat(HYS, '}'); % end IMPLEMENTATION
HYS = strvcat(HYS, '}'); % end SYSTEM

if isempty(findstr(fname, '.hys')),
    fname = [fname '.hys'];
end
fid = fopen(fname, 'w');
if fid<0
    error(['cannot open file "' fname '" for writing!']);
end
for ii=1:size(HYS,1)
    % remove spaces after text (due to strvcat);
    line = deblank(HYS(ii,:));
    fprintf(fid, '%s\n', line);
end
fclose(fid);
return

%------------------------------------------------------------------------
function Gstr = sub_getZupdate(sysStruct, dyn, x, Xstr, Ustr);

[nx,nu,ny,ndyn,nbool,ubool] = mpt_sysStructInfo(sysStruct);

A = sysStruct.A{dyn};
B = sysStruct.B{dyn};
F = sysStruct.f{dyn};

a = A(x,:);
b = B(x,:);
f = F(x);

nx = length(Xstr);
nu = length(Ustr);

Gstr = '';
for ix = 1:nx,
    Gstr = [Gstr sprintf('(%f * %s) + ', a(ix), Xstr{ix})];
end
for iu = 1:nu,
    if ~isempty(find(ubool==iu))
        % this input is a boolean one, if the correspondig coefficient is
        % non-zero, we have a problem, recasting with REAL does not fully work.
        if b(iu)~=0
            error('Affine update equation depends on binary input, this case is not handles yet, sorry');
        end
        %Gstr = [Gstr sprintf('(%f * (REAL %s)) + ', b(iu), Ustr{iu})];
    else
        Gstr = [Gstr sprintf('(%f * %s) + ', b(iu), Ustr{iu})];
    end
end
Gstr = [Gstr sprintf('(%f)', f)];


%------------------------------------------------------------------------
function Gstr = sub_getYupdate(sysStruct, dyn, y, Xstr, Ustr);

[nx,nu,ny,ndyn,nbool,ubool] = mpt_sysStructInfo(sysStruct);

C = sysStruct.C{dyn};
D = sysStruct.D{dyn};
G = sysStruct.g{dyn};

c = C(y,:);
d = D(y,:);
g = G(y);

nx = length(Xstr);
nu = length(Ustr);

Gstr = '';
for ix = 1:nx,
    Gstr = [Gstr sprintf('(%f * %s) + ', c(ix), Xstr{ix})];
end
for iu = 1:nu,
    if ~isempty(find(ubool==iu))
        % this input is a boolean one, if the correspondig coefficient is
        % non-zero, we have a problem, recasting with REAL does not fully work.
        if d(iu)~=0
            error('Affine update equation depends on binary input, this case is not handles yet, sorry');
        end
        %Gstr = [Gstr sprintf('(%f * (REAL %s)) + ', d(iu), Ustr{iu})];
    else
        Gstr = [Gstr sprintf('(%f * %s) + ', d(iu), Ustr{iu})];
    end
end
Gstr = [Gstr sprintf('(%f)', g)];


%------------------------------------------------------------------------
function RG=sub_regionGuards(sysStruct, Guards, idyn)
% returns indices of guardlines defined in 'Guards' which form guarlines for
% dynamics 'idyn' of 'sysStruct'

gX = sysStruct.guardX{idyn};
gU = sysStruct.guardU{idyn};
gC = sysStruct.guardC{idyn};
G = unique([gX gU gC], 'rows');
RG = [];

for ii=1:size(G,1)
    for jj=1:size(Guards,1)
        if G(ii,:)==Guards(jj,:)
            RG = [RG; jj];
        end
    end
end
