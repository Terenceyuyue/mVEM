function [newCtrlStruct]=mpt_removeOverlaps(Partition,Options)
%MPT_REMOVEOVERLAPS Removes overlaps from (a set of) polyhedral partitions with associated linear cost
%
% [newCtrlStruct]=mpt_removeOverlaps(ctrl,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Given 'n' (possibly overlapping) polyhedral partitions with associated linear cost,
% this function detect such overlaps and removes them by picking up regions which
% have the least cost associated to them. If intersection of the overlapping regions
% is non-empty, slicing is introduced to create one non-overlapping polyhedral
% partition.
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl        - Explicit controller a cell aray thereof; with fields:
%      Pfinal - maximum feasible set (i.e. the outer hull of U Pn)
%      Pn     - polytopic regions
%      Fi,Gi  - cells containing control law (u = Fi{j} x + Gi{j})
%      Bi,Ci  - cells containing the linear cost (i.e. J = Bi{j} x + Ci{j})
% Options.verbose   - level of verbosity
% Options.lpsolver  - default solver for LP's
% Options.abs_tol   - absolute tolerance
% Options.rel_tol   - relative tolerance
% Options.Vpruning  - whether or not to prune non-valid intersections by
%                     examining extreme points (0 is default)
% Options.lowmem    - defines memory saving mode
%                       0 - no memory saving - fast computation (default)
%                       1 - heavy memory saving (slow computation)
%
% Note: If Options is missing or some of the fiels are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% newCtrl    - controller which consists of purely non-overlapping regions
%
% see also POLYTOPE/REDUCEUNION, POLYTOPE/UNIQUE
%

% Copyright is with the following author(s):
%
% (C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch

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

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<2,
    Options=[];
end


%====================================================================
% set default options and parameters
Options = mpt_defaultOptions(Options, ...
    'abs_tol', mptOptions.abs_tol, ...
    'rel_tol', mptOptions.rel_tol, ...
    'bbox_tol', 1e4*mptOptions.abs_tol, ...
    'verbose', mptOptions.verbose, ...
    'lpsolver', mptOptions.lpsolver, ...
    'lowmem', 0, ...
    'Vpruning', 0, ...
    'sphratio', 0.1 );

abs_tol = Options.abs_tol;
rel_tol = Options.rel_tol;
bbox_tol = Options.bbox_tol;
lpsolver = Options.lpsolver;
emptypoly = mptOptions.emptypoly;
lowmem = Options.lowmem;
bboxOpt = Options;
bboxOpt.noPolyOutput = 1;
mptctrl_input = 0;
input_is_polytope = 0;


%====================================================================
% decide type of input argument
if isa(Partition, 'polytope'),
    % input is a polytope array, convert it to a dummy controller structure
    input_is_polytope = 1;
    Partition = mpt_dummyCS(Partition);
else
    input_is_polytope = 0;
end 

if ~iscell(Partition)
    if isa(Partition, 'mptctrl')
        if ~isexplicit(Partition)
            error('This function supports only explicit controllers!');
        end
        Partition = struct(Partition);
        ctrlStruct = Partition;
        mptctrl_input = 1;
    end
    if ~mpt_isValidCS(Partition),
        error('Input argument must be either a cell array containing polyhedral partitions or a valid ctrlStruct structure!');
    end
    % convert ctrlStruct into "Partition" cells (see help above)
    
    PPfinal = Partition.Pfinal;
    Options.noSectionCheck = 1;
    ctrlStruct = Partition;
    if ~ctrlStruct.overlaps,
        disp('mpt_removeOverlaps: Partitions do not overlap, nothing to do here...');
        newCtrlStruct = ctrlStruct;
        newCtrlStruct.details.keptParts = 1;
        if mptctrl_input
            newCtrlStruct = mptctrl(newCtrlStruct);
        end
        return
    end
    if ctrlStruct.probStruct.subopt_lev == 1,
        newCtrlStruct = sub_removeOverlaps_MinTime(ctrlStruct);
        newCtrlStruct.details.keptParts = 1;        
        if mptctrl_input
            newCtrlStruct = mptctrl(newCtrlStruct);
        end
        return
    end
    sysStruct = Partition.sysStruct;
    probStruct = Partition.probStruct;
    details = Partition.details;
    Partition = cell(1,length(ctrlStruct.Pn));
    for ii=1:length(ctrlStruct.Pn),
        Partition{ii}.Pn = ctrlStruct.Pn(ii);
        Partition{ii}.Pfinal = Partition{ii}.Pn;
        Partition{ii}.Fi{1} = ctrlStruct.Fi{ii};
        Partition{ii}.Gi{1} = ctrlStruct.Gi{ii};
        Partition{ii}.Ai{1} = ctrlStruct.Ai{ii};
        Bi = ctrlStruct.Bi{ii};
        Partition{ii}.Bi{1} = (Bi(:))';
        Partition{ii}.Ci{1} = ctrlStruct.Ci{ii};
        Partition{ii}.dynamics = ctrlStruct.dynamics(ii);
    end
    if mptctrl_input,
        Partition{1}.sysStruct = ctrlStruct.sysStruct;
        Partition{1}.probStruct = ctrlStruct.probStruct;
        Partition{1}.details = ctrlStruct.details;
    end
else
    Options.noSectionCheck = 0;
    PPfinal = [];
    for ii=1:length(Partition),
        if isa(Partition{ii}, 'mptctrl')
            if ~isexplicit(Partition{ii})
                error('This function supports only explicit controllers!');
            end
            Partition{ii} = struct(Partition{ii});
            mptctrl_input = 1;
        end
            
        if ~mpt_isValidCS(Partition{ii}),
            error('Each element of the cell array has to be a valid ctrlStruct structure!');
        end
        PPfinal = [PPfinal Partition{ii}.Pfinal];
        if ~isfield(Partition{ii},'dynamics')
            Partition{ii}.dynamics = zeros(1,length(Partition{ii}.Pn));
        end
    end
    sysStruct = Partition{1}.sysStruct;
    probStruct = Partition{1}.probStruct;
    details = Partition{1}.details;
end


%====================================================================
% check whether we have quadratic cost terms
npart = length(Partition);
Ai = {};
for i = 1:npart,
    Ai = cat(2, Ai, Partition{i}.Ai);
end
if nnz([Ai{:}])>0,
    error('Quadratic cost terms are not allowed.');
end


%====================================================================
% stack bounding boxes of all partitions' Pfinal polytopes together
% for faster access later
nx = dimension(Partition{1}.Pn);
PartBBoxes = zeros(nx, npart*2);
PartLength = zeros(1, npart);
for ii = 1:npart,
    [d, bmin, bmax] = bounding_box(Partition{ii}.Pfinal, bboxOpt);
    PartBBoxes(:, (ii-1)*2+1:ii*2) = [bmin bmax];
    PartLength(ii) = length(Partition{ii}.Pn);
end
% prepare an intersection table which contains information about which
% partitions do intersect
PartIntersect = ones(npart);
for ipart = 1:npart-1,
    PartBBi = PartBBoxes(:, (ipart-1)*2+1:ipart*2);
    for jpart = ipart+1:npart,
        PartBBj = PartBBoxes(:, (jpart-1)*2+1:jpart*2);
        if any(PartBBi(:,2) + bbox_tol < PartBBj(:,1)) | ...
                any(PartBBj(:,2) + bbox_tol < PartBBi(:,1)),
            % the two partitions for sure do not intersect since their bounding
            % boxes do not intersect
            PartIntersect(ipart, jpart) = 0;
            PartIntersect(jpart, ipart) = 0;
        end
    end
end
% we don't need the bounding boxes anymore, clear them to save memory
clear PartBBoxes



%====================================================================
% look for duplicate regions, i.e. regions which have identical cost, control
% law and regions are the same. if such regions is identified, remove them from
% Partition list (be carefull not to remove both identical regions, just one of
% them!)
nRemoved = 0;
for ipart = 1:npart-1,
    Partition_i = Partition{ipart};
    for jpart = ipart+1:npart,
        if PartIntersect(ipart, jpart)==0,
            % partitions do not intersect, hence they cannot share identical
            % regions
            continue
        end
        Partition_j = Partition{jpart};
        Cj = [Partition_j.Ci{:}];
        if isempty(Cj),
            % partition was already removed, skip to next one
            continue
        end
        ij_equal = [];
        for ireg = 1:PartLength(ipart),
            Fi = Partition_i.Fi{ireg};
            Gi = Partition_i.Gi{ireg};
            Bi = Partition_i.Bi{ireg};
            Ci = Partition_i.Ci{ireg};
            for jreg = 1:PartLength(jpart),
                % check if regions Partition{jpart}(jreg) and
                % Partition{ipart}(ireg) are identical. if so, we remove one of
                % them
                %
                % we could also do this comparison in a subfunctions, but it
                % turned out that doing it in-line is much faster.
                if abs(Ci-Cj(jreg)) <= abs_tol,
                    if sum(abs(Gi-Partition_j.Gi{jreg}), 2) <= abs_tol,
                        if sum(abs(Bi-Partition_j.Bi{jreg}), 2) <= abs_tol,
                            if sum(abs(Fi-Partition_j.Fi{jreg}), 2) <= abs_tol,
                                if Partition_i.Pn(ireg) == Partition_j.Pn(jreg),
                                    ij_equal = [ij_equal; jreg];
                                end
                            end
                        end
                    end
                end
            end
        end
        % only keep those regions which are not identical
        j_unique = unique(ij_equal);
        nRemoved = nRemoved + length(j_unique);
        jkeep = setdiff(1:PartLength(jpart), j_unique);
        Partition_j.Pn = Partition_j.Pn(jkeep);
        Partition_j.Fi = {Partition_j.Fi{jkeep}};
        Partition_j.Gi = {Partition_j.Gi{jkeep}};
        Partition_j.Bi = {Partition_j.Bi{jkeep}};
        Partition_j.Ci = {Partition_j.Ci{jkeep}};
        Partition{jpart} = Partition_j;
        PartLength(jpart) = length(Partition_j.Pn);
    end
end
if nRemoved>0 & Options.verbose > 0,
    disp(sprintf('Removed %d duplicate regions',nRemoved));
end



%====================================================================
% it could be that the "removing of identical regions" procedure identified
% a whole partition to be identical to an another one, thus we need to
% remove such partition(s) from the list
keep_partitions = zeros(1, npart);
for ii=1:npart,
    keep_partitions(ii) = isfulldim(Partition{ii}.Pn);
end
keep_indices = find(keep_partitions);
keptParts = keep_indices;
% keep only non-redundant partitions
Partition = {Partition{keep_indices}};
% keep only non-redundant partitions in the interesection map
PartIntersect = PartIntersect(keep_indices,keep_indices);
% keep only non-redundant partitions in the length vector
PartLength = PartLength(keep_indices);

npart = length(Partition);



%====================================================================
% exit immediatelly if only one partition left
if length(Partition)==1,
    newCtrlStruct = Partition{1};
    newCtrlStruct.details.keptParts = 1;
    if mptctrl_input
        newCtrlStruct = mptctrl(newCtrlStruct);
    end
    return
end



%====================================================================
% start removing overlaps
if Options.verbose>0,
    disp('Removing overlaps...');
end
nR = 0;
Pn = emptypoly;
Fi = {};
Gi = {};
Ai = {};
Bi = {};
Ci = {};
dynamics = []; 
keptPartsIdx = []; 



%====================================================================
% prepare bounding boxes of each region of each partition as big matrices
BBOXES = cell(1, npart);
for ipart = 1:npart,
    B = pelemfun(@bounding_box, Partition{ipart}.Pn, struct('Voutput', 1));
    BBOXES{ipart} = [B{:}];
end


%====================================================================
% prepare H-representation as a cell array (only for Options.lowmem==0)
if lowmem==0,
    HPART = cell(1, npart);
    KPART = cell(1, npart);
    for ipart = 1:npart,
        [HPART{ipart}, KPART{ipart}] = pelemfun(@double, Partition{ipart}.Pn);
    end
end


%====================================================================
% prepare map of intersecting regions
IMAP = {};
for ipart = 1:npart-1,
    for jpart = ipart+1:npart,
        if PartIntersect(ipart, jpart),
            if Options.Vpruning,
                % use pruning based on vertices - doesn't help too much compared
                % to using just bounding boxes and is much more prone to
                % numerical issues
                [IMAP{ipart}{jpart}, Partition{ipart}.Pn, Partition{jpart}.Pn] = ...
                    imap(Partition{ipart}.Pn, Partition{jpart}.Pn, bbox_tol, Options.sphratio, lpsolver);
            else
                % use pruning based on bounding boxes
                IMAP{ipart}{jpart} = bboxmap(BBOXES{ipart}, BBOXES{jpart}, bbox_tol);
            end
        end
    end
end
% we don't need the bounding boxes anymore, clear them to save memory
clear BBOXES



for ipart = 1:npart,

    if Options.verbose < 2 & Options.verbose > -1,
        if mod(ipart, round(npart/3))==0 | ipart==1 | ipart==npart,
            fprintf('hull %d/%d\t', ipart, npart);
        end
    end

    for ireg = 1:PartLength(ipart),
        
        Intersection=[];  % intersection information
        Intersection.Pn = emptypoly;
        Intersection.nR=0; 

        if lowmem == 0,
            Hi = HPART{ipart}{ireg};
            Ki = KPART{ipart}{ireg};
        else
            [Hi, Ki] = double(Partition{ipart}.Pn(ireg));
        end
        
        for jpart = 1:npart,

            if ipart==jpart,
                % do not consider overlaps within of the same partition
                continue;
                
            elseif ~PartIntersect(ipart, jpart),
                % partitions do not intersect, don't bother and skip to next one
                continue
                
            elseif ipart < jpart,
                % extract proper intersection map
                IM = IMAP{ipart}{jpart};
                
            else
                % extract proper intersection map                
                IM = IMAP{jpart}{ipart}';
                
            end
            
            for jreg = 1:PartLength(jpart),

                if IM(ireg, jreg)==0,
                    % no intersection exists between these two regions
                    continue
                end

                if lowmem==0,
                    Hj = HPART{jpart}{jreg};
                    Kj = KPART{jpart}{jreg};
                else                    
                    [Hj, Kj] = double(Partition{jpart}(jreg));
                end
                
                Haux = [Hi; Hj];
                Kaux = [Ki; Kj];
                
                [x, R] = sub_chebyball(Haux, Kaux, nx, rel_tol, abs_tol, lpsolver);
                dosect = (R > abs_tol);
                
                if dosect,

                    % intersection exists
                    Baux = Partition{jpart}.Bi{jreg} - Partition{ipart}.Bi{ireg};
                    Caux = Partition{jpart}.Ci{jreg} - Partition{ipart}.Ci{ireg};
                    % if costs are the same for two regions be careful not to remove both regions
                    % here we keep only the region belonging to the first partition
                    if sum(abs([Baux Caux]),2)<=abs_tol       % is the cost identical?
                        if ipart < jpart,                     % if so, we just remove the region once
                            Intersection.nR=Intersection.nR+1;
                            Intersection.who(Intersection.nR,:)=[jpart jreg];
                            Paux = polytope(Haux, Kaux, 0, 2);
                            Intersection.Pn = [Intersection.Pn Paux];
                        end
                    else
                        % find a part of the region in which cost associated to index 'mm' is lower than the one of index 'jj'
                        H = [Haux; Baux];
                        K = [Kaux; -Caux];
                        [x, R] = sub_chebyball(H, K, nx, rel_tol, abs_tol, lpsolver);
                        if (R >= abs_tol),
                            Paux = polytope(H, K, 0, 2, x, R);
                            Intersection.nR = Intersection.nR+1;
                            Intersection.who(Intersection.nR,:) = [jpart jreg];
                            Intersection.Pn = [Intersection.Pn Paux];
                        end
                    end
                else
                    % mark regions as non-intersecting in the intersection map
                    IM(ireg, jreg)=0;
                    
                end %dosect
            end % jreg
            
            % write updated transition map back
            if ipart < jpart,
                % extract proper intersection map
                 IMAP{ipart}{jpart} = IM;
             else
                % extract proper intersection map                
                IMAP{jpart}{ipart} = IM';
            end
            
        end % jpart
        
        if Intersection.nR > 0
            % get all segments which have the same 'minimal' cost
            Ri = regiondiff(Partition{ipart}.Pn(ireg), Intersection.Pn, Options);
            if ~isfulldim(Ri) % union is covering Ri
                if Options.verbose>1,
                    disp(['      reg ' num2str([ipart ireg]) '     ']);
                end
                continue;
            else
                if Options.verbose>1,
                    disp([' ==== reg ' num2str([ipart ireg]) ' ==== kept']);
                end
                Pn = [Pn Ri];
                for mm=1:length(Ri)
                    nR = nR+1;
                    Fi{nR} = Partition{ipart}.Fi{ireg};
                    Gi{nR} = Partition{ipart}.Gi{ireg};
                    Ai{nR} = Partition{ipart}.Ai{ireg};
                    Bi{nR} = Partition{ipart}.Bi{ireg};
                    Ci{nR} = Partition{ipart}.Ci{ireg};
                    dynamics = [dynamics Partition{ipart}.dynamics(ireg)];
                    if ~isempty(keptPartsIdx),
                        if ( ~any(keptPartsIdx == ipart) ),
                            keptPartsIdx(end+1) = ipart;
                        end
                    else
                        keptPartsIdx = ipart;
                    end
                end
            end
        else
            if Options.verbose>1,
                disp([' ==== reg ' num2str([ipart ireg]) ' ==== kept']);
            end
            nR = nR+1;
            Pn = [Pn Partition{ipart}.Pn(ireg)];
            Fi{nR} = Partition{ipart}.Fi{ireg};
            Gi{nR} = Partition{ipart}.Gi{ireg};
            Ai{nR} = Partition{ipart}.Ai{ireg};
            Bi{nR} = Partition{ipart}.Bi{ireg};
            Ci{nR} = Partition{ipart}.Ci{ireg};
            dynamics = [dynamics Partition{ipart}.dynamics(ireg)];
            if ~isempty(keptPartsIdx),
                if ( ~any(keptPartsIdx == ipart) ),
                    keptPartsIdx(end+1) = ipart;
                end
            else
                keptPartsIdx = ipart;
            end
        end
    end % ireg
end % ipart

if Options.verbose > -1,
    fprintf('\n');
end

if input_is_polytope
    % return only Pn if input was not a controller structure
    newCtrlStruct = Pn;
    return
end

newCtrlStruct = Partition{1};

Pfinal = polytope;
for ii = 1:npart,
    Pfinal = [Pfinal Partition{ii}.Pfinal];
end
details.keptParts = keptParts(keptPartsIdx);
newCtrlStruct.Pfinal = Pfinal;
newCtrlStruct.Pn = Pn;
newCtrlStruct.Fi = Fi;
newCtrlStruct.Gi = Gi;
newCtrlStruct.Ai = Ai;
newCtrlStruct.Bi = Bi;
newCtrlStruct.Ci = Ci;
newCtrlStruct.dynamics = dynamics;
newCtrlStruct.details = details;
newCtrlStruct.overlaps = 0;

if mptctrl_input
    % if input was an MPTCTRL object, return MPTCTRL object
    newCtrlStruct = mptctrl(newCtrlStruct);
end

return


% ---------------------------------------------------------------------------------------
function [xcheb, R]=sub_chebyball(H,K,nx,rel_tol,abs_tol,lpsolver)

if any(K<=-1e9),
    % problem seems to be unbounded, perform normalization later on
    exitflag = 0;
else    
    % use 'rescue' function - resolve an LP automatically if it is infeasible
    % with default solver
    [xopt,fval,lambda,exitflag,how]=mpt_solveLPs([zeros(1,nx) 1],[H, -sqrt(sum(H.*H,2))],...
        K,[],[],[],lpsolver);
end

%if ~strcmp(how,'ok')
if exitflag==1,
    % all ok,
else
    % maybe there is a numerical problem, thus we normalize H and K
    
    [nc,nx]=size(H);
    Anorm=sqrt(sum(H .* H,2));
    ii=find(Anorm<=abs_tol | Anorm<=abs(rel_tol*K));
    nii=length(ii);
    if nii>0
        Anorm(ii)=1;
        H(ii,:)=repmat([1 zeros(1,nx-1)],nii,1);
        % decide if some constraint is always true
        jj=(K(ii)>=-abs_tol);
        K(ii(jj))=Inf;
        % or is it always false
        K(ii(~jj))=-Inf;
    end
    temp=1./Anorm;
    H=H .* temp(:,ones(1,nx));
    K=K .* temp;
    
    % If any boundary is -Inf polytope P is empty
    %--------------------------------------------
    if any(K==-Inf)
        xcheb=zeros(nx,1);
        R=-Inf;
        return;
    end
    
    % Remove rows with Inf boundaries
    %--------------------------------
    ii=(K==Inf);
    H(ii,:)=[];
    K(ii)=[];
    
    if size(H,1)==0
        xcheb=zeros(nx,1);
        R=Inf;
        return;
    end
    
    x0 = [zeros(nx,1); 1000];         % hard-coded initial conditions

    % use 'rescue' function - resolve an LP automatically if it is infeasible
    % with default solver
    [xopt,fval,lambda,exitflag,how]=mpt_solveLPs([zeros(1,nx) 1],[H, -sqrt(sum(H.*H,2))],...
        K,[],[],x0,lpsolver);
    
    if exitflag~=1,
        % solution is not optimal, set R to zero, this will be threated as an
        % infeasible transition later on
        xopt = 0*xopt;
    end
end

xcheb = xopt(1:nx); % center of the ball
R=-xopt(nx+1); % Radius of the ball



% ----------------------------------------------------------------------------
function [ctrlStruct] = sub_removeOverlaps_MinTime(ctrlStruct)

lenPn = length(ctrlStruct.Pn);
Pn = polytope; Fi = {}; Gi = {}; Ai = {}; Bi = {}; Ci = {}; dyn = [];
for ii=lenPn:-1:2,
    count = lenPn-ii+1;
    if count==1 | count==lenPn-1 | mod(count,20)==0,
        if count==lenPn-1,
            count=lenPn;
        end
        disp(['Region ' num2str(count) '/' num2str(lenPn)])
    end
    P = ctrlStruct.Pn(ii) \ ctrlStruct.Pn(1:ii-1);
    if isfulldim(P(1)),
        Pn = [Pn P];
        for jj=1:length(P),
            Fi{end+1} = ctrlStruct.Fi{ii};
            Gi{end+1} = ctrlStruct.Gi{ii};
            Ai{end+1} = ctrlStruct.Ai{ii};
            Bi{end+1} = ctrlStruct.Bi{ii};
            Ci{end+1} = ctrlStruct.Ci{ii};
            dyn = [dyn ctrlStruct.dynamics(ii)];
        end
    end
end
Pn = [Pn ctrlStruct.Pn(1)];
dyn = [dyn ctrlStruct.dynamics(1)];
Pn = fliplr(Pn);
dyn = fliplr(dyn);
Fi{end+1} = ctrlStruct.Fi{1};
Gi{end+1} = ctrlStruct.Gi{1};
Ai{end+1} = ctrlStruct.Ai{1};
Bi{end+1} = ctrlStruct.Bi{1};
Ci{end+1} = ctrlStruct.Ci{1};

ctrlStruct.Pn = Pn;
ctrlStruct.Fi = sub_flip_cell(Fi);
ctrlStruct.Gi = sub_flip_cell(Gi);
ctrlStruct.Ai = sub_flip_cell(Ai);
ctrlStruct.Bi = sub_flip_cell(Bi);
ctrlStruct.Ci = sub_flip_cell(Ci);
ctrlStruct.dynamics = dyn;
ctrlStruct.overlaps = 0;

return


%------------------------------------------------------------------------------
function [fC]=sub_flip_cell(C)

lenC = length(C);
fC=cell(1,lenC);
for ii=1:lenC,
    fC{lenC-ii+1}=C{ii};
end


%------------------------------------------------------------------------------
function map = bboxmap(BP, BQ, bbox_tol)

[dim, nP] = size(BP);
[dim, nQ] = size(BQ);
nP = nP/2;
nQ = nQ/2;
map = ones(nP, nQ);
for preg = 1:nP
    for qreg = 1:nQ,
        if any(BP(:,preg*2)+bbox_tol < BQ(:,(qreg-1)*2+1)) 
            % even bounding boxes of the two polytopes do not intersect
            map(preg, qreg) = 0;
        elseif any(BQ(:,qreg*2)+bbox_tol < BP(:, (preg-1)*2+1)),
            % even bounding boxes of the two polytopes do not intersect
            map(preg, qreg) = 0;
        end
    end
end

%------------------------------------------------------------------------------
function [map, P, Q] = imap(P, Q, bbox_tol, maxsph, lpsolver)

bboxOpt.noPolyOutput = 1;
lenP = length(P);
lenQ = length(Q);
map = ones(lenP, lenQ);

PB = cell(1, lenP);
QB = cell(1, lenQ);
if maxsph>0,
    % we will check intersection of V-representations
    PV = cell(1, lenP);
    QV = cell(1, lenQ);
end

% prepare bounding boxes and extreme points
for preg = 1:lenP,
    [d, Pmin, Pmax] = bounding_box(P(preg), bboxOpt);
    if maxsph>0,
        [V, R, d] = extreme(d);
        PV{preg} = V';
    end
    P(preg) = d;
    PB{preg} = [Pmin Pmax];
end
for qreg = 1:lenQ,
    [d, Qmin, Qmax] = bounding_box(Q(qreg), bboxOpt);
    if maxsph>0,
        [V, R, d] = extreme(d);
        QV{qreg} = V';
    end
    Q(qreg) = d;
    QB{qreg} = [Qmin Qmax];
end

% first do pruning based on bounding boxes
for preg = 1:lenP,
    pb = PB{preg};
    for qreg = 1:lenQ,
        qb = QB{qreg};
        if any(pb(:,2)+bbox_tol < qb(:,1)) | any(qb(:,2)+bbox_tol < pb(:, 1)),
            % even bounding boxes of the two polytopes do not intersect
            map(preg, qreg) = 0;
        end
    end
end

if maxsph<=0,
    return
end

maxsph = lenP*lenQ*maxsph;


dim = dimension(P);
vecP = [PV{1:end}];
indiciesP = cumsum([1 cellfun('prodofsize', PV)/dim]);
vecQ = [QV{1:end}];
indiciesQ = cumsum([1 cellfun('prodofsize', QV)/dim]);
sphcount = 0;

for i=1:lenP
    for j = find(map(i,:))
        if map(i,j)
            sphcount = sphcount + 1;
            if sphcount > maxsph,
                % exit if maximum number of allowed separating haperplanes has
                % been reached
                return
            end
            [a,b,f] = separatinghp(PV{i}, QV{j}, lpsolver, dim);
            if (f == 1) % Found one!
                AA = a*vecQ;
                BB = a*vecP;
                jprunes  = zeros(lenQ,1);
                for prunej = 1:lenQ
                    jprunes(prunej) = all((AA(indiciesQ(prunej):indiciesQ(prunej+1)-1)-b)>=bbox_tol);
                end
                jprunes = logical(jprunes);
                if any(jprunes)
                    iprunes  = zeros(lenP,1);
                    for prunei = i:lenP
                        iprunes(prunei) = all((BB(indiciesP(prunei):indiciesP(prunei+1)-1)-b)<=-bbox_tol);
                    end
                    for prunei = find(iprunes)
                        map(prunei,jprunes) = 0;
                    end
                end
            end

        end
    end
end


%-------------------------------------------------------------------
function [a,b,exitflag] = separatinghp(X,Y,lpsolver,nx)

mx = size(X,2);
my = size(Y,2);

A = [X' -ones(mx,1);-Y' ones(my,1)];
b = [-ones(mx,1);-ones(my,1)];
f = zeros(1,nx+1);
[xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f,A,b,[],[],[],lpsolver);
a = xopt(1:nx)';
b = xopt(end);
% how = ~isequal(how,'infeasible');
% how = how & ~all(abs(xopt)<1e-8) & ~any(xopt==1e9);
