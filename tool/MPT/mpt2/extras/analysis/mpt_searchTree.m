function ctrlST=mpt_searchTree(ctrl,Options)
%MPT_SEARCHTREE Computes a search-tree for a given explicit solution
%
% ctrlST=mpt_searchTree(ctrl)
% ctrlST=mpt_searchTree(ctrl,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes a search-tree for a given explicit solution according to the algorithm
% proposed by Tondel, Johansen and Bemporad (see literature)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl             - Explicit controller (EPXCTRL object)
% Options.abs_tol  - absolute tolerance
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% ctrlST        - updated controller controller object with the search tree
%                 stored in:
% ctrlST.details.searchTree  - Matrix containing the search tree in the
%                              following form:  
%                       Each row consists of the elements Hi Ki c1 c2 which
%                       have the following meaning:
%         Hi Ki - Definition of the actual hyperplane, where the tree is branched
%         c1 c2 - if positive: Numbers of the rows of searchTree, where the two
%                 child nodes can be found.
%                 if negative: Negative numbers of the control laws, active in 
%                 the two child nodes, which are leaf nodes of the tree
%
% ---------------------------------------------------------------------------
% LITERATURE
% ---------------------------------------------------------------------------
% "Evaluation of piecewise affine control via binary search tree"; P. Tondel, 
% T. A. Johansen and A. Bemporad; Automatica, Vol. 39, No. 5, pp. 945-950
%

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2004 Arne Linder, Faculty of Electrical, Information and Media Engineering, 
%          Wuppertal University, alinder@uni-wuppertal.de

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

if nargin<2,
    Options=[];
end

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if isa(ctrl, 'mptctrl')
    if ~isexplicit(ctrl)
        error('This function supports only explicit controllers!');
    end
    ctrlStruct = struct(ctrl);
else
    ctrlStruct = ctrl;
end

if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;  % absolute tolerance
end
if ~isfield(Options,'rel_tol')
    Options.rel_tol=mptOptions.rel_tol;  % relative tolerance
end
if ~isfield(Options,'lpsolver')
    Options.lpsolver=mptOptions.lpsolver;  % relative tolerance
end

if ~mpt_isValidCS(ctrlStruct),
    error('input argument must be a valid controller structure!');
end

if isfield(ctrlStruct,'searchTree') & isfield(ctrlStruct.details,'searchTree'),
    disp('Search tree is already stored in the controller structure, aborting.');
    ctrlStructST = ctrlStruct;
    return
end
    

PA=ctrlStruct.Pn;
Fi=ctrlStruct.Fi;
Gi=ctrlStruct.Gi;
%Pfinal=hull(ctrlStruct.Pfinal);  % took too long to calculate for PWA systems
Pfinal = bounding_box(ctrlStruct.Pn);

if iscell(ctrlStruct.sysStruct.B),
    nu = size(ctrlStruct.sysStruct.B{1},2);
else
    nu = size(ctrlStruct.sysStruct.B,2);
end

% Preparing work: Search which polytopes have the same control law
% and build a table which polytope corresponds with which control law.
% Control laws which are the same in the first step but different in
% all succeeding steps are regarded as identical. So the principle of
% receding horizon control is used to save on-line calculation time.

nx = dimension(PA);
abs_tol = Options.abs_tol;
rel_tol = Options.rel_tol;
lpsolver = Options.lpsolver;

lenH=length(PA);
Hstore = cell(1, lenH);
Kstore = cell(1, lenH);

fprintf('Eliminating redundant control laws\n');
for i=1:lenH
    [Hstore{i},Kstore{i}] = double(PA(i));
    ii=1;
    diffrow=[];
    for iii=1:nu
        diffrow=[diffrow [Fi{i}(iii,:) Gi{i}(iii)]-[Fi{ii}(iii,:) Gi{ii}(iii)]];
    end
    while (diffrow*diffrow')>Options.abs_tol;
        ii=ii+1;
        diffrow=[];
        for iii=1:nu
            diffrow=[diffrow [Fi{i}(iii,:) Gi{i}(iii)]-[Fi{ii}(iii,:) Gi{ii}(iii)]];
        end
    end
    F(i)=ii;
    if mod(i,10)==0 | i==1 | i==lenH,
        fprintf('Law %d/%d        \r',i,lenH);
    end
end
fprintf('\n');
% first search the j hyperplanes

fprintf('Extracting hyperplanes\n');
rowsH=[];
for i=1:lenH
    %rowsH=[rowsH; double(PA(i))];
    rowsH=[rowsH; [Hstore{i} Kstore{i}]];
end;
[lenT colT]=size(rowsH);


% remove the borders of the set of hyperplanes
rowsH_Pfinal=double(Pfinal);
[lenT_Pfinal colT_Pfinal]=size(rowsH_Pfinal);
i=1;
while i<=lenT_Pfinal
    pivot=0;
    biggest=0.0;
    for ii=1:colT_Pfinal % search biggest element in row to avoid division by zero
        if abs(rowsH_Pfinal(i,ii))>biggest
            biggest=abs(rowsH_Pfinal(i,ii));
            pivot=ii;
        end
    end
    ii=1;
    while ii<=lenT  % search for rows that are borders of the whole region
        factor=rowsH(ii,pivot)/rowsH_Pfinal(i,pivot);
        if size(unique(round([rowsH(ii,:); factor*rowsH_Pfinal(i,:)]*1e12)/1e12,'rows'),1)==1 %row is a border
            rowsH(ii,:)=[]; % delete row
            lenT=lenT-1;
        else
            ii=ii+1; % else test next row
        end 
    end
    i=i+1;
end

fprintf('Removing linear dependent hyperplanes\n');
i=1;
while i<lenT    
    pivot=0;
    biggest=0.0;
    for ii=1:colT   % search biggest element in row to avoid division by zero
        if abs(rowsH(i,ii))>biggest
            biggest=abs(rowsH(i,ii));
            pivot=ii;
        end
    end
    ii=i+1;
    while ii<=lenT  % search for linear dependent rows
        factor=rowsH(ii,pivot)/rowsH(i,pivot);
        if size(unique(round([rowsH(ii,:); factor*rowsH(i,:)]*1e12)/1e12,'rows'),1)==1 % row is linear dependent
            rowsH(ii,:)=[]; % delete row
            lenT=lenT-1;
        else
            ii=ii+1; % else test next row
        end 
    end
    if mod(i,10)==0 | i==1 | i==lenT,
        fprintf('Line %d/%d        \r',i,lenT);    
    end
    i=i+1;
end
fprintf('\n');
% Now rowsH contains only the linear indepentent rows of PA.
% These are the hyperplanes dividing the state space.

% next calculate the index sets for every hyperplane
% notation: I(j+) = I(2*j) // I(j-)=I(2*j-1)

I=cell(2*lenT,1);
for i=1:lenT            % all hyperplanes
    if mod(i,round(lenT/4))==0 | i==1 | i==lenT,
        fprintf('hyperplane %d/%d        \r',i,lenT);
    end
    I{2*i}=[];
    I{2*i-1}=[];
    H2=rowsH(i,1:(colT-1));  
    K2=rowsH(i,colT);   
    for ii=1:lenH       % all polytopes
        H = Hstore{ii};
        K = Kstore{ii};
        Hplus=[H; H2];           % with hyperplane i as additional boundary
        Kplus=[K; K2];           % If the result is fulldimensional then at least a part
        Hminus=[H; -H2];         % of PA(ii) is above respective below the hyperplane.
        Kminus=[K; -K2];
        rplus = sub_chebyball(Hplus,Kplus,nx,rel_tol,abs_tol,lpsolver);
        rminus = sub_chebyball(Hminus,Kminus,nx,rel_tol,abs_tol,lpsolver);
        if rplus>abs_tol  % PA(ii) has elements above hyperplane
            I{2*i}=union(I{2*i},[ii]);
        end 
        if rminus>abs_tol % PA(ii) has elements below hyperplane
            I{2*i-1}=union(I{2*i-1},[ii]);
        end
    end
end
fprintf('\n');
% initialize root node
% node contents: N{i,1}=Ii   N{i,2}=Ji   N{i,3}=ji  N{i,4/5}=child nodes N-/N+
nodecounter=1;
leafnodecounter=0;
N{1,1}=1:lenH;
N{1,2}=[];
U=1; 

% now check the unexplored nodes, until there are no more nodes to explore

ctr = 0;

% catch cases when we have only one region
if length(ctrlStruct.Pn)==1
    U = [];
    searchTree = [zeros(1, dimension(ctrlStruct.Pn)) -1 -1];
end

while ~isempty(U)
    ctr = ctr + 1;
    if mod(ctr,20)==0 | ctr==1
        fprintf('iteration %d        \r',ctr);
    end
    
    Node=U(1);

% Definitions in Tondel's dissertation: 
% Nr = Number of Regions (Polytope)  here: lenH
% K  = Number of affine functions F
% L  = Number of hyperplanes         here: lenT

    U=setdiff(U,Node);

% 5
% compute the approximations intersect(I(Jk),I(j+-)) for all remaining hyperplanes
    Test_set=cell(lenT,2); 
    J_remain=setdiff((1:lenT),fix((N{Node,2}+1)/2));
    for j=J_remain %1:(2*lenT) %Hyperplanes already used in that branch can be neglected
        Test_set{j,1}=intersect(N{Node,1},I{2*j-1}); % minus
        Test_set{j,2}=intersect(N{Node,1},I{2*j});   % plus
    end 
    
% sort the hyperplanes by the quantities of F
    lenT_remain=length(J_remain);   % Number of remaining hyperplanes
    FNum=zeros(lenT_remain,1);
    for j=1:lenT_remain     % determine the sets of F for the intersections
        Fplus=[];
        Fminus=[];
        for i=Test_set{J_remain(j),2} %Test_set{2*j}
            Fplus=union(Fplus,F(i));
        end
        for i=Test_set{J_remain(j),1} %Test_set_{2*j-1}
            Fminus=union(Fminus,F(i));
        end
        FNum(j)=max(length(Fplus),length(Fminus));
    end
    [FNum2,ii]=sort(FNum);
    Sorted=J_remain(ii);
    
% 6
% compute the exact index set of union(Jk,j+-) for the first nj elements of the sorted list.
% nj is chosen to include all hyperplanes, which minimize the quantites of F

    % determine nj
    nj=1;
    ii=FNum2(1);
    while (nj<lenT_remain) & (FNum2(nj+1)==ii)
        nj=nj+1;
    end
    
    % compute exact index sets
    [H K]=double(Pfinal);
    for j=N{Node,2}        % first calculate polytope of Jk
        i=fix((j+1)/2);
        H2=rowsH(i,1:(colT-1));
        K2=rowsH(i,colT);
        if (2*i>j)
            H2=-H2;
            K2=-K2;
        end
        H = [H; H2];
        K = [K; K2];
    end
    Iplus=cell(nj,1);
    Iminus=cell(nj,1);
    for j=1:nj
        H2=rowsH(Sorted(j),1:(colT-1));  % now, calculate polytope of Jk and j+-
        K2=rowsH(Sorted(j),colT);
        Hplus=[H; H2];
        Kplus=[K; K2];
        Hminus=[H; -H2];
        Kminus=[K; -K2];
        
        rplus = sub_chebyball(Hplus,Kplus,nx,rel_tol,abs_tol,lpsolver);
        rminus = sub_chebyball(Hminus,Kminus,nx,rel_tol,abs_tol,lpsolver);

        Opt_int.reduce_intersection=0;
        if rplus>abs_tol    % Pplus is fulldim
            for i=N{Node,1} % faster since only polytops not cut away by earlier hyperplanes have to be checked.
                Hi = Hstore{i};
                Ki = Kstore{i};
                rint = sub_chebyball([Hplus;Hi],[Kplus;Ki],nx,rel_tol,abs_tol,lpsolver);
                fulldim = rint > abs_tol;
                if fulldim  % if yes, add number to index list
                    Iplus{j}=union(Iplus{j},i);
                end
            end
        end
        if rminus>abs_tol    % Pminus is fulldim
            for i=N{Node,1}
                Hi = Hstore{i};
                Ki = Kstore{i};
                rint = sub_chebyball([Hminus;Hi],[Kminus;Ki],nx,rel_tol,abs_tol,lpsolver);
                fulldim = rint > abs_tol;
                if fulldim
                    Iminus{j}=union(Iminus{j},i);
                end
            end
        end
    end
    
    % now, search again for minimum of max(F(I+),F(i-)) 
    Fplus=cell(nj,1);
    Fminus=cell(nj,1);
    ii=Inf;        % max possible value
    for j=1:nj     % determine the sets of F for the new index sets
        for i=Iplus{j}
            Fplus{j}=union(Fplus{j},F(i));
        end
        for i=Iminus{j}
            Fminus{j}=union(Fminus{j},F(i));
        end
        if max(length(Fplus{j}),length(Fminus{j}))<ii   % new mininum found
            ii=max(length(Fplus{j}),length(Fminus{j}));
        end
    end
    minimum=[];    % set of all hyperplanes with max(F(I+),F(I-))=min.
    for j=1:nj
        if max(length(Fplus{j}),length(Fminus{j}))==ii
            minimum=union(minimum,j);
        end
    end
    if length(minimum)>1    % more than 1 hyperplane are optimal ==> additional criteria
    ii=Inf;
        for j=minimum
            if max(length(Iplus{j}),length(Iminus{j}))<ii  % new minimum found
                ii=max(length(Iplus{j}),length(Iminus{j}));
                optimal=j;        % if still more than 1 optimal solution ==> don't care
            end
        end
    else
        optimal=minimum(1);
    end
    if isempty(Fplus{optimal})
        Fplus{optimal} = 0;
    end
    if isempty(Fminus{optimal})
        Fminus{optimal} = 0;
    end
    if isempty(Iplus{optimal})
        Iplus{optimal} = 0;
    end
    if isempty(Iminus{optimal})
        Iminus{optimal} = 0;
    end
    
    
% 7 complete the node
    N{Node,3}=Sorted(optimal);
    for i=1:colT
        searchTree(Node,i)=rowsH(Sorted(optimal),i);
    end
    % create N- as child node
    N{Node,4}=nodecounter+1;
    N{nodecounter+1,1}=Iminus{optimal};
    N{nodecounter+1,2}=union(N{Node,2},Sorted(optimal)*2-1);
    % create N+ as child node
    N{Node,5}=nodecounter+2;
    N{nodecounter+2,1}=Iplus{optimal};
    N{nodecounter+2,2}=union(N{Node,2},Sorted(optimal)*2);
    
    % 8 check, if new nodes are leaf nodes
    if length(Fminus{optimal})>1
        U=union(U,nodecounter+1);
        searchTree(Node,colT+1)=nodecounter+1-leafnodecounter;
    else
        N{nodecounter+1,3}=Fminus{optimal};
        N{nodecounter+1,4}=0;
        N{nodecounter+1,5}=0;
        searchTree(Node,colT+1)=-Fminus{optimal};
        leafnodecounter=leafnodecounter+1;
    end
    if length(Fplus{optimal})>1
        U=union(U,nodecounter+2);
        searchTree(Node,colT+2)=nodecounter+2-leafnodecounter;
    else
        N{nodecounter+2,3}=Fplus{optimal};
        N{nodecounter+2,4}=0;
        N{nodecounter+2,5}=0;
        searchTree(Node,colT+2)=-Fplus{optimal};
        leafnodecounter=leafnodecounter+1;
    end

% 9
    % increment nodecounter
    nodecounter=nodecounter+2;
end
fprintf('Algorithm finished at iteration %d        \r\n',ctr);

% Remove the empty lines from searchTree
[lenST colST]=size(searchTree);
i=1;
while (i<=lenST)
    if (searchTree(i,colST-1)==0) & (searchTree(i,colST)==0)
        searchTree=[searchTree(1:(i-1),:);searchTree((i+1):lenST,:)];
        lenST=lenST-1;
    else
        i=i+1;
    end
end

ctrlStructST = ctrlStruct;
ctrlStructST.details.searchTree = searchTree;

ctrlST = mptctrl(ctrlStructST);

% assign variable which contains the search tree in caller's workspace
if ~isempty(inputname(1)) & nargout==0,
    assignin('caller',inputname(1),ctrlST);
end

return



% ---------------------------------------------------------------------------------------
function R=sub_chebyball(H,K,nx,rel_tol,abs_tol,lpsolver)

if all(K>-1e9),
    % use 'rescue' function - resolve an LP automatically if it is infeasible
    % with default solver
    [xopt,fval,lambda,exitflag,how]=mpt_solveLPs([zeros(1,nx) 1],[H, -sqrt(sum(H.*H,2))],...
        K,[],[],[],lpsolver);
else
    how = 'infeasible';
end

if ~strcmp(how,'ok')
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
    %S%H=H .* repmat(temp,1,nx);
    H=H .* temp(:,ones(1,nx));
    K=K .* temp;
    
    % If any boundary is -Inf polytope P is empty
    %--------------------------------------------
    if any(K==-Inf)
        xc=zeros(nx,1);
        R=-Inf;
        lambda=zeros(nc,1);
        return;
    end
    
    % Remove rows with Inf boundaries
    %--------------------------------
    ii=(K==Inf);
    H(ii,:)=[];
    K(ii)=[];
    
    if size(H,1)==0
        xc=zeros(nx,1);
        R=Inf;
        lambda=zeros(nc,1);
        return;
    end
    
    x0 = [zeros(nx,1); 1000];         % hard-coded initial conditions
    
    % use 'rescue' function - resolve an LP automatically if it is infeasible
    % with default solver
    [xopt,fval,lambda,exitflag,how]=mpt_solveLPs([zeros(1,nx) 1],[H, -sqrt(sum(H.*H,2))],K,[],[],x0,lpsolver);
end

R=-xopt(nx+1); % Radius of the ball
