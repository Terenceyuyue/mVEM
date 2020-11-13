function model = bb_mpmilp(Matrices,OriginalModel,model,options,ExploreSpace)

nu = Matrices.nu;
nx = Matrices.nx;

if nargin<5
    % Derive a bounding box
    [global_lower,global_upper] = detect_and_improve_bounds(Matrices,Matrices.lb,Matrices.ub,Matrices.binary_var_index,options);
    Matrices.lb = global_lower;
    Matrices.ub = global_upper;
    A = [eye(nx);-eye(nx)];
    b = [global_upper(end-nx+1:end);-global_lower(end-nx+1:end)];
    ExploreSpace = polytope(A,b);
    model = bb_mpmilp(Matrices,OriginalModel,model,options,ExploreSpace);
    return
end

% Find a feasible integer
vartype = repmat('C',nx+nu,1);
vartype(Matrices.binary_var_index) = 'B';
[xmin,fmin,how,exitflag]=mpt_solveMILP([Matrices.H';zeros(nx,1)],[Matrices.G -Matrices.E],Matrices.W,[Matrices.Aeq Matrices.Beq],Matrices.beq,Matrices.lb,Matrices.ub,vartype,[],[],options.mpt.milpsolver);

while isequal(how,'ok')

    feasible_binary = round(xmin(Matrices.binary_var_index));

    % Solve mpLP for this binary
    lower = Matrices.lb;
    upper = Matrices.ub;
    lower(Matrices.binary_var_index) = feasible_binary;
    upper(Matrices.binary_var_index) = feasible_binary;
    nodeModel = solveNode(Matrices,lower,upper,OriginalModel,[],options);

    % Cut away this binary solution
    cut = zeros(1,nu);
    zv = find(feasible_binary==0);
    ov = find(feasible_binary==1);
    cut(Matrices.binary_var_index(ov)) = 1;
    cut(Matrices.binary_var_index(zv)) = -1;
    Matrices.G = [Matrices.G;cut];
    Matrices.E = [Matrices.E;zeros(1,nx)];
    Matrices.W = [Matrices.W;length(ov)-1];
    
    if ~isempty(nodeModel)
        disp('Real case')
        % Dig into this set
        plot(nodeModel{1}.Pn);
        for i = 1:length(nodeModel{1}.Pn)
            newMatrices = Matrices;
            % Add constraint that we explore current region
            [H,K] = double(nodeModel{1}.Pn(i));
            newMatrices.G = [newMatrices.G;zeros(size(H,1),nu)];
            newMatrices.E = [newMatrices.E;-H];
            newMatrices.W = [newMatrices.W;K];
            % Add upper bound constraint HU <
            newMatrices.G = [newMatrices.G;newMatrices.H];
            newMatrices.E = [newMatrices.E;nodeModel{1}.Bi{i}];
            newMatrices.W = [newMatrices.W;nodeModel{1}.Ci{i}];

            % Solve recusively in this new sub region
            plot(nodeModel{1}.Pn(i),'g')
            newmodel{i} = bb_mpmilp(newMatrices,OriginalModel,[],options,nodeModel{1}.Pn(i));
        end
    
        % Now we have explored all feasible regions.
        % As a second step, we must explore infeasible regions     
        InfeasibleRegions = regiondiff(ExploreSpace,nodeModel{1}.Pfinal);
        for i=1:length(InfeasibleRegions)

            newMatrices = Matrices;
            [H,K] = double(InfeasibleRegions(i));
            newMatrices.G = [newMatrices.G;zeros(size(H,1),nu)];
            newMatrices.E = [newMatrices.E;-H];
            newMatrices.W = [newMatrices.W;K];
            
            % Solve recusively in this new sub region
            plot(InfeasibleRegions(i),'g')
            newmodel{end+1} = bb_mpmilp(newMatrices,OriginalModel,model,options,InfeasibleRegions(i));
        end
        
        % Now, we have the upper bound models nodeModel{1} for the investigated
        % binary varible, and candidate models in newmodel

        % Assume nodeModel is best
        model = nodeModel;
        for i = 1:length(newmodel)
            if ~isempty(newmodel{i})
                model = {rmovlps({nodeModel{1},newmodel{i}{1}},struct('verbose',0))};
            end
        end
    end
      
    % Pick a new integer solution
    [xmin,fmin,how,exitflag]=mpt_solveMILP([Matrices.H';zeros(nx,1)],[Matrices.G -Matrices.E],Matrices.W,[Matrices.Aeq Matrices.Beq],Matrices.beq,Matrices.lb,Matrices.ub,vartype,[],[],options.mpt.milpsolver);
end























        
function model = duabb_mpmilp(Matrices,OriginalModel,model,options)

nu = size(Matrices.G,2);
nx = size(Matrices.E,2);
vartype = repmat('C',nx+nu,1);
vartype(Matrices.binary_var_index) = 'B';

% Derive a bounding box 
[global_lower,global_upper] = detect_and_improve_bounds(Matrices,Matrices.lb,Matrices.ub,Matrices.binary_var_index,options);
Matrices.lb = global_lower;
Matrices.ub = global_upper;

% Find a feasible integer
[xmin,fmin,how,exitflag]=mpt_solveMILP([Matrices.H';zeros(nx,1)],[Matrices.G -Matrices.E],Matrices.W,[Matrices.Aeq Matrices.Beq],Matrices.beq,Matrices.lb,Matrices.ub,vartype,[],[],options.mpt.lpsolver);

while isequal(how,'ok')

    feasible_binary = xmin(Matrices.binary_var_index);

    % Solve mpLP for this binary
    lower = Matrices.lb;
    upper = Matrices.ub;
    lower(Matrices.binary_var_index) = feasible_binary;
    upper(Matrices.binary_var_index) = feasible_binary;
    nodeModel = solveNode(Matrices,lower,upper,OriginalModel,model,options);

    if ~isempty(nodeModel)
        % Dig into this set
        plot(nodeModel{1}.Pn);
        if ~isempty(model)
        model = rmovlps({model{1},nodeModel{1}});
        end
        for i = 1:length(nodeModel{1}.Pn)
            newMatrices = Matrices;
            % Add constraint that we explore current region
            [H,K] = double(nodeModel{1}.Pn(i));
            newMatrices.G = [newMatrices.G;zeros(size(H,1),nu)];
            newMatrices.E = [newMatrices.E;-H];
            newMatrices.W = [newMatrices.W;K];
            % Add upper bound constraint HU <
            newMatrices.G = [newMatrices.G;newMatrices.H];
            newMatrices.E = [newMatrices.E;nodeModel{1}.Bi{i}];
            newMatrices.W = [newMatrices.W;nodeModel{1}.Ci{i}];
            % Cut away this binary solution
            cut = zeros(1,nu);
            zv = find(feasible_binary==0);
            ov = find(feasible_binary==1);
            cut(Matrices.binary_var_index(ov)) = 1;
            cut(Matrices.binary_var_index(zv)) = -1;
            newMatrices.G = [newMatrices.G;cut];
            newMatrices.E = [newMatrices.E;zeros(1,nx)];
            newMatrices.W = [newMatrices.W;length(ov)-1];

            % Solve recusively in this new sub region
            plot(nodeModel{1}.Pn(i),'g')
            newmodel = solveDua(newMatrices,OriginalModel,nodeModel,options);
            if ~isempty(newmodel)
                if isempty(model)
                    model = newmodel;
                else
                 model = {rmovlps({model{1},newmodel{1}})};
                end
            end
        end
    end
    % Cut away this binary solution
    cut = zeros(1,nu);
    zv = find(feasible_binary==0);
    ov = find(feasible_binary==1);
    cut(Matrices.binary_var_index(ov)) = 1;
    cut(Matrices.binary_var_index(zv)) = -1;
    Matrices.G = [Matrices.G;cut];
    Matrices.E = [Matrices.E;zeros(1,nx)];
    Matrices.W = [Matrices.W;length(ov)-1];
    % Pick a new integer solution
    [xmin,fmin,how,exitflag]=mpt_solveMILP([Matrices.H';zeros(nx,1)],[Matrices.G -Matrices.E],Matrices.W,[Matrices.Aeq Matrices.Beq],Matrices.beq,Matrices.lb,Matrices.ub,vartype,[],[],options.mpt.lpsolver);
end

% now solve in the complement regions (where this feasible solution is
% infeasible)
%  OtherRegions = regiondiff(polytope([eye(nx);-eye(nx)],[global_upper(end-nx+1:end);-global_lower(end-nx+1:end)]),model{1}.Pfinal);
%  for i = 1:length(OtherRegions)
%
%  end
