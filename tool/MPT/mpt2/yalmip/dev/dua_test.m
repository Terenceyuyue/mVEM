function sol = dua_test(F,obj,ops,x)

% F   : All constraints
% obj : Objective
% x   : parametric variables
% y   : all binary variables

if isempty(ops)
    ops = sdpsettings;
end
ops.mp.algorithm = 1;
ops.cachesolvers = 0;
ops.mp.presolve=1;
ops.solver = '';

% Expand nonlinear operators only once 
F = expandmodel(F,obj);
ops.expand = 0;

% Gather all binary variables
y = unique([depends(F) depends(obj)]);
n = length(y)-length(x);
y = intersect(y,[yalmip('binvariables') depends(F(find(is(F,'binary'))))]);
y = recover(y);

% Make sure binary relaxations satisfy 0-1 constraints
F = F + set(0 <= y <= 1);

% Bounded search-space
Universe = unitbox(length(x),10);

% recursive, starting in maximum universe
sol = sub_dua(F,obj,ops,x,y,Universe);

% Nice, however, we have introduced variables along the process, so the
% parametric solutions contain variables we don't care about
for i = 1:length(sol)
    for j = 1:length(sol{i}.Fi)
         sol{i}.Fi{j} = sol{i}.Fi{j}(1:n,:);
         sol{i}.Gi{j} = sol{i}.Gi{j}(1:n,:);
    end
end


function sol = sub_dua(F,obj,ops,x,y,Universe)

sol = {};

% Find a feasible point in this region
% ops.verbose = ops.verbose-1;
% intsol = solvesdp(F,obj,ops);
% ops.verbose = ops.verbose+1;

y_feasible = 0;
if 1
    localsol = {[]};
    intsol.problem = 0;
    while isempty(localsol{1}) & (intsol.problem == 0)
        ops.verbose = ops.verbose-1;
        intsol = solvesdp(F,obj,ops);
        ops.verbose = ops.verbose+1;
        if intsol.problem == 0
            if     isequal(y_feasible , round(double(y)))
                1;
            end
            y_feasible = round(double(y))
            'start'
            localsol = solvemp(F+set(y == y_feasible),obj,sdpsettings(ops,'relax',1),x);
            'end'
            if isempty(localsol{1})
                F = F + not_equal(y,y_feasible);
            end
        end
    end
else
    ops.verbose = ops.verbose-1;
    intsol = solvesdp(F,obj,ops);
    ops.verbose = ops.verbose+1;
    y_feasible = round(double(y));
end

%if intsol.problem == 0
if ~isempty(localsol{1})%intsol.problem == 0
    
    % Compute feasible mpLP solution for this binary combination
    %y_feasible = round(double(y));
    %localsol = solvemp(F+set(y == y_feasible),obj,sdpsettings(ops,'relax',1),x);
    
    % YALMIP syntax...
    if isa(localsol,'cell')
        localsol = localsol{1};
    end
    
    % Now we want to find solutions with other binary combinations, in
    % order to find the best one. Cut away the current bionary using
    % overloaded not equal
    F = F + not_equal(y , y_feasible);
    
    % Could be that the binary was feasible, but the feasible space in the
    % other variables is empty/lower-dimensional
    if ~isempty(localsol)        
        % Dig into this solution. Try to find another feasible binary 
        % combination, with a better cost, in each of the regions      
        for i = 1:length(localsol.Pn)
            G = F;
            % Better cost
            G = G + set(obj <= localsol.Bi{i}*x + localsol.Ci{i});
            % In this region
            [H,K] = double(localsol.Pn(i));
            G = G + set(H*x <= K);
            % Recurse
            diggsol{i} = sub_dua(G,obj,ops,x,y,localsol.Pn(i));
        end

        % Create all neighbour regions, and compute solutions in them too
        flipped = regiondiff(union(Universe),union(localsol.Pn));
        flipsol={};
        for i = 1:length(flipped)
            [H,K] = double(flipped(i));
            flipsol{i} = sub_dua(F+ set(H*x <= K),obj,ops,x,y,flipped(i));
        end

        % Just place all solutions in one big cell. We should do some
        % intersect and compare already here, but I am lazy now.
        sol = appendlists(sol,localsol,diggsol,flipsol);       
    end   
end

function sol = appendlists(sol,localsol,diggsol,flipsol)

sol{end+1} = localsol;
for i = 1:length(diggsol)
    if ~isempty(diggsol{i})
        if isa(diggsol{i},'cell')
            for j = 1:length(diggsol{i})
                sol{end+1} = diggsol{i}{j};
            end
        else
            sol{end+1} = diggsol{i};
        end
    end
end
for i = 1:length(flipsol)
    if ~isempty(flipsol{i})
        if isa(flipsol{i},'cell')
            for j = 1:length(flipsol{i})
                sol{end+1} = flipsol{i}{j};
            end
        else
            sol{end+1} = flipsol{i};
        end
    end
end

