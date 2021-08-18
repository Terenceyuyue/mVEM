function [sol,dgn,Z,J] = mpbbmilp(F,obj,options,parametric)

% Figure out binary variables
binary = recover(intersect(getvariables(F),yalmip('binvariables')));

if isempty(options)
    options = sdpsettings('relax',1);
else
    options = sdpsettings(options,'relax',1);
end

% Solve relaxed problem recursively
[sol,dgn,Z,J] = sub_mpmilp(F,obj,options,parametric,binary,[])

function [sol,dgn,Z,J] = sub_mpmilp(F,obj,options,parametric,binary,fixed)

% Solve relaxed problem, extract only binary solution
disp(['-> Fixed variables : ' num2str(fixed)]);

[sol,dgn,Z,J] = solvemp(F,obj,options,parametric,binary);
if isempty(sol)
    return
end

% Find parametrically integer solutions
not_integer = find((sum(abs([sol{1}.Fi{:}]),2) >1e-8) | sum(abs([sol{1}.Gi{:}]-round([sol{1}.Gi{:}])),2)>1e-8);

if ~isempty(not_integer)

    j = not_integer(1);

    sol_down = sub_mpmilp(F+set(binary(j) == 0),obj,options,parametric,binary,[fixed -j]);
    sol_up   = sub_mpmilp(F+set(binary(j) == 1),obj,options,parametric,binary,[fixed j]);
    
    if isempty(sol_down)
        sol = sol_up;
    elseif isempty(sol_up)
        sol = sol_down;
    else
        sol = {rmovlps({sol_down{1},sol_up{1}})};
    end
else
    disp('-> Integer node found')
end



