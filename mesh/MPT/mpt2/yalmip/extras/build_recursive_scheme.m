function model = build_recursive_scheme(model);

model.evaluation_scheme = [];
model.monomials = find(model.variabletype);
model.deppattern = model.monomtable | model.monomtable;

if ~isempty(model.evalMap)

    % Figure out arguments in all polynomials & sigmonials
    for i = 1:length(model.monomials)
        model.monomialMap{i}.variableIndex = find(model.monomtable(model.monomials(i),:));
    end

    remainingEvals  = ones(1,length(model.evalVariables));
    remainingMonoms = ones(1,length(model.monomials));
    model = recursive_call(model,remainingEvals,remainingMonoms);

    % Define the dependency structure (used to speed up jacobian
    % computations etc)
    for i = 1:size(model.monomtable,1)
        k = depends_on(model,i);
        model.deppattern(i,k) = 1;
    end
else
    % Only polynomials
    model.evaluation_scheme{1}.group = 'monom';
    model.evaluation_scheme{1}.variables = 1:nnz(model.variabletype);
end

function r = depends_on(model,k)
if model.variabletype(k)
    vars = find(model.monomtable(k,:));
    r=[];
    for i = 1:length(vars)
        r = [r depends_on(model,vars(i))];
    end
elseif ismember(k,model.evalVariables)
    j = find(k == model.evalVariables);
    r = [];
    for i = 1:length(model.evalMap{j}.variableIndex)
        argument = model.evalMap{j}.variableIndex(i);
        % ???? really add argument in all cases
        r = [r argument depends_on(model,argument)];
    end
else
    r = k;
end

function model = recursive_call(model,remainingEvals,remainingMonoms)

if ~any(remainingEvals) & ~any(remainingMonoms)
    return
end

% Yep, this code can be sped up significantly...

% Extract arguments in first layer
if any(remainingEvals)
    for i = 1:length(model.evalMap)
        composite_eval_expression(i) = any(ismember(model.evalMap{i}.variableIndex,model.evalVariables(find(remainingEvals))));
        composite_eval_expression(i) = composite_eval_expression(i) | any(ismember(model.evalMap{i}.variableIndex,model.monomials(find(remainingMonoms))));
    end
end

if any(remainingMonoms)
    for i = 1:length(model.monomials)
        composite_monom_expression(i) = any(ismember(model.monomialMap{i}.variableIndex,model.monomials(find(remainingMonoms))));
        composite_monom_expression(i) = composite_monom_expression(i) | any(ismember(model.monomialMap{i}.variableIndex,model.evalVariables(find(remainingEvals))));
    end
end

% Bottom layer
if ~isempty(model.monomials) & any(remainingMonoms)
    if ~isempty(find(~composite_monom_expression & remainingMonoms))
        model.evaluation_scheme{end+1}.group = 'monom';
        model.evaluation_scheme{end}.variables = find(~composite_monom_expression & remainingMonoms);
    end
    remainingMonoms = composite_monom_expression & remainingMonoms;
end

% Bottom layer
if ~isempty(model.evalMap) & any(remainingEvals)
    if ~isempty(find(~composite_eval_expression & remainingEvals));
        model.evaluation_scheme{end+1}.group = 'eval';
        model.evaluation_scheme{end}.variables = find(~composite_eval_expression & remainingEvals);
    end
    remainingEvals = composite_eval_expression & remainingEvals;
end

model = recursive_call(model,remainingEvals,remainingMonoms);