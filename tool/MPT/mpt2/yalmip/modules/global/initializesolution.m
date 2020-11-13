function [p,x_min,upper] = initializesolution(p);

x_min = zeros(length(p.c),1);
upper = inf;
if p.options.usex0
    x = p.x0;
    z = evaluate_nonlinear(p,x);
    residual = constraint_residuals(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=p.options.bmibnb.pdtol);
    if relaxed_feasible
        upper = p.f+p.c'*z+z'*p.Q*z;
        x_min = z;
    end
else
    x0 = p.x0;
    p.x0 = zeros(length(p.c),1);
    % Avoid silly warnings
    if ~isempty(p.evalMap)
        for i = 1:length(p.evalMap)
            if (isequal(p.evalMap{i}.fcn,'log') | isequal(p.evalMap{i}.fcn,'log2') | isequal(p.evalMap{i}.fcn,'log10'))
                p.x0(p.evalMap{i}.variableIndex) = (p.lb(p.evalMap{i}.variableIndex) +  p.ub(p.evalMap{i}.variableIndex))/2;
            end
        end
    end
    x = p.x0;
    z = evaluate_nonlinear(p,x);
    residual = constraint_residuals(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=p.options.bmibnb.pdtol);
    if relaxed_feasible
        upper = p.f+p.c'*z+z'*p.Q*z;
        x_min = z;
        x0 = x_min;
    end
    p.x0 = (p.lb + p.ub)/2;
    if ~isempty(p.integer_variables)
        p.x0(p.integer_variables) = round(p.x0(p.integer_variables));
    end
    if ~isempty(p.binary_variables)
        p.x0(p.binary_variables) = round(p.x0(p.binary_variables));
    end
    
    x = p.x0;
    x(isinf(x))=eps;
    x(isnan(x))=eps;
    z = evaluate_nonlinear(p,x);
    residual = constraint_residuals(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=p.options.bmibnb.pdtol);
    if relaxed_feasible & ( p.f+p.c'*z+z'*p.Q*z < upper)
        upper = p.f+p.c'*z+z'*p.Q*z;
        x_min = z;
        x0 = x_min;
    end
    p.x0 = x0;
end

