function dX = apply_recursive_differentiation(model,x,requested,recursivederivativeprecompute);
dX = [];

% Compute all evaluation-based derivatives df(x)
for i = 1:length(model.evaluation_scheme)
    if isequal(model.evaluation_scheme{i}.group,'eval')
        for j = model.evaluation_scheme{i}.variables
            k = model.evalMap{j}.variableIndex;
            if any(requested(model.evalMap{j}.computes))
                z{i,j} = model.evalMap{j}.properties.derivative(x(k));
            end
        end
    end
end

% Apply chain-rule. This code is horrible
if 1
ss = model.deppattern(requested,model.linearindicies);
ss = ss | ss;
ss= sum(ss,1)>1;

for variable = 1:length(model.linearindicies)
    dx = zeros(length(model.c),1);
    dx(model.linearindicies(variable)) = 1;
    if ss(variable)%nnz(model.deppattern(requested,model.linearindicies(variable)))>1
        for i = 1:length(model.evaluation_scheme)
            switch model.evaluation_scheme{i}.group
                case 'eval'
                    for j = recursivederivativeprecompute{variable,i}%model.evaluation_scheme{i}.variables
                        k = model.evalMap{j}.variableIndex;
                        %r = find(dx(k));
                       % if ~isempty(r)                            
                            %if any(requested(model.evalMap{j}.computes))
                                if length(model.evalMap{j}.computes) == 1
                                    dx(model.evalMap{j}.computes) = dx(k)'*z{i,j};
                                else
                                    dx(model.evalMap{j}.computes) = dx(k).*z{i,j};
                                end
                            %end
                       % end
                    end
                case 'monom'
                    computed = model.monomials(model.evaluation_scheme{i}.variables);
                    hh = model.monomtable(computed,:);
%     try
         hh = double(hh | hh)*(dx | dx);
%     catch
%         1
%     end
                    for j = computed(find(hh))
                        if requested(j)
                            dp = 0;
                            monomsj = model.monomtable(j,:);
                            for k = find(dx' & monomsj)
                                monoms = monomsj;
                                monoms(k) = 0;
                                r = model.monomtable(j,k);
                                s = find(monoms);
                                dp = dp + r*x(k)^(r-1)*dx(k)*prod((x(s)').^monoms(s));
                            end
                            dx(j) = real(dp);
                        end
                    end

                otherwise
            end
        end
    end
    dX = [dX dx];
end

































else



% Apply chain-rule. This code is horrible
for variable = 1:length(model.linearindicies)
    dx = zeros(length(model.c),1);
    dx(model.linearindicies(variable)) = 1;
    if nnz(model.deppattern(requested,model.linearindicies(variable)))>1
        for i = 1:length(model.evaluation_scheme)
            switch model.evaluation_scheme{i}.group
                case 'eval'
                    for j = model.evaluation_scheme{i}.variables
                        k = model.evalMap{j}.variableIndex;
                        r = find(dx(k));
                        if ~isempty(r)
                            %any(requested(model.evalMap{j}.computes))
                            if any(requested(model.evalMap{j}.computes))
                                if length(model.evalMap{j}.computes) == 1
                                    dx(model.evalMap{j}.computes) = dx(k)'*z{i,j};
                                else
                                    dx(model.evalMap{j}.computes) = dx(k).*z{i,j};
                                end
                            end
                        end
                    end
                case 'monom'

                    computed = model.monomials(model.evaluation_scheme{i}.variables);
                    hh = model.monomtable(computed,:);
                    hh = double(hh | hh)*(dx | dx);
                    for j = computed(find(hh))
                        if requested(j)
                            dp = 0;
                            monomsj = model.monomtable(j,:);
                            for k = find(dx' & monomsj)
                                monoms = monomsj;
                                monoms(k) = 0;
                                r = model.monomtable(j,k);
                                s = find(monoms);
                                dp = dp + r*x(k)^(r-1)*dx(k)*prod((x(s)').^monoms(s));
                            end
                            dx(j) = real(dp);
                        end
                    end

                otherwise
            end
        end
    end
    dX = [dX dx];
end

end