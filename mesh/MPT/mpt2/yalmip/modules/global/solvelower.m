function [output,cost,psave] = solvelower(p,options,lowersolver,xmin,upper)

psave = p;
removeThese = find(p.InequalityConstraintState==inf);
p.F_struc(p.K.f + removeThese,:) = [];
p.K.l = p.K.l - length(removeThese);

removeThese = find(p.EqualityConstraintState==inf);
p.F_struc(removeThese,:) = [];
p.K.f = p.K.f - length(removeThese);

p_cut = p;

if ~isempty(p.bilinears)
    p_cut.F_struc(1:p.K.f,:)=[];
    p_cut = addBilinearVariableCuts(p_cut);
    p_cut.F_struc = [p.F_struc(1:p.K.f,:);p_cut.F_struc];
end

if ~isempty(p.evalMap)
    p_cut = addEvalVariableCuts(p_cut);
    psave.evalMap = p_cut.evalMap;
end

% **************************************
% SOLVE NODE PROBLEM
% **************************************

% clf;
% sdpvar x t
% aux = sdpvar(5,1);
% plot(projection(polytope([p_cut.F_struc*[1;x;aux]>0,p_cut.f + p_cut.c'*[x;aux] < t, t<2, p_cut.lb < [x;aux] < p_cut.ub]),[1 2]))
% xx = linspace(p.lb(1),p.ub(1),1000);
% hold on;
% plot(xx,sin(cos(xx.^2).^2).^2 + 0.1*(xx-2).^2)
% axis([-5 5 -0.3 2])
if any(p_cut.ub+1e-8<p_cut.lb)
    output.problem=1;
    cost = inf;
else
    % We are solving relaxed problem (penbmi might be local solver)
    p_cut.monomtable = eye(length(p_cut.c));

    if p.solver.lowersolver.objective.quadratic.convex
        % Setup quadratic
        for i = 1:size(p.bilinears,1)
            if p_cut.c(p.bilinears(i,1))
                p_cut.Q(p.bilinears(i,2),p.bilinears(i,3)) = p_cut.c(p.bilinears(i,1))/2;
                p_cut.Q(p.bilinears(i,3),p.bilinears(i,2)) = p_cut.Q(p.bilinears(i,3),p.bilinears(i,2))+p_cut.c(p.bilinears(i,1))/2;
                p_cut.c(p.bilinears(i,1)) = 0;
            end
        end
        
        if ~all(eig(full(p_cut.Q))>-1e-12)
            p_cut.Q = p.Q;
            p_cut.c = p.c;
        end
    end

    fixed = p_cut.lb >= p_cut.ub;
    if nnz(fixed) == length(p.c)
        % All variables are fixed to a bound
        output.Primal = p.lb;
        res = constraint_residuals(p,output.Primal);
        eq_ok = all(res(1:p.K.f)>=-p.options.bmibnb.eqtol);
        iq_ok = all(res(1+p.K.f:end)>=p.options.bmibnb.pdtol);
        feasible = eq_ok & iq_ok;
        if feasible
            output.problem = 0;
        else
            output.problem = 1;
        end
        cost = output.Primal'*p.Q*output.Primal + p.c'*output.Primal + p.f;
    else

        if nnz(fixed)==0
            output = feval(lowersolver,p_cut);
            cost = output.Primal'*p_cut.Q*output.Primal + p_cut.c'*output.Primal + p.f;
            % Minor clean-up
            pp=p;
            output.Primal(output.Primal<p.lb) = p.lb(output.Primal<p.lb);
            output.Primal(output.Primal>p.ub) = p.ub(output.Primal>p.ub);
            x=output.Primal;
            return
        else
            pp = p_cut;
            removethese = fixed;
            if ~isempty(p_cut.F_struc)
                p_cut.F_struc(:,1)=p_cut.F_struc(:,1)+p_cut.F_struc(:,1+find(fixed))*p_cut.lb(fixed);
                p_cut.F_struc(:,1+find(fixed))=[];

                rf = find(~any(p_cut.F_struc,2));
                rf = rf(rf<=(p_cut.K.f + p_cut.K.l));
                p_cut.F_struc(rf,:) = [];
                p_cut.K.l = p_cut.K.l - nnz(rf>p_cut.K.f);                
                p_cut.K.f = p_cut.K.f - nnz(rf<=p_cut.K.f);               
            end
            p_cut.c(removethese)=[];
            if nnz(p_cut.Q)>0
                p_cut.c = p_cut.c + 2*p_cut.Q(find(~removethese),find(removethese))*p_cut.lb(removethese);
                p_cut.Q(:,find(removethese))=[];
                p_cut.Q(find(removethese),:)=[];
            else
                p_cut.Q = spalloc(length(p_cut.c),length(p_cut.c),0);
            end

            if ~isempty(p_cut.binary_variables)
                new_bin = [];
                new_var = find(~fixed);
                for i = 1:length(p_cut.binary_variables)
                    temp = find(p_cut.binary_variables(i) == new_var);
                    new_bin =  [new_bin temp(:)'];
                end
                p_cut.binary_variables = new_bin;
            end
            if ~isempty(p_cut.integer_variables)
                new_bin = [];
                new_var = find(~fixed);
                for i = 1:length(p_cut.integer_variables)
                    temp = find(p_cut.integer_variables(i) == new_var);
                    new_bin =  [new_bin temp(:)'];
                end
                p_cut.integer_variables = new_bin;
            end

            p_cut.lb(removethese)=[];
            p_cut.ub(removethese)=[];
            p_cut.x0(removethese)=[];
            p_cut.monomtable(:,find(removethese))=[];
            p_cut.monomtable(find(removethese),:)=[];
            try
                output = feval(lowersolver,p_cut);
            catch
                1
            end
            x=p.c*0;
            x(removethese)=p.lb(removethese);
            x(~removethese)=output.Primal;
            output.Primal = x;
            cost = output.Primal'*pp.Q*output.Primal + pp.c'*output.Primal + p.f;
        end      
    end
end