function [F,mptmissing] = filter_enumeration(F_xw,Zmodel,x,w,ops)

mptmissing = 0;
if length(F_xw) == 0
    F = [];
    return;
else

    if any(Zmodel.K.q) | any(Zmodel.K.s)
        error('Only polytope uncertainty supported in duality based robustification');
    else
        if isempty(intersect(depends(F_xw),getvariables(w)))
            F = F_xw;
        else
            % FIX : Assumes all uncertainty in all constraints
            K = Zmodel.K;
            A = -Zmodel.F_struc((1+K.f):(K.f + K.l),2:end);
            b =  Zmodel.F_struc((1+K.f):(K.f + K.l),1);

            try
                % Some preprocessing to extract bounds from equality
                % constraints in order to make the uncertainty polytope
                % bounded (required since we are going to run vertex
                % enumeration)
                % We might have x>=0, sum(x)=1, and this code simply extracts
                % the implied bounds x<=1
                [lo,up] = findulb(Zmodel.F_struc(1:K.f + K.l,:),K);
                Zmodel.lb = lo;Zmodel.ub = up;
                Zmodel = presolve_bounds_from_equalities(Zmodel);
                up = Zmodel.ub;
                lo = Zmodel.lb;
                upfi = find(~isinf(up));
                lofi = find(~isinf(lo));
                aux = Zmodel;
                aux.F_struc = [aux.F_struc;-lo(lofi) sparse(1:length(lofi),lofi,1,length(lofi),size(A,2))];
                aux.F_struc = [aux.F_struc;up(upfi) -sparse(1:length(upfi),upfi,1,length(upfi),size(A,2))] ;
                aux.K.l = aux.K.l + length(lofi) + length(upfi);
                K = aux.K;
                A = -aux.F_struc((1+K.f):(K.f + K.l),2:end);
                b =  aux.F_struc((1+K.f):(K.f + K.l),1);
                P = polytope(A,b);
                if ~isbounded(P)
                    error('The uncertainty space is unbounded (could be an artefact of YALMIPs modelling of nonolinear oeprators).')                     
                else
                    vertices = extreme(polytope(A,b))';
                end
            catch
                mptmissing = 1;
                if ops.verbose>0
                    %lasterr
                    disp('You probably need to install MPT (needed for vertex enumeration)')
                    disp('http://control.ee.ethz.ch/~joloef/wiki/pmwiki.php?n=Solvers.MPT')
                    disp('Alternatively, you need to add bounds on your uncertainty.')
                    disp('Trying to switch to dualization approach')
                    %error('MPT missing');
                end
                F = [];
                return
            end
            % The vertex enumeration was done without any equality constraints.
            % We know check all vertices so see if they satisfy equalities.
            if K.f > 0
                Aeq = -Zmodel.F_struc(1:K.f,2:end);
                beq =  Zmodel.F_struc(1:K.f,1);
                feasible = sum(abs(Aeq*vertices - repmat(beq,1,size(vertices,2))),1) < 1e-6;
                vertices = vertices(:,feasible);
                if isempty(feasible)
                    error('The uncertainty space is infeasible.')
                end
            end

            % We know replace all occurances of w with the fixed vertices
            % Doing LP constraints in a vectorized manner saves a lot of time
            F_xw_lp = F_xw(find(is(F_xw,'elementwise')));
            F_xw_socp_sdp = F_xw -  F_xw_lp;
            F = set([]);
            if length(F_xw_lp)>0
                rLP = [];
                for i = 1:size(vertices,2)
                    rLP = [rLP;replace(sdpvar(F_xw_lp),w,vertices(:,i),0)];
                end

                % FIXME: More general detection of silly constraints
                if isa(rLP,'double') & all(rLP>=-eps^0.75)
                    F = set([]);
                else
                    % Easily generates redundant constraints
                    [aux,index] = uniquesafe(getbase(rLP),'rows');
                    try
                        F = set(rLP(index(randperm(length(index)))) >= 0);
                    catch
                        1
                    end
                end
            end

            % Remaining conic stuff
            for j = 1:length(F_xw_socp_sdp)
                for i = 1:size(vertices,2)
                    F = F + set(replace(F_xw_socp_sdp(j),w,vertices(:,i),0));
                end
            end
        end
    end
end
