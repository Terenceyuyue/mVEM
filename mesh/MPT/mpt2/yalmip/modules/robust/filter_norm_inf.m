function [F,feasible] = filter_norm_inf(all_f,all_c_w,all_c_x,all_Q_xw,x,Zmodel,allw,X)

feasible = 1;
lower = Zmodel.lb;
upper = Zmodel.ub;

F = set([]);

% To speed up the construction, compute the ci vectors for all constraints
% in one call ci_basis = [c1 c2 ...]
ci_basis = all_c_w';

collectA = [];
collectb = [];
collectE = [];
collectd = [];
collectC = [];

for i = 1:length(all_f)
    Bi = 2*all_Q_xw(:,length(x)*(i-1)+1:length(x)*i)';
    bi = all_c_x(length(x)*(i-1)+1:length(x)*i);
    
    if (nnz(ci_basis(:,i))==0) & nnz(Bi)==0
        % This constraint row is constant
        F = F + set(X(i)>0);
%        F = F + set(bi(:)'*x + all_f(i) >= 0);
    else
        ci = ci_basis(:,i);

        di = all_f(i);
        % Scale to -1,1 uncertainty
        T = diag(sparse((upper-lower)))/2;
        e = (upper+lower)/2;
        if nnz(Bi) == 0
            if nnz(bi)==0
                % Basically constant + w > 0
                if  (di+e'*ci) - norm(T*ci,1) < 0
                    error('Problem is trivially infeasible');
                    feasible = 0;
                    return
                end
            else
                F = F + set(bi'*x + (di+e'*ci) - norm(T*ci,1) > 0);
            end
        else
            % (bi' + (Bi*w)')*x + (ci'*w + di)
            % (bi' + (Bi*(e+Tw))')*x + (ci'*(e+Tw) + di). |w|<1
            % (bi'+e'*Bi')*x + (di+e'*ci) +(ci'T+x'*Bi*T')*w
            non_zeroBirow = find(sum(abs(Bi'),2));
            zeroBirow = find(sum(abs(Bi'),2) == 0);
            if length(non_zeroBirow)>1
                % This is what we are doing...
                %t = sdpvar(length(non_zeroBirow),1);
                %F = F + set((bi'+e'*Bi')*x + (di+e'*ci) - norm(T(zeroBirow,:)*ci,1)-sum(t) >= 0) + set(-t < T(non_zeroBirow,:)*(ci+Bi'*x) < t);
                % However, to save over-head, we save information and
                % post-pone the generation of a constraint 
                % A*x+b+C*t>0,-t<d+E*x<t
                collectA = [collectA;sparse((bi'+e'*Bi'))];
                collectb = [collectb;(di+e'*ci) - norm(T(zeroBirow,:)*ci,1)];
                collectC = blkdiag(collectC,-sparse(ones(1,length(non_zeroBirow))));
                collectE = [collectE; T(non_zeroBirow,:)*Bi'];
                collectd = [collectd; T(non_zeroBirow,:)*ci];
                %F = F + set((bi'+e'*Bi')*x + (di+e'*ci) - norm(T(zeroBirow,:)*ci,1)-sum(t) >= 0) + set(-t < T(non_zeroBirow,:)*(ci+Bi'*x) < t);
            else
                % There is only one expression involving product between x
                % and w. We explicitly construct the absolut value
                % constraint projection
                F = F + set((bi'+e'*Bi')*x + (di+e'*ci) - norm(T(zeroBirow,:)*ci,1)-T(non_zeroBirow,:)*(ci+Bi'*x) >= 0) ;
                F = F + set((bi'+e'*Bi')*x + (di+e'*ci) - norm(T(zeroBirow,:)*ci,1)+T(non_zeroBirow,:)*(ci+Bi'*x) >= 0) ;
            end
        end
    end
end
if ~isempty(collectA)
    t = sdpvar(size(collectC,2),1);
  %  F = F + set(collectA*x + collectb+collectC*t >= 0);
  %  F = F + set(-t <= collectE*x + collectd <= t);
    F = F + set([-collectA -collectC;collectE -speye(length(t));-collectE -speye(length(t))]*[x;t] + [-collectb;collectd;-collectd] <= 0);
end