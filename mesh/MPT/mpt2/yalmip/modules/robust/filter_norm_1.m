function [F,feasible] = filter_norm_1(all_f,all_c_w,all_c_x,all_Q_xw,x,Zmodel,allw,X)

feasible = 1;
% As a first step, we figure out the radius
r = Zmodel.r;
center = Zmodel.center;
F = set([]);
ci_basis = all_c_w';
lastBici = [];
for i = 1:length(all_f)
    Bi = 2*all_Q_xw(:,length(x)*(i-1)+1:length(x)*i)';
    bi = all_c_x(length(x)*(i-1)+1:length(x)*i);
    if (nnz(ci_basis(:,i))==0) & nnz(Bi)==0
        F = F + set(X(i)>0);
        % This constraint row is constant
    %    ci = ci_basis(:,i);
    %    di = all_f(i);
    %    used = full(any([full(Bi') ci],2));
    %    ci = ci(used);
    %    Bi = Bi(:,used);
        % Shift |w-center|, wtilde = w-center i.e. w=wtilde+center
    %    di = di + ci'*center;
    %    bi = bi + Bi*center;
    %    F = F + set(bi(:)'*x + all_f(i) >= 0);
    else
        ci = ci_basis(:,i);
        di = all_f(i);
        used = find(full(any([full(Bi') ci],2)));
        ci = ci(used);
        Bi = Bi(:,used);
        % Shift |w-center|, wtilde = w-center i.e. w=wtilde+center
        di = di + ci'*center;
        bi = bi + Bi*center;
        if isequal(lastBici,[Bi' ci])
            F = F + set(bi'*x + di - r*s > 0);% + set(-s<Bi'*x+ci<s);
        else
            s = sdpvar(1,1);
            F = F + set(bi'*x + di - r*s > 0) + set(-s<Bi'*x+ci<s);
            lastBici = [Bi' ci];
        end
    end
end