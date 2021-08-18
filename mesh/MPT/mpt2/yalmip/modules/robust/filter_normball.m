function [F,feasible] = filter_normball(F_xw,Zmodel,x,w,allw,norm_p)

% Creates robustified version of the uncertain set of linear inequalities
% s.t A(w)*x <= b(w) for all norm(w,p)<r

feasible = 1;
if length(F_xw) == 0
    F = [];
    return
end

X = sdpvar(F_xw);

% Some pre-calc. Consider uncertainties not in this class as decision
% variables (i.e. other than the ones we are robustifying)
vars = [getvariables(x) depends(X) getvariables(allw)];
vars = setdiff(vars,getvariables(w));
x = recover(vars);

%x = recover(unique([getvariables(x),setdiff(getvariables(allw),getvariables(w))]));

% Create a bilinear decomposition of the constraints

xw = [x;w];
xind = find(ismembc(getvariables(xw),getvariables(x)));
wind = find(ismembc(getvariables(xw),getvariables(w)));
[Qs,cs,fs,dummy,nonquadratic] = vecquaddecomp(X,xw);
all_f = [];
all_c_w = [];
all_c_x = [];
all_Q_xw = [];
for i = 1:length(X)
    Q = Qs{i};
    c = cs{i};
    f = fs{i};
    if nonquadratic
        error('Constraints can be at most quadratic, with the linear term uncertain');
    end
    Q_ww = Q(wind,wind);
    Q_xw = Q(xind,wind);
    Q_xx{i} = Q(xind,xind);
    c_x = c(xind);
    c_w = c(wind);

    all_f = [all_f;f];
    all_c_w = [all_c_w;c_w'];
    all_c_x = [all_c_x;sparse(c_x)];
    all_Q_xw = [all_Q_xw Q_xw'];
end
% Linear uncertain constraint is (Bbetai*x + cdi) >= 0 for all w, or
% (bi' + (Bi*w)')*x + (ci'*w + di).

% Currently 3 special cases implemented
switch norm_p
    case 1
        [F,feasible] = filter_norm_1(all_f,all_c_w,all_c_x,all_Q_xw,x,Zmodel,allw,X);

    case 2
        [F,feasible] = filter_norm_2(all_f,all_c_w,all_c_x,all_Q_xw,x,Zmodel,allw,X,Q_xx);

    case inf
        [F,feasible] = filter_norm_inf(all_f,all_c_w,all_c_x,all_Q_xw,x,Zmodel,allw,X);
        
    otherwise
        error('The detected norm-ball has not been implemented yet')
end