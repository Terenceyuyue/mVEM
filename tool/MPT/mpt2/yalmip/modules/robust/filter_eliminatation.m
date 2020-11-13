function [F_xw,F_robust] = filter_eliminatation(F_xw,w,order)
F_robust = set([]);
Fvars = getvariables(F_xw);
wvars = getvariables(w);
[mt,vt] = yalmip('monomtable');
if any(sum(mt(Fvars,wvars),2)>order)
    for i = 1:length(F_xw)
        Fi = sdpvar(F_xw(i));
        if any(sum(mt(getvariables(Fi),wvars),2)>order)
            [BilinearizeringConstraints,failure] = deriveBilinearizing(Fi,w,order);
            if failure
                error('Cannot get rid of nonlinear uncertainty in uncertain SDP')
            else
                % remove all the violating terms from the expression
                F_xw(i) = clear_poly_dep(F_xw(i),w,order);
                % and add the constraint that actually does this
                F_robust = F_robust + BilinearizeringConstraints;
            end
        end
    end
end
