function diagnostic = callstrul(F,h,options)

% Author Johan Löfberg
% $Id: callstrul.m,v 1.5 2010-01-13 14:18:15 joloef Exp $

F_new = [];

for i = 1:length(F)
    if ~is(F(i),'lmi')
        F_new = F_new + F(i);
    else
        X = sdpvar(F(i));
        [l,m,r]=factors(X);
        if isempty(m)
            F_new = F_new + F(i);
        else
            [L1,R1,A1,M1,negated_cont1,negated_disc1,epsilon1,delta1,numpos1] = preprocess_constraint(X);
            F_new = F_new + assignschur(F(i),'HKM_schur_LR_structure',L1,R1,A1,M1,negated_cont1,negated_disc1,epsilon1,delta1,numpos1);
        end
    end
end

if nargin < 2
    options = sdpsettings('solver','sdpt3','debug',1,'sdpt3.smallblkdim',1);
else
    options.solver = 'sdpt3';
    options.sdpt3.smallblkdim = 1;
end

diagnostic = solvesdp(F_new,h,options);
