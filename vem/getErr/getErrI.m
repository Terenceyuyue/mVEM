function ErrI = getErrI(node,elem,uh,info,pde,kOrder)

chie = getDof(node,elem,pde,kOrder);
chi = uh;
kk = info.kk; DofI = info.DofI;
chie = chie(DofI); chi = chi(DofI);
kk = kk(DofI,DofI);
ErrI = sqrt(abs((chi-chie)'*kk*(chi-chie)));



