function [mDIM,nBLOCK,bLOCKsTRUCT,c,F] = sedumi2sdpa(F_struc,c,K);
%SEDUMI2SDPA Internal function to convert SeDuMi structure to format needed in SDPA

% Author Johan Löfberg
% $Id: sedumi2sdpa.m,v 1.2 2004-07-02 08:17:32 johanl Exp $

start = 1;

% This is a hack. K.f is only available when called from bnb with dynamically added equalities
if K.f>0
    F_struc = [-F_struc(1:K.f,:);F_struc];
    K.l = K.l + K.f;
    K.f = 0;
end
if K.s>0
    nc = length(K.s);
else
    nc=0;
end
if K.l>0
    nc = nc+1;
end

F = cell(nc,size(F_struc,2));

% Linear constraints
bLOCKsTRUCT =[];
nl = 0;
if K.l>0
    F{1,1}=sparse(-F_struc(1:K.l,1));
    for j = 2:size(F_struc,2)      
        F{1,j}=sparse(F_struc(1:K.l,j)); 
    end
    start = K.l+1;
    bLOCKsTRUCT = [-K.l];
    nl = 1;
end

% Semidefinite constraints
if K.s>0
    for i = 1:length(K.s)
        theend = start+power(K.s(i),2)-1;
        F{i+nl,1}=triu(sparse(-reshape(F_struc(start:theend,1),K.s(i),K.s(i))));
        for j = 2:size(F_struc,2)      
            F{i+nl,j}=triu(sparse(reshape(F_struc(start:theend,j),K.s(i),K.s(i))));
        end
        start = theend+1;
        bLOCKsTRUCT = [bLOCKsTRUCT K.s(i)];
    end
end
mDIM = size(F_struc,2)-1;
nBLOCK = length(bLOCKsTRUCT);