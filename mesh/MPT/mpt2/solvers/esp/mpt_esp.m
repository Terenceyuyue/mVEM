function pP = mpt_esp(P,ax)
%
% pP = mpt_esp(P,ax)
%
% Wrapper for the ESP function
%
% Takes and returns MPT polytopes
%

H = double(P);
h = esp(H,ax);
[A,b] = a2s(h);
pP = polytope(A,b);
