function H=H2normalize(H,rem)
%
% nH = Hnormalize(H,rem)
%
% Given an H-rep, computes an H-rep such that the 2-norm of each column is 1
%
% If rem == 1 remove all zero rows
%
% 27-01-2003 - CNJ - created
% 24-10-2003 - CNJ - added rem
%

global zerotol

if(nargin < 2) rem = 0; end;

%[A,b] = aug2std(H);
A = H(:,1:end-1);
a = sqrt(sum(A.*A,2));
n = size(H,2);
i = find(a>zerotol);
H(i,:) = H(i,:)./a(i,ones(1,n));
%nH = H./repmat(a,1,size(H,2));

if(rem == 1)
    H = H(i,:);
end;
