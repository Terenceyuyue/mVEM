function [A,b] = a2s(H)

if(isempty(H))
  A = [];
  b = [];
else
  A = H(:,1:end-1);
  b = H(:,end);
end;

