function sH = esp_shift(H,offset)
%
% sH = esp_shift(H,offset)
%
% Shift the polytope P by offset
%
% 13/02/3003 : CNJ
% 01/04/2004 : CNJ modified for use with ESP

if(isempty(H)) 
  sH = [];
else
  [A,B] = a2s(H);
  sH = [A B-A*offset];
end;
