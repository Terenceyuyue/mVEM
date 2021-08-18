function D = esp_setdiff(A,B)
%
% D = esp_setdiff(A,B)
%
% Fast set difference for interger sets
%
%  2004/04/01
%     Colin Jones, Cambridge Control Laboratory, Cambridge UK
%     cnj22@cam.ac.uk

if(isempty(B)) D = A;  return; end;
if(isempty(A)) D = []; return; end;

D = esp_setdiff_c(sort(A),sort(B));
D = D';
return;


% Slow backup
if(isempty(B)) D = A;  return; end;
if(isempty(A)) D = []; return; end;

D = [];
for i=[1:length(A)]
    if(isempty(find(B - A(i) == 0)))
        D = [D A(i)];
    end;
end;

