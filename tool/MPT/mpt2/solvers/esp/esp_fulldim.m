function [Q,E]= espfulldim(P,ax)
%
% ESP helper function
%
%  Compute the projection of a zero containing full dimensional polytope
%

global zerotol
global verbose

% Storage space for the stack
BUCKETS   = spalloc(10^9,1,1000);

% Projection to 1D is a special case
if(length(ax) == 1) 
  [Q,E] = esp_1D(P,ax);
  return;
end;

H = H2normalize(P);
C = H(:,ax);
D = H(:,esp_setdiff([1:polydim(P)],ax));
b = H(:,end);

% Remove zero rows/columns
zerorows  = find(sum(abs([C D]),2) < zerotol);
zeroxcols = find(sum(abs(C),1)     < zerotol);
zeroycols = find(sum(abs(D),1)     < zerotol);

nonzerorows  = esp_setdiff([1:size(C,1)],zerorows);
nonzeroxcols = esp_setdiff([1:size(C,2)],zeroxcols);
nonzeroycols = esp_setdiff([1:size(D,2)],zeroycols);

if(~isempty(zerorows) | ~isempty(zeroxcols) | ~isempty(zeroycols))
  if(verbose > 2)
    fprintf('All zero rows/cols present\n');
  end;
end;

% Drop the zero columns/rows
C = C(nonzerorows,nonzeroxcols);
D = D(nonzerorows,nonzeroycols);
b = b(nonzerorows);

d = length(ax);
k = size(D,2);
N = size(C,1);

% Find initial facet
[E0,af,bf] = esp_shoot(C,D,b);
Er         = esp_rdg(E0,af,bf,C,D,b);

% Push ridges onto the stack
for i=[1:length(Er)]
  [H,t]       = esp_preptable(Er(i).Er,Er(i).ar,Er(i).br,E0,af,bf);
  [S,BUCKETS] = esp_table('add',H,t,BUCKETS);
end;

% Initialize return variables
G = af';
g = bf;
E(1).E = E0;
E(1).a = af';
E(1).b = bf;

% Search until L is empty
count = 0;
while(esp_table('num',[],[],BUCKETS) ~= 0)
  if((verbose >= 1) & (mod(count,20)==0))
    fprintf('Discovered: %4i  To search %4i\r',length(E),esp_table('num',[],[],BUCKETS));
  end;
  count = count + 1;
  
  L = esp_table('get',[],[],BUCKETS);
  
  % Compute the adjacent facet and its ridges
  [Eadj,aadj,badj] = esp_adj(L.Er,L.ar,L.br,L.Ef,L.af,L.bf,C,D,b);
  Er               = esp_rdg(Eadj,aadj,badj,C,D,b);

  % Push the new ridges onto the stack
  % Remove those that were re-discovered
  Found = 0;
  for i=[1:length(Er)]
    [H,t]       = esp_preptable(Er(i).Er,Er(i).ar,Er(i).br,Eadj,aadj,badj);
    [S,BUCKETS] = esp_table('add',H,t,BUCKETS);
    Found = Found + S;
  end;
  if(Found == 0)
    error('Stepped over a ridge that was not rediscovered');
  end;

  % Store the results in return format
  G = [G;aadj'];
  g = [g;badj];
  i = length(E)+1;
  E(i).E = Eadj;
  E(i).a = aadj';
  E(i).b = badj;
end;

% Return the zero rows/cols
d = length(nonzeroxcols) + length(zeroxcols);
for i=[1:length(E)]
  E(i).E  = nonzerorows(E(i).E);
  t = zeros(d,1);
  t(nonzeroxcols) = E(i).a;
  E(i).a = t;
  
  G(i,:) = E(i).a';
  g(i,1)   = E(i).b;
end;

Q = H2normalize([G g]);
