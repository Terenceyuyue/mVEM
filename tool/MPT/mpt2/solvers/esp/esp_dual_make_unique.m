function E = dual_make_unique(star,C,D,b)
%
% E0min = dual_make_unique(star,C,D,b)
%
% Compute a unique equality set during dual degen
%
% We want to equality set to be unique and for the pre-image of the 
% facet of the projection that contains xstar to be of dimension d-1
%

global zerotol
global DUAL_METHOD

d = size(C,2);
k = size(D,2);
N = size(C,1);

if(strcmp(DUAL_METHOD,'continuous law'))
  % Choose the equality set that is optimal, but minimizes a quadratic cost
  xstar = star(1:d);
  [y,val,flag,output,lambda] = quadprog(eye(k),zeros(k,1),D,b-C*xstar);
  if(flag < 0)
    error('QP error in dual_make_unique');
  end;
  
  E = find(abs(C*xstar+D*y-b) < zerotol);
  E = sort(E);
  E = reshape(E,1,length(E));
elseif(strcmp(DUAL_METHOD,'max region'))
  % Choose the minimal equality set that is optimal

  % Compute the affine hull of the projection
  E0 = find(abs([C D]*star-b) < zerotol);
  [af,bf] = esp_projaff(C(E0,:),D(E0,:),b(E0),1);
  
  % Compute the equality set of the preimage
  hhf = [af' zeros(1,k) bf];
  T = [C D b;hhf;-hhf];
  [Ae,be,Ai,bi,E] = esp_affhull(T);
  E = sort(E(find(E<=N)));
  
  % Check that the equality set defines the same affine hull in the
  % projection
  [at,bt] = esp_projaff(C(E,:),D(E,:),b(E));
  if(norm([af' bf]-[at' bt]) > zerotol)
    error('Affhull returned a bogus eset');
  end;
else
  error('Unknown dual degeneracy handler specified');
end;
