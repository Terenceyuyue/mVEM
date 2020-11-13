function Er = rdg(E,af,bf,C,D,b)
%
% Er = rdg(E,af,bf,C,D,b)
%
% Compute the ridges of the facet P_E
%

global zerotol
global verbose

d = size(C,2);
k = size(D,2);
N = size(C,1);

Ec = esp_setdiff([1:N],E);

Ce = C(E,:);
De = D(E,:);
be = b(E,:);
Cec = C(Ec,:);
Dec = D(Ec,:);
bec = b(Ec,:);

% ESP Report section 5.2
[U,S,V] = svd(De);
r = rank(S,zerotol);
s = S(1:r,1:r);
Vh = V(:,1:r);
Vt = V(:,r+1:end);
Uh = U(:,1:r);
Ut = U(:,r+1:end);

iD = Vh*inv(s)*Uh';

S = Cec-Dec*iD*Ce;
U = Dec*Vt;
t = bec-Dec*iD*be;

Er = [];

% Compute dimension of Pe
if(rank([Ce De],zerotol) < k+1)
  % Section 5.2
  %
  %  Note: this can probably be sped up by removing all of the constraints
  %        in Ef since they're all represented by [af bf]. We'd then have
  %        to add then back in after the temporary projection
  
  % Embed the facet in it's affine hull and compute the projection
  [u,s,v] = svd(af');
  s  = s(1);
  vh = v(:,1);
  vt = v(:,2:end);
  Pt = [S*vt U t-S*vh*bf/s];

  if(verbose > 1)
    fprintf('Dimension of Pe is not d-1 - projecting from %i -> %i\n', polydim(Pt), d-1);
  end;
  
  % Project
  [pt,Et] = esp_helper(Pt,[1:d-1],0);
  
  % Convert facet of this projection to the ridges that we're interessted in
  [At,bt] = a2s(pt);
  At = At*vt';
  for i=[1:length(Et)]
    er = sort([E Ec(Et(i).E)]);
    Er = [Er struct('Er',er,'ar',At(i,:)','br',bt(i))];
  end;
else
    % Handle 2D quickly
    if(d == 2)
      naf = esp_null(af')';
%      [minstar,lambda,how,val] = tom_lp(naf,[S t],[af' bf]);
      [minstar,lambda,how,val] = esp_LP(naf,[S t],[af' bf]);
      if(how ~= 1) error('LP error'); end;
      er = find(abs(S*minstar - t) < zerotol);
      [ar,br] = a2s(H2normalize([S(er(1),:) t(er(1))]));
      er = sort([E Ec(er)]);
      Er = [Er struct('Er',er,'ar',ar','br',br)];
      
%      [maxstar,lambda,how,val] = tom_lp(-naf,[S t],[af' bf]);
      [maxstar,lambda,how,val] = esp_LP(-naf,[S t],[af' bf]);
      if(how ~= 1) error('LP error'); end;
      er = find(abs(S*maxstar - t) < zerotol);
      [ar,br] = a2s(H2normalize([S(er(1),:) t(er(1))]));
      er = sort([E Ec(er)]);
      Er = [Er struct('Er',er,'ar',ar','br',br)];
    else
      X = [1:size(S,1)];
      while(~isempty(X))
        i = X(1);
        X = esp_setdiff(X,i);
        
        if(norm(S(i,:)) < zerotol)
          continue;
        end;
        Si = S(i,:)'; 
        Si = Si./norm(Si);
        if(norm(af - (Si'*af)*Si) > zerotol)
          test = esp_null([af' bf;S(i,:) t(i)],1)'*[S t]';
          test = sum(abs(test),1);
          Qc = find(test > zerotol);
          Q  = find(test <= zerotol);
          X  = esp_setdiff(X,Q);
          
          Hi = [S(Qc,:) -ones(length(Qc),1) t(Qc,:)];
          Hi = [Hi;zeros(1,d) -1 1];
          
          He = [af' 0 bf];
          He = [He;S(i,:) 0 t(i)];
          
          He = He(find(sum(abs(He(:,1:end-1)),2) > zerotol),:);
          Hi = Hi(find(sum(abs(Hi(:,1:end-1)),2) > zerotol),:);
          
%          [star,lambda,how,val] = tom_lp([zeros(1,d) 1],Hi,He);
          [star,lambda,how,val] = esp_LP([zeros(1,d) 1],Hi,He);
          if(how ~= 1) error('LP error'); end;
          tau = star(end);

          if(tau < -zerotol)
            er = sort([E Ec(Q)]);
            [ar,br] = a2s(H2normalize([S(i,:) t(i)]));
            Er = [Er struct('Er',er,'ar',ar','br',br)];
          end;
        end;
      end;
    end;
  end;
  
