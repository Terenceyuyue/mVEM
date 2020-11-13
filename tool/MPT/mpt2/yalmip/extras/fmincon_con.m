function [g,geq,dg,dgeq,xevaled] = fmincon_con(x,model,xevaled)

% Early bail for linear problems
g = [];
geq = [];
dg = [];
dgeq = [];
if model.linearconstraints
    xevaled = [];
    return
end

if nargin<3
    xevaled = zeros(1,length(model.c));
    xevaled(model.linearindicies) = x;
    xevaled = apply_recursive_evaluation(model,xevaled);
end

if model.nonlinearinequalities
    g = full(model.Anonlinineq*xevaled(:)-model.bnonlinineq);
end

if model.nonlinearequalities
    geq = full(model.Anonlineq*xevaled(:)-model.bnonlineq);
end

dgAll_test = [];
if 0%all(model.variabletype<=2) & isempty(model.evalMap)
%    Let us set-up a structure for bilinear entities
    is = [];
    js = [];
    gns = [];
    inds = [];
    for i = 1:length(model.variabletype)
        switch model.variabletype(i)
            case 0
                is = [is;i];
                js = [js;i];
                gns = [gns;1];
                inds = [inds;0];
            case 1
                r = find(model.monomtable(i,:));
                is = [is;i];
                js = [js;r(1)];
                gns = [gns;1];
                inds = [inds;r(2)];
                is = [is;i];
                js = [js;r(2)];
                gns = [gns;1];
                inds = [inds;r(1)];
            case 2
                r = find(model.monomtable(i,:));
                is = [is;i];
                js = [js;r(1)];
                gns = [gns;2];
                inds = [inds;r(1)];
        end
    end
    z = [1;xevaled(model.linearindicies)];
    dZ = sparse(is,js,gns.*z(inds+1));
    dgAll_test = [model.Anonlineq;model.Anonlinineq]*dZ;
    dgAll_test = dgAll_test(:,model.linearindicies);
end

if nargout == 2
    return
elseif ~isempty(dgAll_test) & isempty(model.evalMap)
    dgAll = dgAll_test;
elseif isempty(model.evalMap) & (model.nonlinearinequalities | model.nonlinearequalities)
    allA = [model.Anonlineq;model.Anonlinineq];
    dgAll = [];
    n = length(model.c);
    linearindicies = model.linearindicies;
    mtNonlinear = model.monomtable(model.nonlinearindicies,:);
    xevaled = zeros(1,n);
    xevaled(linearindicies) = x;
    X = repmat(xevaled,size(mtNonlinear,1),1);
    % FIXME: This should be vectorized
    
    if 0
        allDerivemt = [];
        news = [];
        c = [];
        for i = 1:length(linearindicies)
            for j = 1:size(mtNonlinear,1)
                if mtNonlinear(j,linearindicies(i))
                    s=mtNonlinear(j,:);
                    c = [c;s(linearindicies(i))];
                    s(linearindicies(i)) = s(linearindicies(i))-1;
                    allDerivemt = [allDerivemt;s(linearindicies(:)')];
                    news = [news;j i];
                end
            end
        end
    end
    
     news = model.fastdiff.news;
     allDerivemt = model.fastdiff.allDerivemt;
     c = model.fastdiff.c;
    
    %
%     
%     dxx = [];
%     for i = 1:length(linearindicies)
%         mt = mtNonlinear;
%         oldpower = mtNonlinear(:,linearindicies(i));
%         mt(:,linearindicies(i)) = mt(:,linearindicies(i))-1;
%         [ii,jj,ss] = find(mt);
%         
%         dxevaledNonLinear = prod(X.^mt,2);
%         dxevaledNonLinear = dxevaledNonLinear(:)'.*oldpower';dxevaledNonLinear(isnan(dxevaledNonLinear))=0;
%         dx = zeros(1,n);
%         dx(linearindicies(i)) = 1;
%         dx(model.nonlinearindicies) = dxevaledNonLinear;
%         dxx = [dxx;dx];
%         dgAll = [dgAll allA*dx'];
%     end
%     
    
    zzz = c.*prod(repmat(x(:)',length(c),1).^allDerivemt,2);
    newdxx = spalloc(length(linearindicies),max(linearindicies),length(linearindicies));
    for i = 1:length(c)
    %    newdxx(news(i,2),model.nonlinearindicies(news(i,1)))=c(i)*prod(X(1,:).^allDerivemt(i,:));
        newdxx(news(i,2),model.nonlinearindicies(news(i,1)))=zzz(i);%c(i)*prod(x(:)'.^allDerivemt(i,:));
    end
    for i = 1:length(linearindicies)
        newdxx(i,linearindicies(i)) = 1;
    end
    dgAll = allA*newdxx';
    %if norm(dgAllnew-dgAll,1)>0
    %    1
    %end
    
    
else
    allA = [model.Anonlineq;model.Anonlinineq];
    requested = any(allA',2);
    [i,j,k] = find((model.deppattern(find(requested),:)));
    requested(j) = 1;
    dx = apply_recursive_differentiation(model,xevaled,requested,model.Crecursivederivativeprecompute);
%     dx2 = apply_oldrecursive_differentiation(model,xevaled,requested);
%     
%     if norm(dx(requested,:)-dx2(requested,:))>1e-14
%         disp('DIFFERENTIATION WRONG')
%         error('DIFFERENTIATION WRONG')
%     end
    
    dgAll = allA*dx;
end

%  if ~isempty(dgAll_test)
%     if norm(dgAll_test-dgAll,inf)>1e-15
%         disp('Hmm')
%     end
%  end

if  model.nonlinearequalities
    dgeq = dgAll(1:size(model.Anonlineq,1),:)';
end
if model.nonlinearinequalities
    dg = dgAll(size(model.Anonlineq,1)+1:end,:)';
end