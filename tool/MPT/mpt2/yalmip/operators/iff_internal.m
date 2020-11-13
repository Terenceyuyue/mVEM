function varargout = iff_internal(varargin)
X = varargin{1};
Y = varargin{2};
switch class(varargin{1})

    case 'sdpvar'

        if length(X)>1
            error('IMPLIES not implemented for this case');
        end

        switch class(Y)
            case 'sdpvar'               % X <--> Y
                varargout{1} = set(Y == X);

            case {'lmi','constraint'}
                Y=set(Y,[],[],1);
                switch settype(Y)
                    case 'elementwise'  % X <--> Y(:)>=0
                        varargout{1} = binary_iff_lp(X,-sdpvar(Y));
                    case 'equality'     % X <--> Y(:)==0
                        varargout{1} = binary_iff_lp(X,sdpvar(Y)) + binary_iff_lp(X,-sdpvar(Y));
                    otherwise
                        error('IFF not implemented for this case');
                end

            otherwise
                error('IFF not implemented for this case');
        end

    case {'lmi','constraint'}

        if isa(X,'constraint')
            X = set(X,[],[],1); % FIX: passes one to avoid pruning infeasible constraints
        end
        switch class(Y)
            case 'sdpvar'
                switch settype(X)
                    case 'elementwise'
                        varargout{1} = binary_iff_lp(Y,-sdpvar(X));
                    case 'equality'
                        varargout{1} = binary_iff_lp(Y,X) + binary_iff_lp(1-Y,-X)
                    otherwise
                        error('IFF not implemented for this case');
                end

            case {'lmi','constraint'} % F(X) <--> F(Y)
                d = binvar(1,1);
                varargout{1} = iff_internal(X,d)+iff_internal(Y,d);

            otherwise
                error('IFF not implemented for this case');
        end

    otherwise
        error('IFF not implemented for this case');
end

function F = binary_iff_lp(X,f)
[M,m,infbound] = derivebounds(f);
if infbound
    warning('You have unbounded variables in IFF leading to a lousy big-M relaxation.');
end
eps = 1e-5;

% X == 1    <=>   f<=0
[nf,mf]=size(f);
if nf*mf==1
    % X == 1 implies f<=0
    F = [f <= M*(1-X)];

    % X == 0 implies f>=0
    F = [F, f >= m*X];

    % f < -eps implies X==1
    F = [F, f >= -eps + (m+eps)*X];

    % f > eps implies X == 0
    F = [F, f <= eps+(M-eps)*(1-X)];

else
    f = reshape(f,nf*mf,1);

    di = binvar(nf*mf,1);
    % di=0 means the ith hypeplane is violated
    % X=1 means we are in the polytope
    F  = set(f <= M*(1-X)) + set(f>=eps+(m-eps).*di)+set(X>=sum(di)-length(di)+1) + set(X <= di);

    % Add some cuts for c < a'x+b < d
    [bA] = getbase(f);
    b = bA(:,1);
    A = bA(:,2:end);
    S = zeros(0,length(di));
    for i = 1:length(b)
        j = findrows(abs(A),abs(A(i,:)));
        j = j(j > i);
        if length(j)==1
            S(end+1,[i j]) = 1;
        end
    end
    if size(S,1) > 0
        % Add cut cannot be outside both constraints
        F = F + set(S*di >= 1);
    end
end
