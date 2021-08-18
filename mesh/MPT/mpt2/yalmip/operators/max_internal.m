function varargout = max_internal(varargin)

switch class(varargin{1})
    case 'double'
        varargout{1} = max(varargin{:});
    case 'char'
        extstruct.var = varargin{2};
        extstruct.arg = {varargin{3:end}};
        [F,properties,arguments]=max_model([],varargin{1},[],extstruct);
        varargout{1} = F;
        varargout{2} = properties;
        varargout{3} = arguments;
    otherwise
end

function [F,properties,arguments]=max_model(X,method,options,extstruct)
switch method
    case 'graph'
        arguments=[];
        if length(extstruct.arg) == 2-1
            basis = getbase(extstruct.arg{1});
            inf_row = find(basis(:,1) == -inf);
            if length(inf_row)>0
                extstruct.arg{1}(inf_row) = [];
            end
            F = set(extstruct.var - extstruct.arg{1});
            arguments = extstruct.arg{1}(:);
        else
            arguments=[];
            F = set([]);
            for j = 1:(length(extstruct.arg)-1*0)
                basis = getbase(extstruct.arg{j});
                inf_row = find(basis(:,1) == -inf);
                if length(inf_row)>0
                    extstruct.arg{j}(inf_row) = [];
                end
                F = F + set(extstruct.var - extstruct.arg{j});
                arguments = [arguments;extstruct.arg{j}(:)];
            end
        end
        properties = struct('convexity','convex','monotonicity','increasing','definiteness','none');
    case 'exact'
        arguments = [];
        F = set([]);
        t = extstruct.var;
        for j = 1:(length(extstruct.arg)-1*0) % MAX(x,y)
            basis = getbase(extstruct.arg{j});
            inf_row = find(basis(:,1) == -inf);
            if length(inf_row)>0
                extstruct.arg{j}(inf_row) = [];
            end
            X = extstruct.arg{j};
            X = reshape(X,length(X),1);
            if prod(size(X)) == 1
                F = F + set(X == t);
            elseif (prod(size(X)) == 2) & ((nnz(basis(1,:))==0) | (nnz(basis(2,:))==0))
                % Special case to test a particular problem max(0,y), so we
                % keep it since it is optimized
                if (nnz(basis(2,:))==0)
                    X = [0 1;1 0]*X;
                end
                [M,m] = derivebounds(X);
                d = binvar(1,1);
                F = F + set(m(1) <= t      <= m(1) + (M(2)-m(1))*d);
                F = F + set(0    <= t-X(2) <= (m(1)-m(2))*(1-d));
            else
                [M,m] = derivebounds(X);
                n = length(X);
                d = binvar(n,1);
                F = F + set(sum(d)==1);
                F = F + set(-(max(M)-min(m))*(1-d) <= t-X <= (max(M)-min(m))*(1-d));
                kk = [];
                ii = [];
                for i = 1:n
                    k = [1:1:i-1 i+1:1:n]';
                    ii = [ii;repmat(i,n-1,1)];
                    kk = [kk;k];
                    Mm = M(k)-m(i);
                end
                xii = extsubsref(X,ii);
                dii = extsubsref(d,ii);
                xkk = extsubsref(X,kk);
                F = F + set(xkk <= xii+(M(kk)-m(ii)).*(1-dii));
            end
            arguments = [arguments;extstruct.arg{j}(:)];
        end
        properties = struct('convexity','convex','monotonicity','increasing','definiteness','none','model','integer');

    otherwise
        F = [];
        return
end