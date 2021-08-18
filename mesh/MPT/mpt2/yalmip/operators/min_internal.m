function varargout = min_internal(varargin)

switch class(varargin{1})
    case 'double'
        varargout{1} = min(varargin{:});
    case 'char'
        extstruct.var = varargin{2};
        extstruct.arg = {varargin{3:end}};
        [F,properties,arguments]=min_model([],varargin{1},[],extstruct);
        varargout{1} = F;
        varargout{2} = properties;
        varargout{3} = arguments;
    otherwise
end

function [F,properties,arguments] = min_model(X,method,options,extstruct)
switch method
    case 'graph'
        arguments=[];
        F = set([]);
        for j = 1:(length(extstruct.arg)-1*0)
            basis = getbase(extstruct.arg{j});
            inf_row = find(basis(:,1) == inf);
            if length(inf_row)>0
                extstruct.arg{j}(inf_row) = [];
            end
            F = F + set(extstruct.arg{j} - extstruct.var);
            arguments= [arguments;extstruct.arg{j}(:)];
        end
        properties = struct('convexity','concave','monotonicity','increasing','definiteness','none');
    case 'exact'
        arguments = [];
        F = set([]);
        t = extstruct.var;
        for j = 1:(length(extstruct.arg)-1*0) % MAX(x,y)
            basis = getbase(extstruct.arg{j});
            inf_row = find(basis(:,1) == inf);
            if length(inf_row)>0
                extstruct.arg{j}(inf_row) = [];
            end
            X = extstruct.arg{j};
            X = reshape(X,length(X),1);
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
            F = F + set(xii <= xkk+(M(ii)-m(kk)).*(1-dii));
            arguments = [arguments;extstruct.arg{j}(:)];
        end
        properties = struct('convexity','concave','monotonicity','increasing','definiteness','none','model','integer');

    otherwise
        F = [];
        return
end