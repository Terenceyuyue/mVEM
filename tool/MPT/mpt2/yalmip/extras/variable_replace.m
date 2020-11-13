function Z = variable_replace(X,Y,W)

% Check so that Y is a simple unit variable
Ybase = getbase(Y);
Yvariables = getvariables(Y);
Xbase = getbase(X);
Xvariables = getvariables(X);
[i,j,k] = find(Ybase);
if ~isequal(sort(i),1:length(i))
end
if ~isequal(sort(j),2:(length(i)+1))
end
if ~all(k == 1)
end

[mt,variabletype] = yalmip('monomtable');
if all(variabletype(Xvariables) == 0)
    % Simple linear replacement
    v = 1;
    for i = 1:length(Xvariables)
        XisinY = find(Xvariables(i) == Yvariables);
        if ~isempty(XisinY)
            v = [v;W(XisinY)];
        else
            v = [v;recover(Xvariables(i))];
        end
    end
    Z = Xbase*v;
    Z = reshape(Z,size(X,1),size(X,2));
else
    Z = nonlinearreplace(X,Y,W);
    return
%    error('Nonlinear replacement not supported yet')
end

% This has not been tested (copied from variable_replace) so it is placed
% in a catch to be safe.
try
    Xvariables = getvariables(Z);
    extvar = yalmip('extvariables');
    Xext = find(ismember(Xvariables,extvar));
    if ~isempty(Xext)
        %We must dig down in extended operators to see if they use the replaced
        %set of variables
        for i = 1:length(Xext)
            extstruct = yalmip('extstruct',Xvariables(Xext(i)));
            anychanged = 0;
            for j = 1:length(extstruct.arg)
                if isa(extstruct.arg{j},'sdpvar')
                    XinY = find(ismembc(getvariables(extstruct.arg{j}),Yvariables));
                    if ~isempty(XinY)
                        anychanged = 1;
                        extstruct.arg{j} = replace(extstruct.arg{j},Y,W);
                    else
                    end
                end
            end
            if anychanged
                Zi = yalmip('define',extstruct.fcn,extstruct.arg{:});
                Xvariables(Xext(i)) = getvariables(Zi);
            end
        end
        % And now recover this sucker
        Z = struct(Z);
        Z.lmi_variables = Xvariables;
        % Fucked up order (lmi_variables should be sorted)
        if any(diff(Z.lmi_variables)<0)
            [i,j]=sort(Z.lmi_variables);
            Z.basis = [Z.basis(:,1) Z.basis(:,j+1)];
            Z.lmi_variables = Z.lmi_variables(j);
        end
        Z = sdpvar(Z.dim(1),Z.dim(2),[],Z.lmi_variables,Z.basis);
    end
catch
end
