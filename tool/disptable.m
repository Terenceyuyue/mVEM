function disptable(colname,varargin)

ncol = length(colname);
if length(varargin)~=2*ncol
    error('length(varargin) = 2*ncol !');
end

var = cell(1,ncol);
for s = 1:ncol
    if isempty(varargin{2*s})
        var{s} = num2str(varargin{2*s-1});
        continue;
    end
     var{s} = num2str(varargin{2*s-1}, varargin{2*s});
end

if size(var{1},1)==1  % one-row table
    disp(array2table(str2double(var), 'VariableNames', colname))
else
    disp( table(var{:}, 'VariableNames', colname) );
end