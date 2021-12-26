function tf = mycontains(str, obj)
%%MYCONTAINS the same input and output as contains
%
% See also: strfind
%
% Copyright (C) Terence Yu

matlabversion = version;
location = strfind(matlabversion, 'R');
if str2double(matlabversion(location+1:location+4)) > 2016
    tf = contains(str,obj);
    return; % the remaining code will be neglected.
end

% self code
if ~iscell(str), str = {str}; end
nstr = length(str); tf = true(1,nstr);
for s = 1:nstr
    if isempty(strfind(str{s},obj)), tf(s) = false; end
end