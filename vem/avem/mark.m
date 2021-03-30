function elemMarked = mark(elem,eta,theta,method)
%Mark marks element.
%
% Copyright (C) Long Chen, modified by Terence Yu.

NT = size(elem,1); isMark = false(NT,1);
if ~exist('method','var'), method = 'Dorfler'; end  
% default marking is Dofler bulk strategy
switch strcmpi(method,'max')
    case true
        isMark(eta>theta*max(eta))=1;
    case false
        [sortedEta,id] = sort(eta.^2,'descend');
        x = cumsum(sortedEta);
        isMark(id(x < theta*x(end))) = 1;
        isMark(id(1)) = 1;
end
elemMarked = find(isMark==true);