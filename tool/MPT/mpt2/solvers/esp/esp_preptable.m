function [H,I] = preptable(Er,ar,br,Ef,af,bf)
%
% Create appropriate structure and hash value for ridge data
%

persistent HASH

if(isempty(HASH))
    HASH = ceil(10000*rand(2000,1));
end;

I = struct('lEr',length(Er),'Er',Er,'ar',ar,'br',br,'Ef',Ef,'af',af,'bf',bf);
%H = sum(HASH(Er(1:5)));
H = sum(HASH(sort(Er)));

