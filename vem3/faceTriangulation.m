function [Tri,index3,TriLocal,elemfLocal] = faceTriangulation(elemf)
%Give the triangulation of the face of a polyhedra
% elemf: all the faces
%
% Copyright (C)  Terence Yu.

% number of faces
nf = size(elemf,1);
% numbers of vertices on each face
faceLen = cellfun('length',elemf); 

% triangulation
nTri = sum(faceLen-2);
Tri = zeros(nTri,3);
id = 0;
for s = 1:nf
    face = elemf{s}; face = face(:);
    nt = faceLen(s)-2;
    Tri(id+1:id+nt,:) = [face(1)*ones(nt,1), face(2:end-1), face(3:end)];
    id = id + nt;
end

% global index of vertices
index3 = unique(horzcat(elemf{:}));

% local index of Tri
Idx = min(index3):max(index3);
Idx(index3) = 1:length(index3);
TriLocal = Idx(Tri);

% local index of elemf
elemfLocal = cellfun(@(face) Idx(face), elemf, 'UniformOutput', false); 