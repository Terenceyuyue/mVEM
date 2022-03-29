function [Tri,index3,TriLocal,elemfLocal] = faceTriangulation(elemf)
%Give the triangulation of the face of a polyhedra
% elemf: all the faces
%
% Copyright (C)  Terence Yu.

% number of faces
nf = size(elemf,1);
% numbers of vertices on each face
faceLen = cellfun('length',elemf); 
% global index of vertices
[index3,~,totalid] = unique(horzcat(elemf{:})');
index3 = index3(:)';
% local index of elemf
elemfLocal = mat2cell(totalid', 1, faceLen)';
% triangulation
nTri = sum(faceLen-2);
Tri = zeros(nTri,3);  TriLocal = zeros(nTri,3);
id = 0;
for s = 1:nf
    face = elemf{s}; face = face(:);
    faceLocal = elemfLocal{s}; faceLocal = faceLocal(:);
    nt = faceLen(s)-2;
    Tri(id+1:id+nt,:) = [face(1)*ones(nt,1), face(2:end-1), face(3:end)];
    TriLocal(id+1:id+nt,:) = [faceLocal(1)*ones(nt,1), faceLocal(2:end-1), faceLocal(3:end)];
    id = id + nt;
end