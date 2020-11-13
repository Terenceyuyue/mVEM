function [Ae,be,Ai,bi,E] = esp_affhull(P)
%
% [Ae,be,Ai,bi,E] = esp_affhull(P)
%
% Compute the affine hull of P and the equality set E
%

% Compute the affine hull in a single high-dim LP
[Ae,be,Ai,bi,E] = esp_affhull_one(P);

% Do it one LP low-dim at a time
%[Ae,be,Ai,bi,E] = esp_affhull_checkall(P);
