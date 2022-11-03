function showsolution(node,elem,u,varargin)
%Showsolution displays the solution corresponding to a mesh given by [node,elem] in 2-D.
% Example:
%
%  showsolution(node,elem,u, 'facecolor',[0.5 0.9 0.45], 'FaceAlpha',0.6, 'linewidth',1)
%
% Copyright (C) Terence Yu.

data = [node,u];
if ~iscell(elem)
    h = patch('Faces', elem,...
        'Vertices', data,...
        'FaceColor', 'interp',...
        'CData', u);
else
    max_n_vertices = max(cellfun(@length, elem));
    padding_func = @(vertex_ind) [vertex_ind,...
        NaN(1,max_n_vertices-length(vertex_ind))];  % function to pad the vacancies
    tpad = cellfun(padding_func, elem, 'UniformOutput', false);
    tpad = vertcat(tpad{:});
    h = patch('Faces', tpad,...
        'Vertices', data,...
        'FaceColor', 'interp',...
        'CData', u);
end
axis square; 
sh = 0.015;
xlim([min(node(:,1)) - sh, max(node(:,1)) + sh])
ylim([min(node(:,2)) - sh, max(node(:,2)) + sh])
zlim([min(u) - sh, max(u) + sh])
xlabel('x'); ylabel('y'); zlabel('u');

view(3); grid on; % view(150,30);

if nargin>3
    set(h,varargin{:});
end

