function showmesh(node,elem)
%Showmesh displays a mesh in 2-D and 3-D.

if ~iscell(elem)
    h = patch('Faces', elem, 'Vertices', node);    
else
    max_n_vertices = max(cellfun(@length, elem));
    padding_func = @(vertex_ind) [vertex_ind,...
        NaN(1,max_n_vertices-length(vertex_ind))];  % function to pad the vacancies
    tpad = cellfun(padding_func, elem, 'UniformOutput', false);
    tpad = vertcat(tpad{:});
    h = patch('Faces', tpad, 'Vertices', node);
end

dim = size(node,2);
if dim==3
    view(3); set(h,'FaceAlpha',0.4); % transparency
end

set(h,'facecolor',[0.5 0.9 0.45],'edgecolor','k');
axis equal; axis tight;

