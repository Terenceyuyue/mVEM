% load sample mesh
load meshdata1.mat
% for consistent result
rng(0);
% display option
opt_mesh = struct('facecolor','none');
opt_mark = struct('facecolor','r');

nrIter = 8;
timetable = zeros(nrIter, 2);
nrElem = zeros(nrIter, 1);

tic
for i = 1 : nrIter
    marked = randi(length(elem), round([0.8*length(elem), 1]));
    
    % refine
    tic
    [newnode1, newelem1] = PolyMeshRefine(node, elem, marked);
    timetable(i, 1) = toc
    
    tic
    [newnode2, newelem2] = PolyMeshRefine2(node, elem, marked);
    timetable(i, 2) = toc
    
    elem = newelem2;
    node = newnode2;
    
    nrElem(i) = length(newelem2);
    

%     % display refined mesh with element numbering
%     subplot(1,3,2);
%     showmesh(node, elem, opt_mesh);
%     text(cellfun(@(v) mean(node(v,1)), elem).', cellfun(@(v) mean(node(v,2)), elem).', num2str((1:length(elem)).'), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Clipping', 'on')
%     axis off;
% 
%     % display refined mesh with node numbering
%     subplot(1,3,3);
%     showmesh(node, elem, opt_mesh);
%     text(cellfun(@(v) mean(node(v,1)), elem).', cellfun(@(v) mean(node(v,2)), elem).', num2str((1:length(elem)).'), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Clipping', 'on')
%     text(node(:,1), node(:,2), num2str((1:size(node,1)).'), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Clipping', 'on');
%     axis off;
%     drawnow;
end

figure;
loglog(nrElem, timetable(:,1), '-x');
hold on;
loglog(nrElem, timetable(:,2), '--o');
legend('Prev','New');
xlabel('The number of resulting elements');
ylabel('Execution Time');