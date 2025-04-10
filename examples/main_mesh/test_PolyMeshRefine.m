clc;clear;close all

load meshdata1.mat
rng(0);
opt_mesh = struct('facecolor','none');
opt_mark = struct('facecolor','r');

toBeMarked = {[7], [11, 19, 20]};
for i = 1 : length(toBeMarked)
    
    marked = toBeMarked{i};
    figure;
    
    % display current mesh and marked elements with node numbering     
    subplot(1,3,1);
    showmesh(node, elem,opt_mesh);
    hold on;
    showmesh(node, elem(marked),opt_mark);
    text(cellfun(@(v) mean(node(v,1)), elem).', cellfun(@(v) mean(node(v,2)), elem).', num2str((1:length(elem)).'), ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Clipping', 'on')
    axis off

    % refine
    [node, elem] = PolyMeshRefine(node, elem, marked);

    % display refined mesh with element numbering
    subplot(1,3,2);
    showmesh(node, elem, opt_mesh);
    text(cellfun(@(v) mean(node(v,1)), elem).', cellfun(@(v) mean(node(v,2)), elem).', num2str((1:length(elem)).'), ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Clipping', 'on')
    axis off;

    % display refined mesh with node numbering
    subplot(1,3,3);
    showmesh(node, elem, opt_mesh);
    text(cellfun(@(v) mean(node(v,1)), elem).', cellfun(@(v) mean(node(v,2)), elem).', num2str((1:length(elem)).'), ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Clipping', 'on')
    text(node(:,1), node(:,2), num2str((1:size(node,1)).'), ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Clipping', 'on');
    axis off;
    
    drawnow;
end