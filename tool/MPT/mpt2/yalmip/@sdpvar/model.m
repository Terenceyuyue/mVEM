function [properties,F,arguments,fcn]=model(X,method,options,extstruct)
%MODEL  Extracts nonlinear operator models
%
% [properties,F] = model(x)
%
% MODEL returns the constraints needed to model a variable related to an
% extended operator such as min, max, abs, norm, geomean, ...
%
% Examples :
%
% sdpvar x y;
% t = min(x,y);
% [properties,F] = epigraph(t)
% Gives (F = set(t<x) + set(t<y))
%
% sdpvar x y
% t = max(norm([x;y],1+y))
% [properties,F] = epigraph(t)
% Gives (F = set(u<t) + set(1+y<t))
% where u is the variable modelling norm([x;y])

% Author Johan Löfberg
% $Id: model.m,v 1.62 2009-05-05 07:20:46 joloef Exp $


extvar = getvariables(X);
arguments   = cell(1,length(extvar));
properties  = cell(1,length(extvar));

if nargin<2
    method = 'graph';
end

if nargin < 3
    options = [];
end

if nargin<4
    extstruct = yalmip('extstruct',extvar);
elseif isempty(extstruct)
    extstruct = yalmip('extstruct',extvar);
end

if isempty(extstruct)
    error('This is not a nonlinear operator variable');
end

fcn = extstruct.fcn;
try
    [F,properties,arguments] = feval(fcn,method,extstruct.var,extstruct.arg{1:end-1});
catch
    error(['Failed when trying to create a model for the "' extstruct.fcn '" operator']);
end

% These fields are not official, or not required yet
if ~isempty(properties)
    if ~iscell(properties)
        properties = {properties};
    end
    for i = 1:length(properties)
        properties{i}.name = fcn;
        if ~isfield(properties{i},'derivative')
            properties{i}.derivative = [];
        end
        if ~isfield(properties{i},'models')
            properties{i}.models = getvariables(extstruct.var);
        end
        if ~isfield(properties{i},'convexhull')
            properties{i}.convexhull = [];
        end
        if ~isfield(properties{i},'bounds')
            properties{i}.bounds = [];
        end
        if ~isfield(properties{i},'domain')
            properties{i}.domain = [-inf inf];
        end
        if ~isfield(properties{i},'range')
            if strcmpi(properties{i}.definiteness,'positive')
                properties{i}.range = [0 inf];
            elseif strcmpi(properties{i}.definiteness,'negative')
                properties{i}.range = [-inf 0];
            else
                properties{i}.range = [-inf inf];
            end
        end
        if ~isfield(properties{i},'model')
            properties{i}.model = 'unspecified';
        end
    end
end

% Normalize the callback expression and check for some obsoleted stuff
if ~isempty(properties)
    if isequal(properties{1}.model,'callback')
        F_normalizing = NormalizeCallback(method,extstruct.var,extstruct.arg{:});
%         if ~isempty(F_normalizing)
%              for i = 1:length(properties)
%                  if ~isinf(properties{i}.domain(1))
%                     F_normalizing = F_normalizing + [extstruct.arg{end} >= properties{i}.domain(1)];
%                  end
%                  if ~isinf(properties{i}.domain(2))
%                     F_normalizing = F_normalizing + [extstruct.arg{end} <= properties{i}.domain(2)];
%                  end                 
%              end
%         end
        F = F + F_normalizing;
    end
    if length(extstruct.computes)>1
        for i = 1:length(properties)
        properties{i}.models = extstruct.computes;
        end
    end
    for i = 1:length(properties)
        if ~any(strcmpi(properties{i}.convexity,{'convex','concave','none'}))
            disp('More cleaning, strange convextiy returned...Report bug')
            error('More cleaning, strange convextiy returned...Report bug')
        end
    end
end

% This is useful in MPT
if ~isempty(F)
    F = tag(F,['Expansion of ' extstruct.fcn]);
end

if ~isempty(properties)
%    properties = properties{1};
end