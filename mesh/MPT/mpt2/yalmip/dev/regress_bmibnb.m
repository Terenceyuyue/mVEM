function regress_bmibnb


ops = sdpsettings('solver','bmibnb');                 % Global solver
ops = sdpsettings(ops,'bmibnb.lowersolver','glpk');   % Lower solver
ops = sdpsettings(ops,'bmibnb.uppersolver','penbmi'); % Local solver
ops = sdpsettings(ops,'bmibnb.lpsolver','glpk');      % LP solver
ops = sdpsettings(ops,'verbose',0);      % LP solver
ops = sdpsettings(ops,'penbmi.P0',0.01);     
i = 0;

location = fileparts(mfilename('fullpath'));
files = dir(location);
searchfor = [mfilename '_'];
testthese = {};
for i = 1:length(files)
    [cc,vv,ff,gg] = fileparts(files(i).name);
    if isequal(ff,'.m')
        if strfind(files(i).name,searchfor)
            testthese{end+1} = files(i);
        end
    end
end

% Now test them
for i = 1:length(testthese)
    [cc,functionname,vv,gg] = fileparts(testthese{i}.name);
    fail = feval(functionname,ops);
    regressreport([ '''' functionname ''''],fail)
end

function regressreport(text,fail)

switch fail
    case 0
        disp(['No problems in ' text]);
    case 1
        disp(['Objective wrong in ' text]);
    case 2
        disp(['Infeasible solution in ' text]);
    otherwise
end
    
function fail =  getfail(problem,obj,objgoal,infeas)
fail = 0;
if problem == 0
    if (obj-objgoal)>1e-3
        fail = 1;
    elseif max(infeas)<-1e-6
        fail = 2;
    end
else
    fail = 3;
end

