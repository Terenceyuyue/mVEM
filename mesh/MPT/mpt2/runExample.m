%RUNEXAMPLE Demonstrates MPT control routines
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% Allows you to choose a sample dynamical system and to solve a given
% optimization problem.
%
% Note that files describing system dynamics are stored in the
% mpt\examples
% directory. If an error occures, please change your working directory
% to that location, or add mpt\examples to your path.
%
% Available dynamical systems are stored in the following m-files:
%
% Double_Integrator.m - 2nd order LTI system with 1 control input
% ThirdOrder.m        - 3rd order LTI system with 2 control inputs
% FourthOrder         - 4th order LTI system with 2 control inputs
% pwa2d.m             - 2nd order PWA system, 1 control input, 2 PWA dynamics
% pwa_DI.m            - 2nd order PWA system, 1 control input, 4 PWA dynamics
% pwa_sincos.m        - 2nd order PWA system, 1 control input, 2 PWA dynamics
% pwa3d.m             - 3rd order PWA system, 1 control input, 2 PWA dynamics
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% none
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% none
%

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------


global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end
Options=[];


%-------------------------------------------
%GET USER INPUT
%-------------------------------------------
fprintf('\n\n')
disp('CHOOSE SYSTEM')
disp('-------------')
disp('0 - Double Integrator     (2 states, 1 input, LTI)');
disp('1 - 3rd order LTI system  (3 states, 2 inputs, LTI)');
disp('2 - 4th order LTI system  (4 states, 2 inputs, LTI)');
disp('3 - 2nd order PWA System  (2 states, 1 input, 2 dynamics)');
disp('4 - PWA periodical system (2 states, 1 input, 2 dynamics)');
disp('5 - PWA Double Integrator (2 states, 1 input, 4 dynamics)');
disp('6 - 3rd order PWA system  (3 states, 1 input, 2 dynamics)');
fprintf('\n')
p=input('Choose Dynamical System (0,1,2,3,4,5,6): ');

unc=0;   % defualt value - no uncertainty
if p~=2 & p~=6,
    fprintf('\n\n')
    disp('UNCERTAINTY DEFINITION')
    disp('----------------------')
    disp('0 - No uncertainty (assume nominal dynamics)');
    if p~=2    % no additive uncertainty for the 4D LTI example
        disp('1 - Additive uncertainty');
    end
    if p<=1,
        disp('2 - Parametric uncertainty (only for LTI systems)');
    end
    fprintf('\n')
    unc=input('Choose Uncertainty (0,1,2): ');
end

try
    if p==0,
        if unc==0,                   
            Double_Integrator;       % nominal system
        elseif unc==1,               
            Double_Integrator_addU;  % additive uncertainty
        else
            Double_Integrator_parU;  % parametric uncertainty
        end
    elseif p==1,
        if unc==0,                   
            ThirdOrder;              % nominal systems
        elseif unc==1,               
            ThirdOrder_addU;         % additive uncertainty
        else
            ThirdOrder_parU;         % parametric uncertainty
        end
    elseif p==2
        if unc==0,
            FourthOrder;             % nominal system
        else
            error('You need to specify additive or parametric uncertainty manually for this system.');
        end
    elseif p==3,
        if unc==0,
            pwa2d;                   % nominal system
        elseif unc==1
            pwa2d_addU;              % additive uncertainty
        else
            error('Parametric uncertainty not allowed for PWA systems!');
        end
    elseif p==4,
        if unc==0,
            pwa_sincos;              % nominal system
        elseif unc==1
            pwa_sincos_addU;         % additive uncertainty
        else
            error('Parametric uncertainty not allowed for PWA systems!');
        end
    elseif p==5,
        if unc==0,
            pwa_DI;                  % nominal system
        elseif unc==1
            pwa_DI_addU;             % additive uncertainty
        else
            error('Parametric uncertainty not allowed for PWA systems!');
        end
    else
        if unc==0,
            pwa3d;                   % nominal system
        elseif unc==1
            pwa3d_addU;              % additive uncertainty
        else
            error('Parametric uncertainty not allowed for PWA systems!');
        end
    end
catch
    fprintf('\n');
    disp('Error: Examples not available!');
    disp('Please cd to mpt\examples or add this directory to your path and run this file again.')
    return
end

fprintf('\n\n')
disp('CHOOSE OBJECTIVE')
disp('----------------')
if unc==0 | p>=3,
    disp('  1 - Linear cost objective (1-norm)')
    disp('  2 - Quadratic cost objective')
    disp('Inf - Linear cost objective (Inf-norm)');
    fprintf('\n');
    norm=input('Please choose norm in cost function (1/2/Inf): ');
else
    disp('Forcing quadratic cost objective in case of uncertainty');
    norm = 2;
end
   
probStruct.norm=norm;

if mptOptions.qpsolver==-1
    % means that no QP solver is available, switch to linear cost function
    probStruct.norm=1;
    if norm==2,
        fprintf('\n');
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        disp('You did not specify any QP solver, switching to linear cost objective');
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    end
end

fprintf('\n\n')
disp('CHOOSE OPTIMIZATION OBJECTIVE')
disp('-----------------------------')
if p<=2 | norm~=2                   
    if unc==0 & (p>2 | norm==2),  % infinite time solution only available for the nominal case and 1-norms
        disp('0 - Infinite Horizon Cost Optimal Solution')
    end
    if unc==0 | p<=2
        disp('1 - Finite Horizon Solution')
    end
end
if ~(norm==1 & p<=2)    % reduced-complexity algorithm for LTI systems with 1-norm not implemented yet
    disp('2 - Minimum Time Solution')
    disp('3 - Low Complexity Suboptimal Solution')
end
if p<=2,
    % simplex controller can be found for LTI systems
    disp('4 - Simplex controller');
end

fprintf('\n');
b=input('Please choose optimization objective: ');
        
if(iscell(sysStruct.A))
    sysType=1;  %PWA
else
    sysType=0;  %LTI
end

if(b==0)
   disp('Computing Infinite Horizon Solution...')
   probStruct.N=Inf;
   probStruct.subopt_lev=0;
elseif(b==1)
   fprintf('\n')
   probStruct.N=input('Please Specify Prediction Horizon N: ');
   probStruct.subopt_lev=0;
   disp('Computing Finite Horizon Solution...')
elseif(b==2)
   disp('Computing Minimum Time Solution...')
   probStruct.subopt_lev=1; 
elseif b==3,
   disp('Computing Suboptimal Solution...')
   probStruct.subopt_lev=2;
   probStruct.N = 1;
else
    disp('Computing simplex controller...');
end


%-------------------------------------------
%START CONTROLLER COMPUTATION
%-------------------------------------------
fprintf('\n\n')
% compute the control law
if b==4,
    [ctrlStruct] = mpt_simplexContr(sysStruct,probStruct,Options);
else
    [ctrlStruct]=mpt_control(sysStruct,probStruct,Options);
end
 
 
%-------------------------------------------
%PLOT/DISPLAY OUTPUT
%-------------------------------------------
fprintf('\n\n')
Options.newfigure=1;     % by this we tell the plot function to open a new figure window rather then overwrite the current figure

disp('Plotting polyhedral partition...')
% plot the polyhedral partition
mpt_plotPartition(ctrlStruct,Options);
% you can also alternatively use: plotc(ctrlStruct,Options); but you will loose nice coloring

title(['Controller partition with ' num2str(length(ctrlStruct.Pn)) ' regions'])


if p~=1 & p~=2 & p~=6,                        % mpt_plotU can only handle systems with one manipulated variable
    % plot value of control moves
    fprintf('\n');
    mpt_plotU(ctrlStruct,Options);
    fprintf('\n');
    a=input('Do you want to analyze the solution by plotting the trajectories? (1/0): ');
    if a,
        disp('Click on the newly opened figure to pick up the initial state')
        disp('Right-click to abort')
        mpt_plotTrajectory(ctrlStruct,Options);
    end
end
