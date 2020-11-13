function sys = semivar(varargin)
%SEMIVAR Create symbolic semicontinuous variable
%
%   SEMIVAR works exactly as SDPVAR, with the only difference that
%   the elements in the variable automatically will be constrained
%   to be semicontinous (x = 0 or [L < x < U]
%
%   Thwe lower and upper bounds are defined separately.
%
%   See also SDPVAR, BINVAR, INTVAR, BINARY, INTEGER

% Author Johan Löfberg
% $Id: semivar.m,v 1.6 2007-08-02 19:33:16 joloef Exp $

all_variable_names = 1;
i = 1;
while all_variable_names & i<=nargin
    all_variable_names = all_variable_names & isvarname(varargin{i});
    i = i + 1;
end
if all_variable_names
    for k = 1:nargin
        varname = varargin{k};
        assignin('caller',varname,semivar(1,1));
    end
else
    sys = sdpvar(varargin{:}); 
    if isa(sys,'cell')
        for i = 1:length(sys)
            yalmip('setsemicontvariables',[yalmip('semicontvariables') getvariables(sys{i})]);
        end
    else
        yalmip('setsemicontvariables',[yalmip('semicontvariables') getvariables(sys)]);
    end
end


% 
% if isa(varargin{1},'char')
%     switch varargin{1}
%         case {'graph','milp','exact'}
%             t = varargin{2};
%             d = binvar(3,1);
%             [M,m] = derivebounds(t);
%             if m>0
%                 F = set(t >= 1);
%             elseif M<0
%                 F = set(t <= -1);
%             else
%                 F = set(t >= 1+(1-d(1))*(m-1)) + set(t < -1+(M+1)*(1-d(2)));
%                 F = F + set(m.*(1-d(3)) <= t <= M.*(1-d(3)));
%                 F = F + set(sum(d) == 1);
%             end
%             varargout{1} = F;
%             varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
%             varargout{3} = 1;
%             return
%         otherwise
%             % Command line creation
%             all_variable_names = 1;
%             i = 1;
%             while all_variable_names & i<=nargin
%                 all_variable_names = all_variable_names & isvarname(varargin{i});
%                 i = i + 1;
%             end
%             if all_variable_names
%                 for k = 1:nargin
%                     varname = varargin{k};
%                     assignin('caller',varname,semicont(1,1));
%                 end
%                 return
%             end
%     end
% end
% 
% % Create dummy SDPVAR object
% sys = sdpvar(varargin{:});
% n = length(getvariables(sys));
% % Create new variables
% y = [];
% for i = 1:n
%     y = [y;yalmip('define',mfilename,rand)];
% end
% %  Replace variables to the semicontinous variables
% varargout{1} = reshape(getbase(sys)*[1;y],size(sys,1),size(sys,2));