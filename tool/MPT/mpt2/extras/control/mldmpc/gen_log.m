function gen_log(S, log);

%===============================================================================
%
% Title:       gen_log.m                                              
%
% Version:     2.0
%                                                                       
% Project:     Control of MLD systems
%                                                                       
% Purpose:     Generate log-file of all variables
%                                                                        
% Author:      Tobias Geyer
%
% History:     date        subject        
%              2002.02.08  initial version  
%              2002.05.07  displays logical outputs 
%              2002.05.15  display lower and upper bounds on variables
%              2002.08.26  FT: modification in the input parameters
%                              update to MLD struct 2
%              2002.11.27  handles now also the case with X or Y empty
%              2003.11.19  first public relase
%                      
% Input:       S: structure containing the MLD-model
%              log: structure containing the log-data of all variables (d, z, u, x, y)
%                   obtained by runing 'hybsim.m' i.e.: log.x has the dimensions 
%                   nx times nm (nx=number of states, nm=number of time-steps)
%
% Output:      logfile systemname.log on disk
%
% Contact:     Tobias Geyer
%              Automatic Control Laboratory                
%              ETH Zentrum,
%              Zurich, Switzerland
%
%              geyer@control.ee.ethz.ch
%
%              Comments and bug reports are highly appreciated                 
%
%===============================================================================
%
% Legal note:   This program is free software; you can redistribute it and/or
%               modify it under the terms of the GNU General Public License as 
%               published by the Free Software Foundation; either version 2.1 of
%               the License, or (at your option) any later version. 
%
%               This library is distributed in the hope that it will be useful,
%               but WITHOUT ANY WARRANTY; without even the implied warranty of
%               MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%               Lesser General Public License for more details.
% 
%               You should have received a copy of the GNU Lesser General Public
%               License along with this library; if not, write to the 
%               Free Software Foundation, Inc., 
%               59 Temple Place, Suite 330, 
%               Boston, MA  02111-1307  USA
%
%===============================================================================


% adapt this part according to your needs:
% -------------------------------------------------
m_start = 1;		% first time-step to be shown
m_end   = 200;		% last time-step 
m_delta = 20;		% number of time-steps per line
perc = 0.2;         % suggest lower and upper bounds 'perc' lower/higher than the min/max
% -------------------------------------------------


m_end = min( m_end, size(log.D, 2) );

tab_pos = 46;		% column with first data entry

% filename
fname = [S.name '.log'];
fid = fopen(fname, 'wt');
fprintf(fid, ['logfile of ' S.name '.hys\n\n']);

% extract different kinds:
info.X = syminfo(S, '', '', 'x');
info.U = syminfo(S, '', '', 'u');
info.D = syminfo(S, '', '', 'd');
info.Z = syminfo(S, '', '', 'z');
info.Y = syminfo(S, '', '', 'y');
info.P = syminfo(S, '', '', 'p');


% write d, z, u, x and y
% ------------------------------------------------------------------------------
m = m_start;
while (m <= min(length(log.Z(1,:)), m_end)) & (log.Z(1,m) ~= NaN) 
	m0 = m;
	m1 = min(m0+m_delta-1, m_end);
	m1 = min(m1, length(log.Z));
	disp(sprintf('m=%i', m))
	
	str = '  i: t: var.name:       range:            m =';
	min_len = tab_pos;
	for n=m0:m1
		while length(str) < min_len, str = [str ' ']; end;
		min_len = min_len + 7;
		str = [str num2str(n)];
	end;
	fprintf(fid, [str '\n']);

	fprintf(fid, '--------------------------------------------------------------------------------\n');
	fid = write_Rvar(fid, info.D, log.D, 0, tab_pos, m0, m1); fprintf(fid, '\n');
	fid = write_Rvar(fid, info.Z, log.Z, 0, tab_pos, m0, m1); fprintf(fid, '\n');
	fid = write_RBvar(fid,S, info.U, log.U, 0, tab_pos, S.nur, m0, m1); fprintf(fid, '\n');
	if isfield(log, 'X'), fid = write_RBvar(fid,S, info.X, log.X, 0, tab_pos, S.nxr, m0, m1); fprintf(fid, '\n'); end;
	if isfield(log, 'Y'), fid = write_RBvar(fid,S, info.Y, log.Y, 0, tab_pos, S.nyr, m0, m1); fprintf(fid, '\n'); end;
	fprintf(fid, '\n\n');
	
	m = m+m_delta;
end;
	
% write lower and upper bounds of z-variables
% ------------------------------------------------------------------------------
fprintf(fid, 'lower and upper bounds of real variables:\n');
fprintf(fid, 'in simulation: minimal and maximal values of the variable during the simulation\n');
fprintf(fid, 'suggested: bounds in simulation -/+ %2i percent or -/+ 0.1 whatever is small/larger\n', perc*100);
fprintf(fid, '           afterwards, the bounds are rounded\n');
fprintf(fid, 'actual: bounds as currently set in MLD model\n');
fprintf(fid, '\n');
fprintf(fid, '  i  t        var.name   in simulation         suggested       actual \n');
fprintf(fid, '---------------------------------------------------------------------------------\n');
fid = write_bounds(fid, info.Z, log.Z, perc); fprintf(fid, '\n');
fid = write_bounds(fid, info.U, log.U, perc); fprintf(fid, '\n');
if isfield(log, 'X'), fid = write_bounds(fid, info.X, log.X, perc); fprintf(fid, '\n'); end;
if isfield(log, 'Y'), fid = write_bounds(fid, info.Y, log.Y, perc); fprintf(fid, '\n'); end;

fclose(fid);

% convert .log to .ps
% eval(['!a2ps -4 ' S.name '.log -o ' S.name '.ps']);

	
	
	
% ------------------------------------------------------------------------------
% write lower and upper bounds of real variable
function fid = write_bounds(fid, VarInfo, VarLog, perc)
	for i=1:length(VarInfo)
		for j=1:length(VarInfo)
			v = VarInfo{j};
			if (v.index == i) & (v.type == 'r')
                x = VarLog(i,:);
                
                % suggested lower bound
                sug_LB = min(x) - abs(max(x)-min(x))*perc;  
                sug_LB = min( sug_LB, min(x)-0.1 );
                sug_LB = floor(sug_LB*100) / 100;
                
                % suggested upper bound
                sug_UB = max(x) + abs(max(x)-min(x))*perc;  
                sug_UB = max( sug_UB, max(x)+0.1 );
                sug_UB = ceil(sug_UB*100) / 100;
                
                % write the bounds
				str = sprintf('%3i  %1s %15s   %2.5f  %0.5f', v.index, v.kind, v.name, min(x), max(x));
                while length(str) < 47, str = [str ' ']; end;
                str = [str sprintf('%0.2f  %0.2f ', sug_LB, sug_UB)];
                while length(str) < 63, str = [str ' ']; end;
                str = [str sprintf('%0.2f  %0.2f', v.min, v.max)];
                                
                % check for violations
                if min(x) < v.min, str = [str '   error: lower bound violated']; end;
                if max(x) > v.max, str = [str '   error: upper bound violated']; end;
                if sug_LB < v.min, str = [str '   warning: suggested lower bound violated']; end;
                if sug_UB > v.max, str = [str '   warning: suggested upper bound violated']; end;
                
                fprintf(fid, [str '\n']);
            end;
		end;
	end;
return

% write real variable
function fid = write_Rvar(fid, VarInfo, VarLog, index_off, tab_pos, m0, m1)
	for i=1:length(VarInfo)
		for j=1:length(VarInfo)
			v = VarInfo{j};
			if v.index == i
				min_len = tab_pos;
				str = sprintf('%3i  %1s %15s  %2.5f %9.5f', v.index+index_off, v.kind, v.name, v.min, v.max);
				for m=m0:m1
					while length(str) < min_len, str = [str ' ']; end;
					min_len = min_len + 7;
					x = VarLog(i,m);
					if x == round(x)
						str = [str sprintf('%1i', round(x))];
					else
						str = [str sprintf('%1.3f', x)];
					end;
				end;
				fprintf(fid, [str '\n']);
			end;
		end;
	end;
return

% write real and binary variable
function fid = write_RBvar(fid, S, VarInfo, VarLog, index_off, tab_pos, nc, m0, m1)

	% show real variable
	for i=1:length(VarInfo)
		for j=1:length(VarInfo)
			v = VarInfo{j};
			if (v.index == i) & (v.type == 'r')
				min_len = tab_pos;
				str = sprintf('%3i  %1sr%15s  %2.5f %9.5f', v.index+index_off, v.kind, v.name, v.min, v.max);
				for m=m0:m1
					while length(str) < min_len, str = [str ' ']; end;
					min_len = min_len + 7;
					x = VarLog(i,m);
					if x == round(x)
						str = [str sprintf('%1i', round(x))];
					else
						str = [str sprintf('%1.3f', x)];
					end;
				end;
				fprintf(fid, [str '\n']);
			end;
		end;
	end;
	fprintf(fid, '\n');

	% show boolean variable
	for i=1:length(VarInfo)
		for j=1:length(VarInfo)
			v = VarInfo{j};
			if (v.index == i) & (v.type == 'b')
				min_len = tab_pos;
				str = sprintf('%3i  %1sb%15s  %2.5f %9.5f', v.index+index_off, v.kind, v.name, v.min, v.max);
				for m=m0:m1
					while length(str) < min_len, str = [str ' ']; end;
					min_len = min_len + 7;
					x = VarLog(nc+i,m);
					if x == round(x)
						str = [str sprintf('%1i', round(x))];
					else
						str = [str sprintf('%1.3f', x)];
					end;
				end;
				fprintf(fid, [str '\n']);
			end;
		end;
	end;
return	
