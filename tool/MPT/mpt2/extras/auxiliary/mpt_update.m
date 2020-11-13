function mpt_update
%MPT_UPDATE Checks for updates of the Multi-Parametric Toolbox
%
% verstr = mpt_version
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Connects to http://control.ee.ethz.ch/~mpt/ and checks for latest available
% version of MPT. Notifies the user if new updates are available.
%
% No data are send to the server.
%
% Note: Requires Matlab 6.5 or newer
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
%
%
% see also MPT_INIT, MPT_VERSION

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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


% this check will only succeed with Matlab 6.5 and Matlab 7
if exist('urlread','file'),
    
    % get version of local MPT
    ver_local = mpt_init('version');
    % try to connect to get update info
    try
        % fetch info file
        [info,status] = urlread('http://control.ee.ethz.ch/~mpt/downloads/update.nfo');
        if status == 0,
            disp('mpt_update: cannot connect to update site, check your internet connection');
            return
        end
        
        % strip version information
        [ver, rest] = strtok(info, '|');

        % strip release date
        [reldate, rest] = strtok(rest, '|');
        
        % strip download location
        [dload_loc, rest] = strtok(rest, '|');
        
        % parse version strings
        [major, minor, revision] = sub_parseVersion(ver);
        [major_loc, minor_loc, rev_loc] = sub_parseVersion(ver_local);
        
        newupdate = 0;
        
        % check version number
        
        if major > major_loc,
            newupdate = 1;
        elseif (major == major_loc) & (minor > minor_loc),
            newupdate = 1;
        elseif (major == major_loc) & (minor == minor_loc),
            if ischar(rev_loc) & ischar(revision),
            elseif ischar(rev_loc) & ~ischar(revision),
                newupdate = 1;
            elseif ~ischar(rev_loc) & ~ischar(revision),
                if revision > rev_loc,
                    newupdate = 1;
                end
            end
        end
        
        if newupdate
            fprintf('New updates are available...\n');
            fprintf('Local version: %s\n', ver_local');
            fprintf('Version available for download: %s\n', ver);
            fprintf('Release date: %s\n', reldate);
            fprintf('Download location: %s\n', dload_loc);
        else
            fprintf('Your installation of MPT is up to date.\n');
        end
        fprintf('\n');
        return
    catch
        disp('Problems connecting to MPT site...');
        disp('Please check for updates manually at http://control.ee.ethz.ch/~mpt/');
        fprintf('\n');
        return
    end
else
    disp('New version notification requires Matlab 6.5 or newer...');
    disp('Please check for updates manually at http://control.ee.ethz.ch/~mpt/');
    fprintf('\n');
end
    

function [major, minor, revision] = sub_parseVersion(verstr),

% find position of dots in string
dotpos = strfind(verstr, '.');

if length(dotpos)==1,
    % version number in form Major.Minor
    ver_num = sscanf(verstr, '%d.%d');
    if length(verstr)>3,
        revision = verstr(4:end);
        if double(revision(1))>=48 & double(revision(1))<=57,
            revision = str2num(revision(1));
        end
    else
        revision = 0;
    end
else
    % version number in form Major.Minor.Revision
    ver_num = sscanf(verstr, '%d.%d.%d');
    revision = ver_num(3);
end

major = ver_num(1);
minor = ver_num(2);
