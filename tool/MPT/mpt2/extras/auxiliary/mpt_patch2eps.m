function  mpt_patch2eps(name)
% Exporting 3-dimensional transparent patches to eps for use in LaTeX
% for the current figure
%
% USAGE:  fjc_patch2eps(name)
%      
%   name  - name-string of the eps-file
%
%   This produces one file 'name.eps' with the real figure
%   and one file 'name_axis.eps' with the corresponding axis.
%
%   Warning: that the grid will be lost.
%
%   in LaTeX use something like:
%     \usepackage{graphicx}
%      ...
%     \begin{figure}
%     \newlength{\picturehight}
%     \setlength{\picturehight}{8cm}
%     \includegraphics[height=\picturehight]{name.eps}
%     \vskip -1.003\picturehight  
%     \includegraphics[height=\picturehight]{name_axes.eps}
%     \end{figure}
%

% (C) 2004 Frank J. Christophersen
%          Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
%
% fjc, 2004-09-02

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



%make figure active
fh = figure(gcf);
hh = get(get(fh,'Children'),'Children');


% export the data with a white background
set(gcf,'PaperPositionMode','auto')
axis off
pstr = [name,'.eps'];
print('-depsc2','-loose',pstr);


% export the axes with transparent background
axis on
set(gca,'xlimmode','manual','xtickmode','manual','xticklabelmode','manual')
set(gca,'ylimmode','manual','ytickmode','manual','yticklabelmode','manual')
set(gca,'zlimmode','manual','ztickmode','manual','zticklabelmode','manual')
set(gca,'color','none')
set(gcf,'color','none')
set(gcf,'inverthardcopy','off') 
set(hh,'visible','off') % hide the data

astr = [name,'_axis','.eps'];
print('-depsc2','-loose',astr);


% in order to make the figure appear again
set(gca,'color','white')
set(gcf,'color',[ 0.800 0.800 0.800 ])
set(gcf,'inverthardcopy','on')
set(hh,'visible','on') % show the data
figure(gcf)


%% END
