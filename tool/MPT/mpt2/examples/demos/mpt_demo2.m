function slide=mpt_demo2
%MPT_DEMO2 Tour through visualization capabilities of the toolbox
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Tour through visualization capabilities of the toolbox
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
%(C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%         kvasnica@control.ee.ethz.ch
%(C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%         grieder@control.ee.ethz.ch

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

if nargout<1,
  playshow mpt_demo2
else
    blue='b';
    red='r';
    slide(1).code={
        'clc; cla; axis([-4 4 -2 2]); Options.newfigure=0;',
        'axis([-4 4 -2 2]); grid off; axis off; text(-2.5,0.5,''Visualization tools'',''FontSize'',16);'};
    slide(1).text={
        'The Topic of this tutorial is to show some basic and advanced features of the MPT toolbox in terms of visualization.',
    };

    slide(2).code={
        'Options.newfigure=0; axis([-2 2 -2 2]);',
        'P=polytope([eye(2); -eye(2)], [1;1;1;1]);',
        'plot(P,Options);'
        'axis([-2 2 -2 2]);' };
    slide(2).text={
        'Let us first create a single polytope:',
        '',
        '>> P = polytope([eye(2); -eye(2)], [1;1;1;1])',
        '',
        'The basic command for plotting is:',
        '>> plot(P)',
        ''
        'This function plots polytope P in default color.'
    };

    slide(3).code={
        'cla; Options.newfigure=0;',
        'plot(P,''b'',Options);',
        'axis([-2 2 -2 2]);'};
    slide(3).text={
        'To plot a polytope in a different color, use the following call:',
        '',
        '>> plot(P, ''b'')',
        ''
        'This will plot the polytope in blue color.'
    };

    slide(4).code={
        'cla; Options.newfigure=0;',
        'Options.wire=1;',
        'plot(P,''b'',Options);',
        'axis([-2 2 -2 2]);'};
    slide(4).text={
        'To draw a polytope in wireframe (i.e. no color fill), the external Options structure has to be passed:',
        '',
        '>> Options.wire = 1',
        '>> plot(P, Options)',
        ''
        'You can see the wireframe on your screen now.'
    };

    slide(5).code={
        'cla; Options.newfigure=0;',
        'Options.wire=0;',
        'P1 = polytope([-4 2; -2 4; 0 2; 0 -2; -2 -4; -4 -2]);',
        'P2 = polytope([eye(2); -eye(2)], [1.5; 3; 0; 0.5]);',
        'P3 = polytope(rand(5,2)*4);',
        'plot(P,P1,P2,P3,Options);',
        'axis([-5 5 -5 5]);'};
    slide(5).text={
        'The plot function also accepts multiple polytope arguments as shown by this example:',
        '',
        '>> P1 = polytope([-4 2; -2 4; 0 2; 0 -2; -2 -4; -4 -2]);',
        '>> P2 = polytope([eye(2); -eye(2)], [1.5; 3; 0; 0.5]);',
        '>> P3 = polytope(rand(5,2)*4);',
        '>> plot(P, P1, P2, P3)',
        'or',
        '>> S = [P P1];',
        '>> plot(S, P2, P3);',
        '',
        'Color used for plotting is cycled through the default colormap.'
    };


    slide(6).code={
        'cla; Options.newfigure=0;',
        'Options.wire=1;',
        'plot(P,P1,P2,P3,Options);',
        'axis([-5 5 -5 5]);'};
    slide(6).text={
        'Or the same in the wireframe mode',
        '',
        '>> P1 = polytope([-4 2; -2 4; 0 2; 0 -2; -2 -4; -4 -2]);',
        '>> P2 = polytope([eye(2); -eye(2)], [1.5; 3; 0; 0.5]);',
        '>> P3 = polytope(rand(5,2)*4);',
        ''
        '>> Options.wire = 1;'
        ''
        '>> plot(P, P1, P2, P3, Options)'
    };

    slide(7).code={
        'cla; Options.newfigure=0;',
        'Options.wire=0;',
        'plot(P,''b'',P1,''r'', [P2 P3], ''g'',Options);',
        'axis([-5 5 -5 5]);'};
    slide(7).text={
        'When plotting multiple polytopes (or polytope arrays), you can specify each color individually:',
        '',
        '>> plot(P, ''b'', P1, ''r'', [P2 P3], ''g'')',
        ''
        'This plots P in blue, P1 in red and the polytope array [P2 P3] in green color.'
    };

    slide(8).code={
        'cla; Options.newfigure=0;',
        'Options.wire=0;',
        'P3D = polytope([eye(3); -eye(3)], [1;1;1;1;1;1]);',
        'plot(P3D,Options);',
        'axis([-2 2 -2 2 -2 2]);'
    };
    slide(8).text={
        'Plotting of polytopes in higher dimensions:',
        '',
        'The MPT toolbox provides a possibility to visualize a polytope in 3D. To show it, we create a unit cube:'
        '',
        '>> P3D = polytope([eye(3); -eye(3)], [1;1;1;1;1;1]);',
        '>> plot(P3D);',
        ''
        'Facets of the polytope are transparent by default.'
    };

    slide(9).code={
        'cla; Options.newfigure=0;',
        'Options.wire=0;',
        'Q1 = polytope(rand(8,3));',
        'Q2 = polytope(rand(10,3));',
        'Q3 = polytope(rand(9,3));',
        'Q4 = polytope(rand(11,3));',
        'Q5 = polytope(rand(13,3));',
        'S1 = [Q1 Q2 Q5]; S2 = [Q4 Q3];',
        'plot(S1,S2,Options);',
    };
    slide(9).text={
        'Again, the plot command can handle multiple input arguments also in 3D, as shown by the following example:',
        '',
        '>> Q1 = polytope(rand(8,3));',
        '>> Q2 = polytope(rand(10,3));',
        '>> Q3 = polytope(rand(9,3));',
        '>> Q4 = polytope(rand(11,3));',
        '>> Q5 = polytope(rand(13,3));',
        '',
        '>> S1 = [Q1 Q2 Q5]; S2 = [Q4 Q3];',
        '>> plot(S1, S2);',
        ''
        'Individual coloring can be specified by the user.'
    };

    slide(10).code={
        'cla; Options.newfigure=0;',
        'Options.wire=0;',
        'Options.xsection=2; Options.xvalues=0.5;',
        'P3D = polytope([eye(3); -eye(3)], [1;1;1;1;1;1]);',
        'plot(P3D,Options); axis([-2 2 -2 2]);'
    };
    slide(10).text={
        'It is sometimes useful (and also necessary for dimensions > 3) to plot a section through a polytope.',
        '',
        'This can be done by specifying the following option(s):',
        '',
        '>> Options.xsection = 2',
        '>> Options.xvalues = 0.5',
        'Says that we want to cut a given polytope along the 2nd dimension at value 0.5.',
        '',
        '>> plot(P3D, Options);',
        '',
        'If Options.xvalues is omitted, zero is assuming as default value.'
    };

    slide(11).code={
        'cla; clear Options; Options.newfigure=0;',
        'Options.wire=0;',
        'load pa4d',
        'plot(pa4d,Options);'
    };
    slide(11).text={
        'To demonstrate helpfulness of making cuts, let us consider a sample polytope array consisting of polytopes in 4D',
        '',
        '>> plot(pa4d);',
        '',
        ['Since the partition is in 4D, the cut through x4=0 is made automatically since no auxiliary Options parameter is given.'...
                'Sections of polytopes of dimension > 3 are plotted in 3D by default']
    };

    slide(12).code={
        'cla; Options.newfigure=0;',
        'Options.wire=0; Options.xsection=[3 4]; Options.xvalues=[-5 10];',
        'load pa4d',
        'plot(pa4d,Options);'
    };
    slide(12).text={
        'But the user can also specify 2-dimensional cuts by providing some extra arguments:',
        '',
        '>> Options.xsection = [3 4];',
        '>> Options.xvalues = [-5 10];',
        '',
        '>> plot(pa4d, Options);',
        '',
        'The cut through x3 = -5 and x4 = 10 is now depicted on your screen. If Options.xvalues is not given, zero values are assumed by default.'
    };

    %end
    slide(13).code={
        'cla; axis([-2 2 -2 2]); grid off; axis off; text(-0.5,0.5,''The end'',''FontSize'',16); title(''''); ',
    };
    slide(13).text={
        ''
    };

end
