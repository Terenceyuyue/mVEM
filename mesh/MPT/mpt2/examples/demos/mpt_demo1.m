function slide=mpt_demo1
%MPT_DEMO1 Explains basic manipulation with polytopes
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Explains basic manipulation with polytopes
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
  playshow mpt_demo1
else
    %intro

    slide(1).code={
        'clc; cla; axis([-4 4 -2 2]); Options.newfigure=0;',
        'axis([-4 4 -2 2]); grid off; axis off; text(-4,0.5,''Basic manipulation with polytopes'',''FontSize'',16);'};
    slide(1).text={
        'This tour presents some of the basic functions',
        'for polytope manipulation',
        '',
        '  - Creation of polytope objects',
        '  - Basic manipulation of polytopes',
        '  - Computational geometry tools',
        '  - Visualization',
        '  - Accessing internal information stored in the polytope object',
    };

    %creating polytopes

    slide(2).code={
        '' };
    slide(2).text={
        'A Polytope is defined as a bounded intersection of a finite number of hyperplanes.',
        '',
        'In general, a polytope object can be initialized in two ways:'
        '- providing the H-representation (i.e. intersecting hyperplanes)',
        '- by computing convex hull of extreme points (vertices)',
        '',
        ['The polytope is internally always converted to an H-representation.'...
        'However, if the vertices are provided by the user, or computed on purpose,'...
        'they are also stored in the internal structure and are available on request.'],
    };

    %H representation

    slide(3).code={
        'disp(''Polytope P1:'');'
        'P1=polytope([eye(2); -eye(2)],[1;1;1;1])',
        'plot(P1,Options); axis([-2 2 -2 2]);' };
    slide(3).text={
        'Let us now create a polytope using the H-representation:',
        '',
        '>> H = [eye(2); -eye(2)];',
        '>> K = [1; 1; 1; 1];';
        '>> P1 = polytope(H, K)',
        '',
        'The H and K matrices essential state that:',
        'x1<=1, x2<=1, x1>=-1, x2>=-1',
        '',
        'The created polytope is plotted on your screen.'
    };

    %V representation

    slide(4).code={
        'fprintf(''\nPolytope P2:\n'');'
        'P2=polytope([0.2 -0.2; 0.2 0.4; -0.3 0.4; -0.3 -0.5])',
        'plot(P2,''b'',Options); axis([-2 2 -2 2]);' };
    slide(4).text={
        'Now we define a polytope having the following vertices:',
        'V1 = [0.2 -0.2], V2=[0.2 0.4], V3=[-0.3 0.4], V4=[-0.3 -0.5]',
        '',
        'This can be done by typing:'
        '>> V = [0.2 -0.2; 0.2 0.4; -0.3 0.4; -0.3 -0.5]',
        '>> P2 = polytope(V)',
        '',
        'The resulting polytope is stored in the workspace and plotted on your screen.'
    };

    %extreme and hull

    slide(5).code={
        'cla; plot(P1,Options); axis([-2 2 -2 2])',
        'clc; fprintf(''\nExtreme points of polytope P1:\n'');'
        'vert=extreme(P1)',
        'fprintf(''\nConvex hull of vertices as H-representation polytope:\n'');',
        'p=hull(vert)'};
    slide(5).text={
        'It is always possible to go from V to H representation and vice versa.',
        '',
        '>> vert = extreme(P1)',
        'Computes extreme points (vertices) of polytope P1',
        ''
        '>> p = hull(vert)',
        'Calculates convex hull of vertices given by matrix "vert"'
        '',
        'see help extreme for more details.'
    };

    %creating arrays

    slide(6).code={
        'clc; fprintf(''\nPolytope array S:\n'');',
        'P3 = polytope([eye(2); -eye(2)],[1.5; 1.5; 0; 0])',
        'S=[P1 P2 P3], plot(S,Options); axis([-2 2 -2 2]);' };
    slide(6).text={
        'Polytopes can be concatenated into arrays:',
        '',
        '>> P3 = polytope([eye(2); -eye(2)], [1.5; 1.5; 0; 0])',
        '>> S = [P1 P2 P3];',
        '>> S2 = [P3 S P2];',
        '',
        'Same functions can be used for polytope arrays as well as single polytopes, e.g.:',
        '',
        '>> plot(S)',
        '',
        'plots polytopes that are stored in the array S.'
    };

    %indexing of arrays

    slide(7).code={
        'clc; cla; fprintf(''\nIndexing a polytope array:\n'');',
        'R1=S(1), fprintf(''\n'');',
        'R2=S(1:3), fprintf(''\n'');',
        'R3=S([3 2]), fprintf(''\n'');',
        'R4=S(end), fprintf(''\n'');'};

    slide(7).text={
        'Subindexing can be used to access individual polytopes stored in an array S=[P1 P2 P3]',
        '',
        '>> R1=S(1)',
        'Returns the first element of polytope array S.',
        ''
        '>> R2=S(1:3)',
        'Returns polytopes stored at positions 1 to 3 in S.',
        ''
        '>> R3=S([3 2])',
        'Returns polytopes stored at positions 3 and 2. It is the same as typing R3=[S(3) S(2)].',
        ''
        '>> R4=S(end)',
        'Returns the final element of S. It is identical to R4=S(length(S)).',
        '',
        'See the workspace for results.'
    };


    %relational operators
    slide(8).code={
        'cla; axis off; plot([P1 P2],Options); axis([-2 2 -2 2]);',
        'clc; fprintf(''\nP2==P1''); P2==P1',
        'fprintf(''\nP2<P1''); P2<P1',
        'fprintf(''\nP2<=P1''); P2<=P1',
    };
    slide(8).text={
        'Basic manipulation with polytopes (and polytopes arrays) includes.',
        '',
        'One can use a set of overloaded operators to perform operations on polytopes:',
        '',
        'P2 == P1   Tests if two polytopes are equal (returns 0 (false) in our case)',
        'P2 < P1  Tests if P1 is a strict subset of P2 (returns 1 (true) in out case)',
        'P2 <= P1  Tests if P1 is a subset of P2 (also true for the two depicted polytopes)',
        '',
        'Polytope P1 is plotted in cyan and P2 in red color.',
        '',
        'You can find more about relational operators by typing:',
        '>> help mpt/polytope',
    };


    %computational geometry
    slide(9).code={
        'cla;clc; axis off',
    };
    slide(9).text={
        ['MPT toolbox provides a set of routines from the field of computational geometry'...
                ' that can be performed on polytopes:'],
        '',
        '- Intersection',
        '- Convex union',
        '- Convex hull',
        '- Convex envelope',
        '- Minkowski addition'
        '- Pontryagin difference',
        '- Set difference',
    };

    %intersection
    slide(10).code={
        'clc; cla;Options.wire=1; plot(P1,P3,''k'',Options);',
        'hold on; Options.wire=0; plot(P1&P3,Options); axis([-2 2 -2 2])',
        'hold off'
    };
    slide(10).text={
        'Intersection of polytopes and/or polytope arrays'
        '',
        '>> I = P1 & P3',
        'Computes intersection of polytopes P1 and P3',
        ''
        'P1 and P3 are depicted in black wireframe, the intersection in red color.'
    };

    %convex union
    slide(11).code={
        'cla;clc; Options.wire=0;',
        'q1=polytope([-1 0; 0 1; 0 -1]);',
        'q2=polytope([0 1; 1 0; 0 -1]);',
        'plot(q1|q2,Options); axis([-2 2 -2 2])',
        'hold on;',
        'Options.wire=1; Options.wirecolor=''k'';',
        'plot(q1,q2,Options); axis([-2 2 -2 2]);',
        'hold off'
    };
    slide(11).text={
        'If union of two polytopes is convex, it can be computed using the overloaded | operator.'
        '',
        'First we define some polytopes:',
        '>> Q1=polytope([-1 0; 0 1; 0 -1]);',
        '>> Q2=polytope([0 1; 1 0; 0 -1]);',
        ''
        'And then we can compute the union:'
        '>> U = Q1 | Q2',
        '',
        'Polytopes Q1 and Q2 are plotted in black wireframe while U in red.',
        '',
        'Note that if the union were not convex, polytope array [Q1 Q2] would be returned.'
    };

    %convex hull
    slide(12).code={
        'cla;clc;Options.wire=0;',
        'plot(hull([P1 P2 P3]),Options); axis([-2 2 -2 2])',
        'hold on;',
        'Options.wire=1;',
        'plot(P1,P2,P3,Options); axis([-2 2 -2 2]);',
        'Options.wire=0; hold off'
    };
    slide(12).text={
        'The convex hull of n polytopes can be computed as follows:'
        '',
        '>> S = [P1 P2 P3];',
        '>> H = hull(S)',
        '>> plot(H)',
        '',
        ['On the figure, polytopes P1, P2 and P3 are plotted in wireframe, the convex hull is depicted in red. '...
                'Note that the hull function automatically eliminates non-extreme points.']
    };


    %minkowski addition
    slide(13).code={
        'cla;clc;Options.wire=0;',
        'plot(P1+P2,Options); axis([-2 2 -2 2])',
        'hold on;',
        'Options.wire=1;',
        'plot(P1,P2,Options); axis([-2 2 -2 2]);',
        'Options.wire=0; hold off'
    };
    slide(13).text={
        ['Let us assume the two polytopes as plotted on your screen. Then, the Minkowski sum of these'...
                'two polytopes can be calculated by typing:'],
        '',
        '>> Madd = P1 + P2',
        '>> plot(Madd)',
        '',
        'Note: Polytopes P1 and P2 are plotted in wireframe, the Minkowski sum is depicted in red.'
    };

    %pontryagin difference
    slide(14).code={
        'cla;clc;Options.wire=0;',
        'P1 = polytope([eye(2); -eye(2)], [1;1;1;1]);',
        'P2 = polytope([eye(2); -eye(2)], [1;1;1;1]*0.2);',
        'plot(P1-P2,Options); axis([-2 2 -2 2])',
        'hold on;',
        'Options.wire=1;',
        'plot(P1,P2,Options); axis([-2 2 -2 2]);',
        'Options.wire=0; hold off'
    };
    slide(14).text={
        'Similarly, one can compute Pontryagin difference by using the overloaded minus (-) operator:',
        '',
        '>> P1 = polytope([eye(2); -eye(2)], [1;1;1;1]);',
        '>> P2 = polytope([eye(2); -eye(2)], [1;1;1;1]*0.2);',
        '>> Pdiff = P1 - P2',
        '>> plot(Pdiff)',
        '',
        'Note: Polytopes P1 and P2 are plotted in wireframe, the Pontryagin difference is depicted in red.'
    };

    %pontryagin difference with sets
    slide(15).code={
        'cla;clc;Options.wire=0;',
        'P1 = polytope([-4 2; -2 4; 0 2; 0 -2; -2 -4; -4 -2]);',
        'P2 = -P1;',
        'P3 = polytope([eye(2); -eye(2)], [1;1;1;1]*0.5);',
        'plot([P1 P2]-P3,Options); axis([-5 5 -5 5])',
        'hold on;',
        'Options.wire=1; Options.wirewidth=2;',
        'Options.wirecolor=''r''; plot(P1,Options); axis([-5 5 -5 5]);',
        'Options.wirecolor=''b''; plot(P2,Options); axis([-5 5 -5 5]);',
        'Options.wirecolor=''k''; plot(P3,Options); axis([-5 5 -5 5]);',
        'Options.wire=0; hold off'
    };
    slide(15).text={
        'Pontryagin difference also accepts polytope arrays, as shown in this example:',
        '',
        '>> P1 = polytope([-4 2; -2 4; 0 2; 0 -2; -2 -4; -4 -2]);',
        '>> P2 = -P1',
        'Mirrors polytope P1 around the origin',
        '',
        '>> P3 = polytope([eye(2); -eye(2)], [1;1;1;1]*0.5);',
        '>> Pdiff = [P1 P2] - P3',
        '>> plot(Pdiff)',
        '',
        'In this case, the solution consists of a finite number of polytopes.'
        'Note: Polytope P1 is in red wireframe, P2 in blue wireframe, P3 in black and the Pontryagin difference in full color.'
    };

    %set difference

    slide(16).code={
        'clc;cla; Options.wire=0;',
        'P1=polytope([-3 3; 0 3; -3 -3; 0 -3]);',
        'P2=polytope([0 0; 5 0; 5 -3; 0 -3]);',
        'P3=polytope([eye(2); -eye(2)], [1;1;1;1]);',
        'P4=polytope([-1 -2; 2 -3; -1 -3]);',
        'plot(P1,P2,Options); hold on',
        'Options.wirecolor=''k''; Options.wire=1;',
        'plot(P3,P4,Options); hold off;'
    };
    slide(16).text={
        'To show the computation of set differences, let us define following polytopes:',
        '',
        '>> P1=polytope([-3 3; 0 3; -3 -3; 0 -3]);',
        '>> P2=polytope([0 0; 5 0; 5 -3; 0 -3]);',
        '>> P3=polytope([eye(2); -eye(2)], [1;1;1;1]);',
        '>> P4=polytope([-1 -2; 2 -3; -1 -3]);',
        'P1 and P2 are depicted in full color and P3 and P4 are drawn in wireframe.',
        '',
        'The set difference between set S1 = [P1 P2] and set S2 = [P3 P4] can be then computed as follows:',
        '',
        '>> S1 = [P1 P2]',
        '>> S2 = [P3 P4]',
        '>> Sdiff = S1 \ S2',
        '',
        'The set difference is displayed on the next slide...'
    };

    slide(17).code={
        'clc;cla; Options.wire=0;',
        'P1=polytope([-3 3; 0 3; -3 -3; 0 -3]);',
        'P2=polytope([0 0; 5 0; 5 -3; 0 -3]);',
        'P3=polytope([eye(2); -eye(2)], [1;1;1;1]);',
        'P4=polytope([-1 -2; 2 -3; -1 -3]);',
        'plot([P1 P2]\[P3 P4],Options);'};
    slide(17).text={
        '...continuation of the previous slide'
        '',
        '>> S1 = [P1 P2]',
        '>> S2 = [P3 P4]',
        '>> Sdiff = S1 \ S2',
    };


    %

    slide(18).code={
        'clc; cla; axis off'};
    slide(18).text={
        'The rest of this tutorial will be devoted to the introduction to the internal structure of a polytope object.'
        '',
        ['We will show which functions can be used to access internal (private) data of a polytope. This routines'...
                ' are mainly useful in connection with control-design and analysis routines of MPT toolbox.']
    };

    %double

    slide(19).code={
        'clc; fprintf(''\nH-representation of polytope P:\n'');',
        'P = polytope([eye(2); -eye(2)], [1;1;1;1]);',
        '[H,K]=double(P1)' };
    slide(19).text={
        'The H-representation of a polytope can be accessed using the overloaded double function:',
        '',
        '>> P = polytope([eye(2); -eye(2)], [1;1;1;1])',
        '>> [H,K] = double(P)',
        '',
        'This command extracts the H-representation of polytope P and stores it into matrices H and K.'
    };

    %isnormal, isminrep

    slide(20).code={
        'clc; fprintf(''\nStatus of the H-representation for polytope P1:\n'');',
        'isnorm=isnormal(P), isred=isminrep(P)' };
    slide(20).text={
        ['When a polytope object is created, its H-representation is automatically reduced and normalized.'...
            'This behavior, however, can be suppressed by the user. To obtain information about normalization and reduction,'...
            ' you can use:'],
        '',
        '>> isnorm = isnormal(P)',
        '>> isred  = isminrep(P)',
        '',
        'The first command returns 1 (true) if the polytope is in normalized H-representation.',
        'The second command will return 1 (true) if the polytope is in minimal H-representation (i.e. if redundant constraints were removed)'
    };

    %chebyball

    slide(21).code={
        'clc; fprintf(''\nCenter and radius of Chebyshevs ball:\n'');', 
        '[xCheb, Rcheb] = chebyball(P1)' };
    slide(21).text={
        ['When you create a polytope, center and radius of the largest ball that can be subscribed'...
                ' into that polytope (Chebyshev''s ball) is computed automatically. To access this information'...
                ' you can use the following command:'],
        '',
        '>> [xCheb, Rcheb] = chebyball(P1)',
        '',
        'You see output of the command in your workspace.'
    };

    %isinside
    slide(22).code={
        'clc; fprintf(''\nChecking if a given point lies in a polytope:\n'');',
        'isin = isinside(P, [0.1; 0.2])' };
    slide(22).text={
        ['To find out if a given point lies in the interior of some polytope, one can use this function:'],
        '',
        '>> isin = isinside(P,[0.1; 0.2])',
        '',
        'Returns 1 (true) if a points (x1=0.1, x2=0.2) lies in interior of polytope P.'
    };

    %end
    slide(23).code={
        'cla; axis([-2 2 -2 2]); grid off; axis off; text(-0.5,0.5,''The end'',''FontSize'',16);',
    };
    slide(23).text={
        'Please consult',
        '',
        '>> help mpt/polytope',
        '',
        'for more details.'
    };

end
