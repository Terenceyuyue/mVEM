function Ucombs = mpt_allcombs(U)
%MPT_ALLCOMBS All combinations of discrete-valued inputs
%
% Ucombs = mpt_allcombs(U)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns all possible combinations of discrete-valued inputs
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% U       - cell array of discrete values of the system inputs
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% Ucombs  - combinations of inputs (row-wise)
%

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

if ~iscell(U),
    Ucombs = U(:);
    return
end

Ucombs = [];
switch length(U)
    case 1
        Ucombs = U{1}(:);
        return;
    case 2
        for i1=1:length(U{1}),
            for i2=1:length(U{2}),
                Ucombs = [Ucombs; U{1}(i1) U{2}(i2)];
            end
        end
    case 3
        for i1=1:length(U{1}),
            for i2=1:length(U{2}),
                for i3=1:length(U{3}),
                    Ucombs = [Ucombs; U{1}(i1) U{2}(i2) U{3}(i3)];
                end
            end
        end
    case 4
        for i1=1:length(U{1}),
            for i2=1:length(U{2}),
                for i3=1:length(U{3}),
                    for i4=1:length(U{4}),
                        Ucombs = [Ucombs; U{1}(i1) U{2}(i2) U{3}(i3) U{4}(i4)];
                    end
                end
            end
        end
    case 5
        for i1=1:length(U{1}),
            for i2=1:length(U{2}),
                for i3=1:length(U{3}),
                    for i4=1:length(U{4}),
                        for i5=1:length(U{5}),
                            Ucombs = [Ucombs; U{1}(i1) U{2}(i2) U{3}(i3) U{4}(i4) U{5}(i5)];
                        end
                    end
                end
            end
        end
    case 6
        for i1=1:length(U{1}),
            for i2=1:length(U{2}),
                for i3=1:length(U{3}),
                    for i4=1:length(U{4}),
                        for i5=1:length(U{5}),
                            for i6=1:length(U{6}),
                                Ucombs = [Ucombs; U{1}(i1) U{2}(i2) U{3}(i3) U{4}(i4) U{5}(i5) U{6}(i6)];
                            end
                        end
                    end
                end
            end
        end
    case 7
        for i1=1:length(U{1}),
            for i2=1:length(U{2}),
                for i3=1:length(U{3}),
                    for i4=1:length(U{4}),
                        for i5=1:length(U{5}),
                            for i6=1:length(U{6}),
                                for i7=1:length(U{7}),
                                    Ucombs = [Ucombs; U{1}(i1) U{2}(i2) U{3}(i3) U{4}(i4) U{5}(i5) U{6}(i6) U{7}{i7}];
                                end
                            end
                        end
                    end
                end
            end
        end
    case 8
        for i1=1:length(U{1}),
            for i2=1:length(U{2}),
                for i3=1:length(U{3}),
                    for i4=1:length(U{4}),
                        for i5=1:length(U{5}),
                            for i6=1:length(U{6}),
                                for i7=1:length(U{7}),
                                    for i8=1:length(U{8}),
                                        Ucombs = [Ucombs; U{1}(i1) U{2}(i2) U{3}(i3) U{4}(i4) U{5}(i5) U{6}(i6) U{7}{i7} U{8}(i8)];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    case 9
        for i1=1:length(U{1}),
            for i2=1:length(U{2}),
                for i3=1:length(U{3}),
                    for i4=1:length(U{4}),
                        for i5=1:length(U{5}),
                            for i6=1:length(U{6}),
                                for i7=1:length(U{7}),
                                    for i8=1:length(U{8}),
                                        for i9=1:length(U{9}),
                                            Ucombs = [Ucombs; U{1}(i1) U{2}(i2) U{3}(i3) U{4}(i4) U{5}(i5) U{6}(i6) U{7}{i7} U{8}(i8) U{9}(i9)];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    case 10
        for i1=1:length(U{1}),
            for i2=1:length(U{2}),
                for i3=1:length(U{3}),
                    for i4=1:length(U{4}),
                        for i5=1:length(U{5}),
                            for i6=1:length(U{6}),
                                for i7=1:length(U{7}),
                                    for i8=1:length(U{8}),
                                        for i9=1:length(U{9}),
                                            for i10=1:length(U{10}),
                                                Ucombs = [Ucombs; U{1}(i1) U{2}(i2) U{3}(i3) U{4}(i4) U{5}(i5) U{6}(i6) U{7}(i7) U{8}(i8) U{9}(i9) U{10}(i10)];
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

    otherwise
        error('Sorry, maximum of 10 inputs are supported. If you got this warning, write an email to mpt@control.ee.ethz.ch and we will help you');
end