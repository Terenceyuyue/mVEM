%  h = fourier(H,ax,tol,qtol)
%
%  C implementation of the fourier elimination algorithm
%
% Inputs:
%  H           = [A b] - Ax < b defines a polyhedron
%  ax          = dimensions of A onto which we want to project
%  tol  [1e-6] = numerical tolerance before two numbers are considered equal
%  qtol [1e-2] = angle in degrees between two vectors before they are
%                considered equal
%
% Outputs:
%  h = [A b] - Ax < b defines the projection of the polyhedron
%
%  
%  Author: Colin Jones
%  Email : cnj22@cam.ac.uk
%  Date  : April 17, 2004

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
