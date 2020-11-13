% Multi-Parametric Toolbox Polytope library
% Version 2.6.3 (R2008a) 14-Jan-2008
%
%
% Constructor and data accessing methods
%   polytope     - Default constructor for polytope objects
%   double       - Function used to access internal properties of the given polytope
%   display      - Displays details about the given polytope
%                
%   isbounded    - Checks if a polytope is bounded
%   isconvex     - Checks if a polytope array forms a convex union
%   isfulldim    - Checks if a polytope is full dimensional
%   isinside     - Checks if a given point lies in the interior of a given polytope
%   isminrep     - Checks if a given polytope is in minimal representation
%   isnormal     - Checks if a given polytope is in normalized description
%                
%   dimension    - Returns dimension of the given polytope
%   nconstr      - Returns number of constraints that form an H-representation of a polytope
%   length       - Returns number of elements in a polytope array
%   size         - Returns size of the given polytope object
%                
%   get          - Get polytope properties
%   set          - Used to modify internal properties of a given polytope object
%
% Computational geometry functions
%   bounding_box - Compute a bounding box for a given polytope
%   domain       - Computes polytope that is mapped to an another polytope using affine map
%   envelope     - Computes envelope of n polytopes
%   extreme      - Calculates extreme points of a given polytope
%   facetvoronoi - Calculates a "voronoi diagram for facets of a polytope"
%   hull         - Convex hull of n polytopes
%   intersect    - Intersection of 2 polytopes
%   plus         - Minkowski sum
%   minus        - Pontryagin difference
%   mldivide     - Set difference
%   mtimes       - Polytope multiplication
%   projection   - Projection of a polytope or a polytope array
%   range        - Affine transformation of a polytope
%   regiondiff   - Region difference
%   regiondiffXU - Computes region difference in lifted XU space
%   slice        - Orthogonal cut through polytope(s)
%   triangulate  - Calculates triangulation of arbitrary polytopes
%   union        - convex union computation
%
% Overloaded relational operators
%   eq           - Checks if two polytopes are equal                        (==)
%   ge           - Checks if polytope P is a superset of polytope Q         (>=)
%   gt           - Checks if polytope P is a strict superset of polytope Q  (>)
%   le           - Checks if polytope P is a subset of polytope Q           (<=)
%   lt           - Checks if polytope P is a strict subset of polytope Q    (<)
%   ne           - Checks if two polytopes are not equal                    (~=)
%
% Basic polytope functions
%   chebyball    - Computes center and radius of the largest ball inscribed in a polytope
%   normalize    - Normalizes a given polytope
%   reduce       - Reduces the polytope by removing redundant inequalities
%   isredundant  - Check if a constraint is redundant
%   plot         - Plots polytope(s)
%
% Other overloaded operators and functions
%   and          - Intersection of n polytopes                              (&)
%   or           - Convex union of n polytopes                              (|)
%   plus         - Minkowski sum                                            (+)
%   minus        - Pontryagin difference                                    (-)
%   mldivide     - Set difference                                           (\)
%   mtimes       - Polytope multiplication                                  (*)
%   uminus       - Unary minus for a polytope (flips given polytope around the origin)
%   uplus        - Unary plus for a polytope (reduces the given polytope)
%   horzcat      - Concatenates polytopes into a polytope array             [,]
%   vertcat      - Concatenates polytopes into a polytope array             [;]
%   end          - Returns the last element in a given polytope array
%   fliplr       - Flips array of polytopes
%   flipud       - Flips array of polytopes
%   subsasgn     - Indexed assignments for polytope objects
%   subsref      - Indexed referencing for polytope objects
%
% Other functions
%   distribpoints - Distribute n points such that distances between them are maximized
%   facetcircle   - Returns largest circle inside facet 'ind' of polytope P
%   pelemfun      - Execute an arbitrary function on each element of a polyarray
%   reduceunion   - Removes some redundant elements from a polytope array
%   unique        - Keeps only unique elements of a polytope array
%   volume        - Calculates volume of a polytope
%
%
% see help mpt
%
% Authors: Michal Kvasnica, Pascal Grieder, Mato Baotic
% kvasnica@control.ee.ethz.ch, grieder@control.ee.ethz.ch, baotic@control.ee.ethz.ch
%
% For support, write to: mpt@control.ee.ethz.ch

