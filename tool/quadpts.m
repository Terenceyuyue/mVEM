function [lambda,weight] = quadpts(order)
%% QUADPTS quadrature points in 2-D.
%
% [lambda,weight] = quadpts(order) return quadrature points with given
% order (up to 9) in the barycentric coordinates.
%
% The output lambda is a matrix of size nQ by 3, where nQ is the number of
% quadrature points. lambda(i,:) is three barycentric coordinate of the
% i-th quadrature point and lambda(:,j) is the j-th barycentric coordinate
% of all quadrature points. The x-y coordinate of the p-th quadrature point
% can be computed as 
%
%     pxy = lambda(p,1)*node(elem(:,1),:) ...
%         + lambda(p,2)*node(elem(:,2),:) ... 
%         + lambda(p,3)*node(elem(:,3),:);
%
% The weight of p-th quadrature point is given by weight(p). See
% verifyquadpts for the usage of qudrature rules to compute integrals over
% triangles.
% 
% References: 
%
% * David Dunavant. High degree efficient symmetrical Gaussian
%    quadrature rules for the triangle. International journal for numerical
%    methods in engineering. 21(6):1129--1148, 1985. 
% * John Burkardt. DUNAVANT Quadrature Rules for the Triangle.
%    http://people.sc.fsu.edu/~burkardt/m_src/dunavant/dunavant.html
% 
% See also quadpts1, quadpts3, verifyquadpts
%
% Order 6 - 9 is added by Huayi Wei, modify by Jie Zhou
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if order>9
    order = 9;
end
switch order
    case 1     % Order 1, nQuad 1
        lambda = [1/3, 1/3, 1/3];
        weight = 1;
    case 2     % Order 2, nQuad 3
        lambda = [2/3, 1/6, 1/6; ...
                  1/6, 2/3, 1/6; ...
                  1/6, 1/6, 2/3];
        weight = [1/3, 1/3, 1/3];
    case 3     % Order 3, nQuad 4
        lambda = [1/3, 1/3, 1/3; ...
                  0.6, 0.2, 0.2; ...
                  0.2, 0.6, 0.2; ...
                  0.2, 0.2, 0.6];
        weight = [-27/48, 25/48, 25/48, 25/48];
    case 4     % Order 4, nQuad 6
        lambda = [0.108103018168070, 0.445948490915965, 0.445948490915965; ...
                  0.445948490915965, 0.108103018168070, 0.445948490915965; ...
                  0.445948490915965, 0.445948490915965, 0.108103018168070; ...
                  0.816847572980459, 0.091576213509771, 0.091576213509771; ...
                  0.091576213509771, 0.816847572980459, 0.091576213509771; ...
                  0.091576213509771, 0.091576213509771, 0.816847572980459];
        weight = [0.223381589678011, 0.223381589678011, 0.223381589678011, ...
                  0.109951743655322, 0.109951743655322, 0.109951743655322];
    case 5     % Order 5, nQuad 7
        alpha1 = 0.059715871789770;      beta1 = 0.470142064105115;
        alpha2 = 0.797426985353087;      beta2 = 0.101286507323456;
        lambda = [   1/3,    1/3,    1/3; ...
                  alpha1,  beta1,  beta1; ...
                   beta1, alpha1,  beta1; ...
                   beta1,  beta1, alpha1; ...
                  alpha2,  beta2,  beta2; ...
                   beta2, alpha2,  beta2; ...
                   beta2,  beta2, alpha2];
        weight = [0.225, 0.132394152788506, 0.132394152788506, 0.132394152788506, ...
            0.125939180544827, 0.125939180544827, 0.125939180544827];
    case 6        
        A =[0.249286745170910  0.249286745170910  0.116786275726379
            0.249286745170910  0.501426509658179  0.116786275726379
            0.501426509658179  0.249286745170910  0.116786275726379
            0.063089014491502  0.063089014491502  0.050844906370207
            0.063089014491502  0.873821971016996  0.050844906370207
            0.873821971016996  0.063089014491502  0.050844906370207
            0.310352451033784  0.636502499121399  0.082851075618374
            0.636502499121399  0.053145049844817  0.082851075618374
            0.053145049844817  0.310352451033784  0.082851075618374
            0.636502499121399  0.310352451033784  0.082851075618374
            0.310352451033784  0.053145049844817  0.082851075618374
            0.053145049844817  0.636502499121399  0.082851075618374];
        lambda = [A(:,[1,2]), 1 - sum(A(:,[1,2]),2)];
        weight = A(:,3)';
    case 7
        A =[0.333333333333333  0.333333333333333 -0.149570044467682
            0.260345966079040  0.260345966079040  0.175615257433208
            0.260345966079040  0.479308067841920  0.175615257433208
            0.479308067841920  0.260345966079040  0.175615257433208
            0.065130102902216  0.065130102902216  0.053347235608838
            0.065130102902216  0.869739794195568  0.053347235608838
            0.869739794195568  0.065130102902216  0.053347235608838
            0.312865496004874  0.638444188569810  0.077113760890257
            0.638444188569810  0.048690315425316  0.077113760890257
            0.048690315425316  0.312865496004874  0.077113760890257
            0.638444188569810  0.312865496004874  0.077113760890257
            0.312865496004874  0.048690315425316  0.077113760890257
            0.048690315425316  0.638444188569810  0.077113760890257];
        lambda = [A(:,[1,2]), 1 - sum(A(:,[1,2]),2)];
        weight = A(:,3)';
    case 8
        A =[0.333333333333333  0.333333333333333  0.144315607677787
            0.081414823414554  0.459292588292723  0.095091634267285
            0.459292588292723  0.081414823414554  0.095091634267285
            0.459292588292723  0.459292588292723  0.095091634267285
            0.658861384496480  0.170569307751760  0.103217370534718
            0.170569307751760  0.658861384496480  0.103217370534718
            0.170569307751760  0.170569307751760  0.103217370534718
            0.898905543365938  0.050547228317031  0.032458497623198
            0.050547228317031  0.898905543365938  0.032458497623198
            0.050547228317031  0.050547228317031  0.032458497623198
            0.008394777409958  0.263112829634638  0.027230314174435
            0.008394777409958  0.728492392955404  0.027230314174435
            0.263112829634638  0.008394777409958  0.027230314174435
            0.728492392955404  0.008394777409958  0.027230314174435
            0.263112829634638  0.728492392955404  0.027230314174435
            0.728492392955404  0.263112829634638  0.027230314174435];
        lambda = [A(:,[1,2]), 1 - sum(A(:,[1,2]),2)];
        weight = A(:,3)';
        
        case 9
        A =[0.333333333333333  0.333333333333333  0.097135796282799
            0.020634961602525  0.489682519198738  0.031334700227139
            0.489682519198738  0.020634961602525  0.031334700227139
            0.489682519198738  0.489682519198738  0.031334700227139
            0.125820817014127  0.437089591492937  0.07782754100474
            0.437089591492937  0.125820817014127  0.07782754100474
            0.437089591492937  0.437089591492937  0.07782754100474
            0.623592928761935  0.188203535619033  0.079647738927210
            0.188203535619033  0.623592928761935  0.079647738927210
            0.188203535619033  0.188203535619033  0.079647738927210
            0.910540973211095  0.044729513394453  0.025577675658698
            0.044729513394453  0.910540973211095  0.025577675658698
            0.044729513394453  0.044729513394453  0.025577675658698
            0.036838412054736  0.221962989160766  0.043283539377289
            0.036838412054736  0.741198598784498  0.043283539377289
            0.221962989160766  0.036838412054736  0.043283539377289
            0.741198598784498  0.036838412054736  0.043283539377289
            0.221962989160766  0.741198598784498  0.043283539377289
            0.741198598784498  0.221962989160766  0.043283539377289];
        lambda = [A(:,[1,2]), 1 - sum(A(:,[1,2]),2)];
        weight = A(:,3)';
end
%% Verification
% The order of the quadrature rule is verified by the function
% verifyquadpts. See <matlab:ifem('verifyquadpts') verifyquadpts>.
