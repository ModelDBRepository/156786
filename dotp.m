%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the normalized dot product between two vectors
%
% Licurgo de Almeida
% 08/30/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = dotp(v1,v2)

d = dot(v1,v2);
n1 = norm(v1);
n2 = norm(v2);
D = 1 - (d / (n1 * n2));