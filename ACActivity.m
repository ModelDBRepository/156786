%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function gives the average activity of the most active cells
%
% Licurgo de Almeida
% 09/26/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,N] = ACActivity(Cell,param)

F = ExtractActivity(Cell,param);
I = F < mean(F) + (2 * std(F)); %criterion adopted in the paper.
F(I) = [];
if ~isempty(F)
    A = mean(F);
else
    A = 0;
end

if nargout > 1
    N = length(F);
end