%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function extracts the average frequency of a group of cells
%
% Licurgo de Almeida
% 05/05/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = ExtractActivity(Cell,param)
if iscell(Cell)
    A = zeros(length(Cell),1);
    for ii = 1:length(Cell)
        A(ii) = sum(Cell{ii}.S) / (param.tsim / 1000);
    end
else
    A = sum(Cell,2) / (param.tsim / 1000);
end