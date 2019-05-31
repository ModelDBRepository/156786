%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the sparseness of the odor representation
%
% Licurgo de Almeida
% 05/03/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = OdorSparseness(V)

if iscell(V)
    Activity = zeros(length(V),1);
    for ii = 1:length(Activity)
        Activity(ii) = sum(V{ii}.S);
    end
    S = OdorSparseness(Activity);
else
    N = length(V);
    S = (1 - ((sum(V / N)^2) / sum((V.^2) / N))) / (1 - 1 / N);
end