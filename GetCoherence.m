%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the coherence in the network activity. The
% function is based on the method proposed in Wang and Buzsáki, 1996
%
% Licurgo de Almeida
% 06/22/2011
% method = L for Linster; B for Buzsáki; N for Null test.
% flag_removelow =  true: Remove cells with low activity (default = false);
% flag_fixtime = true: time is fixed (default = true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function netcoh = GetCoherence(NS,method,flag_removelow,flag_fixtime)
if nargin < 4
    flag_fixtime = true;
end
if nargin < 3
    flag_removelow = false;
end

dt = 0.5;

if iscell(NS)
    aux_NS = zeros(length(NS),length(NS{1}.S));
    for ii = 1:length(NS)
        aux_NS(ii,:) = NS{ii}.S;
    end
    netcoh = GetCoherence(aux_NS,method,flag_removelow,flag_fixtime);
else
    if flag_removelow == true
        S = sum(NS,2);
        I = S <= mean(S);
        NS(I,:) = [];
    end
    if flag_fixtime == true
        twindow = 2; %ms
    else
        time = size(NS,2) * dt;
        twindow = 0.1 * (1 / (mean(sum(NS,2)) / time));

    end
    
    sizeblock = twindow / dt;
    
    switch lower(method)
        case 'l'
            netcoh = Linster(sizeblock,NS);
        case 'b'
            netcoh = Buzsaki(sizeblock,NS);
        case 'n'
            netcoh = NullTest(sizeblock,NS);
    end
end
end

function netcoh = Linster(sizeblock,NS)
netcoh = 0;

for ii = 1:size(NS,1) - 1
    for jj  = ii + 1:size(NS,1)
        sv1 = sum(NS(ii,:));
        sv2 = sum(NS(jj,:));
        if sv1 <= sv2
            mins = find(NS(ii,:));
            maxs = find(NS(jj,:));
        else
            maxs = find(NS(ii,:));
            mins = find(NS(jj,:));
        end
        if ~isempty(mins)
            sumcoh = 0;
            for kk = 1:length(mins)
                sumcoh = sumcoh + sum(abs(maxs - mins(kk)) <= sizeblock);
            end
            netcoh = netcoh + sumcoh;
        end
    end
end
netcoh = netcoh / sum(sum(NS));
end

function netcoh = Buzsaki(sizeblock,NS)
netcoh = 0;
countcoh = 0;
for ii = 1:size(NS,1) - 1
    for jj  = ii + 1:size(NS,1)
        sv1 = sum(NS(ii,:));
        sv2 = sum(NS(jj,:));
        if sv1 <= sv2
            mins = find(NS(ii,:));
            maxs = find(NS(jj,:));
        else
            maxs = find(NS(ii,:));
            mins = find(NS(jj,:));
        end
        if ~isempty(mins)
            sumcoh = 0;
            for kk = 1:length(mins)
                sumcoh = sumcoh + sum(abs(maxs - mins(kk)) <= sizeblock);
            end
            countcoh = countcoh + 1;
            netcoh = netcoh + sumcoh / sqrt(length(maxs) * length(mins));
        end
    end
end
netcoh = netcoh / countcoh;
end

function netcoh = NullTest(sizeblock,NS)
netreal = Linster(sizeblock,NS);
randNS = rand(size(NS));
A = sum(NS,2);
A = A / size(NS,2);
for ii = 1:length(A)
    randNS(ii,:) = randNS(ii,:) <= A(ii);
end
netnull = Linster(sizeblock,randNS);
netcoh = 1 - (netnull / netreal);
end
