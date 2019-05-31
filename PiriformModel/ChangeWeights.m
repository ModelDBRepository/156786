%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates a new weight matrix between cells
%
% Licurgo de Almeida
% 03/28/2011
%
% This function receives the following parameters:
% Cell: type of cell (Pyramidal, Feedforward, etc)
% param: list of parameters from the model
% plotweights: if true, the weight matrix is plotted
% Taus: is used for test of different parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W,Conn] = ChangeWeights(Cell,param,plotweights,Tlearning)
x = 1;
% Set parameteres
param.Tau11 = 50 * x; %ms
param.Tau01 = 250 * x; %ms
param.Tau10 = 250 * x; %ms

param.Taupost = 2; %ms
param.TauNMDAf = 7; %ms
param.TauNMDAr = 1; %ms
param.Tdelay = 1; %ms


Conn = zeros(length(Cell),length(Cell{1}.PyrCon));
for ii = 1:length(Cell)
    Conn(ii,:) = Cell{ii}.PyrCon;
end
Savefile = strcat(param.outputPath,'LearnData'); %pyr-pyr learned weights
load(Savefile,'W');


S = zeros(length(Cell),length(Cell{1}.S));

for ii = 1:length(Cell)
    S(ii,:) = Cell{ii}.S;
end


[ipost,bglu] = WeightCurves(S,param);

if nargin < 4
    Tlearning = size(ipost,2);
else
    Tlearning = round(Tlearning / param.dt);
end

W = CalculateWeights(W,Conn,ipost,bglu,param,Tlearning);
Wsource = param.Mitralsource(1:5);

Savefile = strcat(param.outputPath,'LearnData'); %pyr-pyr learned weights
save(Savefile,'W','Wsource');

Wplot = W ./ max(max(W));

if plotweights == true
    PlotWeights(Wplot,Wsource);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This functions creates the curves for the synaptic weight changes
%
% Licurgo de Almeida
% 03/29/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ipost,bglu] = WeightCurves(S,param)

t0 = zeros(size(S,1),1);
t0 = t0 - 10000000000;

ipost = zeros(size(S));
bglu = ipost;

for ii = 1:size(S,2)
    t = ii * param.dt;
    I = S(:,ii) == 1;
    t0(I) = t;
    ipost(:,ii) = ((t - t0) ./ param.Taupost) .* exp(1 - (t - t0) ./ param.Taupost);
    
    bglu(:,ii) = exp(-((t - t0) ./ param.TauNMDAf)) .* (1 - exp(-((t - t0) ./ param.TauNMDAr)));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the changes in the synaptic weights between two
% cells
%
% Licurgo de Almeida
% 03/30/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = CalculateWeights(W,Conn,ipost,bglu,param,Tlearning)

delay = round(param.Tdelay / param.dt);
b = circshift(bglu,[0 -delay]);
Conn = Conn';
i = ipost ./ param.Tau01;
bb = b ./ param.Tau10;

for ii = 1:size(Conn,1)
    auxW = W(ii,:)';
    for jj = 1:Tlearning - delay
        auxW = auxW + param.dt * (((ipost(ii,jj) .* b(:,jj) .* Conn(:,ii)) ./ param.Tau11) .* (1 - auxW)...
            + (i(ii,jj) + bb(:,jj) .* Conn(:,ii)) .* (0 - auxW));
    end
    W(ii,:) = auxW';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function plots the weight matrix between cells
%
% Licurgo de Almeida
% 03/30/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotWeights(W,Wsource)

sizeW = 10;

[I,J] = find(W);

auxW = ceil(W * sizeW);
figure;
hold on;
cla;
for ii = 1:length(I)
    plot(I(ii),J(ii),'sk','markersize',auxW(I(ii),J(ii)),'markerfacecolor','k');
end
xlabel('neuron','fontsize',12);
ylabel('neuron','fontsize',12);
if strcmp(Wsource,'AChOn')
    title('Mitral input = ACh ON','fontsize',12);
elseif strcmp(Wsource,'AChOf')
    title('Mitral input = ACh OFF','fontsize',12);
end
    
end