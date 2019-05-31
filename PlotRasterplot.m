%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function plots a rasterplot of the neurons
% The main function of this program is piriformmain.m
%
% Licurgo de Almeida
% 03/15/2011
%
% The parameters here are:
% - Cell: Group of cell where we're taking the activity;
% - param: set of network parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotRasterplot(Cell,param)

scrsz = get(0,'ScreenSize');
figH = figure;
set(figH,'position',[0,400,scrsz(3)-0.4*scrsz(3),scrsz(4)-0.6*scrsz(4)]);
cla;
hold on;
figtitle = [Cell{1}.label,' cells'];
title(figtitle,'fontsize',16);
for ii = 1:length(Cell)
    J = find(Cell{ii}.S);
    for jj = 1:length(J)
        spkx = [J(jj),J(jj)] .* param.dt;
        spky = [ii,ii + 0.9];
        line(spkx,spky,'color','k','LineWidth',1);
    end
end

axis([0,param.tsim + param.dt,0,length(Cell) + 2]);
xlabel('time (ms)','fontsize',14);
ylabel('neuron','fontsize',14);