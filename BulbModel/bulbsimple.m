%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% bulbsimple.m is a simplified version of bulbmain.m and only produces
% Mitral cell output
%
% Licurgo de Almeida
% 07/12/2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starting program

function [Mitral param] = bulbsimple(mitfilename,odorant,test)

% Set new rand seed.
s = RandStream.create('mt19937ar','seed',sum(100*clock));
RandStream.setDefaultStream(s);

% Echo data
echodata = false;

% Input parameters
param.outputPath = 'cholinmod/';
param.dt = 0.5;
param.tsim = 4000;
param.tinit = 1;
param.tfinal = 4000;
param.nMitral = 50;
param.flagsavemitral = true;
param.mitfilename = mitfilename;
param.Odorant = odorant;
param.Test = test;
param.refrac = 20;

param.f = 22; 
param.delta = 6.9;
param.A = 9.2;
param.r = 25;



% Create cells
Mitral = cell(param.nMitral,1);

inputmat = CreateAffinity(param);

for ii = 1:param.nMitral
    Mitral{ii}.label = 'Mitral';
    Mitral{ii}.input = inputmat(ii);
end

% Create Mitral activity

time = (param.dt:param.dt:param.tsim) / 1000;
R = ((1 - cos(2 * pi * (param.f) * time)) / 2).^param.delta;
Y = R / (sum(R) / (param.tsim / 1000));
Y = Y * param.r;
A = zeros(param.nMitral,length(time));
Z = rand(size(A));
for ii = 1:param.nMitral
    A(ii,:) = Y * Mitral{ii}.input;
end
Miraster = A >= Z;
Miraster = Miraster * 1.0;
time = round(param.refrac / param.dt);
T = size(Miraster,2);

for ii = 1:T
    I = Miraster(:,ii) == 1;
    if ii + time <= T
        Miraster(I,ii + 1:ii + time) = 0;
    else
        Miraster(I,ii + 1:end) = 0;
    end
end


for ii = 1:param.nMitral
    Mitral{ii}.S = Miraster(ii,:);
end


if param.flagsavemitral == true
    SaveMitralFile(Miraster,param)
end

if echodata == true
    S = OdorSparseness(Mitral)
    F = mean(ExtractActivity(Mitral,param))
    C = GetCoherence(Mitral,'n',false,true)
end
end

function inputmat = CreateAffinity(param)

x = 1:param.nMitral;
y = normpdf(x,round(param.nMitral / 2),param.A);
y = y / max(y);
inputmat = SetAffinity(y,param.Odorant,param.nMitral);

end

function Input = SetAffinity(V,C,N)
V = V';
if C == round(N / 2)
    Input = V;
else
    Input = circshift(V,(C - round(N / 2)));
end
end

function SaveMitralFile(Miraster,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function saves the spike data from mitral cells, so we can use it on
% the piriform cortex program.
% The main function of this program is neurogenesismain.m
%
% Licurgo de Almeida
% 02/21/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MitralMat = Miraster;

cd ..

sfile = strcat('PiriformModel/',param.outputPath,param.mitfilename,num2str(param.Odorant,'%02.0f'),num2str(param.Test,'%02.0f'));
save(sfile,'MitralMat');

cd BulbModel
end