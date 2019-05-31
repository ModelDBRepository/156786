%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% piriformmain.m is the main function of this program. It reads data from
% a input and call all other functions
%
% Licurgo de Almeida
% 02/25/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Functions help
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To plot the activity distribution use PlotActivityBar(TypeofCell,true)
% To plot the activity of a given cell (or cells) and their respective
% inputs use PlotInputOutput(CellType,Inputs,Cellnumber,param)
% To plot a rasterplot use the PlotRasterplot(Cell,param)
% To plot an histogram with the number of pre-sinaptic spikes per
% postsinaptic spikes use HistPreSpikes(CellV,CellS,param,Pretime)
% To run piriformmain in a batch of tests use RunManyTests.m
% To test sparseness use S = OdorSparseness(Cell);
% To implement learning use [W,Conn] =
% ChangeWeights(Cell,param,plotweights)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Set new rand seed.
% s = RandStream.create('mt19937ar','seed',sum(100*clock));
% RandStream.setDefaultStream(s);

function [Mitral,Feedforward,Pyramidal,Feedback,param] = piriformmain(MitFile,PreSet)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reading parameters from input file and creating neurons
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(MitFile(1:5))
    case 'achon'
        param.inputFile = 'cholincortexON.txt';
    case 'achof'
        param.inputFile = 'cholincortexOFF.txt';
end

% Open the input file for reading
try
    fid1 = fopen(param.inputFile,'r');
    if fid1 == -1
        msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
        return;
    end
catch
    msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
    return;
end

str = fgetl(fid1);

% Get network parameters.
while ~strcmpi(str,'neurons')
    switch lower(str)
        case '' % It's possible to use blank lines to organize the
            % network parameters
        otherwise
            param = SetNetworkParameters(param,str);
    end
    
    str = fgetl(fid1);
end

% These parameters are changed outside the input file.
param.Mitralsource = MitFile;
param.flagPreSetConn = PreSet;

% Create cells
[Feedforward Pyramidal Feedback] = CreateCells(param);

str = fgetl(fid1);

% Get cell parameters.
while ~strcmpi(str,'end')
    
    switch lower(str)
        case '' % It's possible to use blank lines to organize the
            % neuronal parameters
        case 'feedforward'
            celltype = 'feedforward';
        case 'pyramidal'
            celltype = 'pyramidal';
        case 'feedback'
            celltype = 'feedback';
        otherwise
            switch celltype
                case 'feedforward'
                    Feedforward = SetNeuronParameters(Feedforward,param.nPyramidal,str);
                case 'pyramidal'
                    Pyramidal = SetNeuronParameters(Pyramidal,param.nPyramidal,str);
                case 'feedback'
                    Feedback = SetNeuronParameters(Feedback,param.nFeedback,str);
            end
    end
    
    str = fgetl(fid1);
end

fclose(fid1); % Close input file


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reading Mitral cell data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mitral = ReadMitralData(param);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Neuronal activity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Feedforward Pyramidal Feedback param] = NeuroActivity(Feedforward,Pyramidal,Feedback,Mitral,param);

if nargout < 5
    SaveFiles(Feedforward,Pyramidal,Feedback,param);
end
toc;

if param.savePyr == true
    SavePyramidalFile(Pyramidal,MitFile);
end

end

function param = SetNetworkParameters(param,str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the parameters for the different neurons
% The main function of this program is piriformmain.m
%
% Licurgo de Almeida
% 02/25/2011
% Information not related with the parameters of different neuros.
% * Path = path where we save input and output file
% * dt = timestep (in ms)
% * tsim = simulation time (in ms)
% * nPyramidal = number of pyramindal cells (and feedforward inhibitory neurons)
% * nFeedback = number of feedback inhibitory cells
% * CChanceMitPyr = Chance of connection between Mitral cells and Pyramidal
% cells.
% * CChanceMitFfo = Chance of connection between Mitral cells and
% Feedforward cells.
% * CChancePyrPyr = Chance of connection between two different Pyramidal
% cells.
% * CChancePyrFba = Chance of connection between Pyramidal and Feedback
% cells.
% * CChanceFbaPyr = Chance of connection between Feedback and Pyramidal
% cells.
% * CChanceFfoPyr = Chance of connection between Feedforward cells and Pyramidal
% cells.
% * NoiseFfo = true if we have noise in the Feedforward neurons
% * NoisePyr = true if we have noise in the Pyramidal cells
% * NoiseFba = true if we have noise in the Feedback cells
% * PreSetConn = true if we want to use always the same set of connections
% * Mitralsource = source of the Pyramidal inputs.
% * NoiseParam = true if we want a slighly variation on the parameters
% * Noiselevel = percetage of variation in each parameter
% * savePyr = true if we want to save Pyramidal cells for the HDB program
% * flagAHP = if true, AHP is ON
% * spikeV = Spike potential
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str_aux = 1;
% Find parameter name
while str(str_aux) ~= ' '
    str_aux = str_aux + 1;
end

ParName = str(1:str_aux - 1); % name of the parameter
ParValue = str(str_aux + 1:end); % value of the parameter

switch lower(ParName)
    case 'path'
        param.outputPath = ParValue; %path
    case 'dt'
        param.dt = str2num(ParValue); %step
    case 'tsim'
        param.tsim = str2num(ParValue); %simulation time
    case 'npyramidal'
        param.nPyramidal = str2num(ParValue); %number of Pyramidal cells (and Feedforward)
    case 'nfeedback'
        param.nFeedback = str2num(ParValue); %number of Feedback cells
    case 'cchancemitpyr'
        param.CChanceMitPyr = str2num(ParValue); %chance of connection between Mitral and Pyramidal cells
    case 'cchancemitffo'
        param.CChanceMitFfo = str2num(ParValue); %chance of connection between Mitral and Feedforward cells
    case 'cchancepyrpyr'
        param.CChancePyrPyr = str2num(ParValue); %chance of connection between two Pyramidal cells
    case 'cchancepyrfba'
        param.CChancePyrFba = str2num(ParValue); %chance of connection between Pyramidal and Feedback cells
    case 'cchancefbapyr'
        param.CChanceFbaPyr = str2num(ParValue); %chance of connection between Feedback and Pyramidal cells
    case 'cchanceffopyr'
        param.CChanceFfoPyr = str2num(ParValue); %chance of connection between Feedforward and Pyramidal cells
    case 'noiseffo'
        param.flagnoiseffo = str2num(ParValue); %flag noise Feedforward neurons
    case 'noisepyr'
        param.flagnoisepyr = str2num(ParValue); %flag noise Pyramidal cells
    case 'noisefba'
        param.flagnoisefba = str2num(ParValue); %flag noise Feedback neurons
    case 'presetconn'
        param.flagPreSetConn = str2num(ParValue); %flag preset connections
    case 'mitralsource'
        param.Mitralsource = ParValue; %name of the file source for Pyramidal information
    case 'noiseparam'
        param.NoiseParam = str2num(ParValue); %flag noise on the parameters
    case 'noiselevel'
        param.NoiseLevel = str2num(ParValue); %noise value
    case 'savepyr'
        param.savePyr = str2num(ParValue); %save pyramidal data
    case 'flagahp'
        param.flagAHP = str2num(ParValue); %flag use of AHP
    case 'spikev'
        param.SpikeV = str2num(ParValue); %Spike voltage
        
        % New parameters must be added here...
        
    otherwise
        disp(['parameter ' ParName ' does not exist']);
        
end
end

function [Feedforward Pyramidal Feedback] = CreateCells(param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates neurons based on the input file.
%
% Licurgo de Almeida
% 02/25/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pyramidal = cell(param.nPyramidal,1);
Feedforward = cell(param.nPyramidal,1);

for ii = 1:param.nPyramidal
    Pyramidal{ii}.input = []; %no input for now
    Pyramidal{ii}.label = 'Pyramidal';
    Feedforward{ii}.input = [];
    Feedforward{ii}.label = 'Feedforward'; %the number of pyramidal and
    % feedforward cells is always the same
end

Feedback = cell(param.nFeedback,1);

for ii = 1:param.nFeedback
    Feedback{ii}.input = []; %no input for now
    Feedback{ii}.label = 'Feedback';
end
end

function N = SetNeuronParameters(N,ncells,str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the parameters for the different neurons
% The main function of this program is piriform.m
%
% Licurgo de Almeida
% 02/28/2011
% Neurons cah present the following set of parameters:
% * tau = charging time constant of the neuron (ms). ms (and nos seconds) is
% the basic time unit in this program.
% * tauAHP = tau of the AHP potential.
% * Fthresh = Average firing threshold (in V)
% * Vrest = resting potential
% * Vhyper = hyperpolarization potential
% * EAMPA = AMPA's Nernst potential.
% * EGABA = GABA's Nernst potential.
% * EAHP = AHP's Nernst potential.
% * tauAMPA1 = AMPA's rising tau.
% * tauAMPA2 = AMPA's falling tau.
% * tauGABA1 = GABA's rising tau.
% * tauGABA2 = GABA's falling tau.
% * tauCA1 = CA's rising tau (used for the AHP).
% * tauCA2 = CA's falling tau (used for AHP).
% * wAMPA = excitatory synaptic weight
% * wAMPAPY = excitatory synaptic weight for Pyramidal association fibers
% (for Pyramidal cells only)
% * wGABA = inhibitory synaptic weight
% * wGABAFF = inhibitory synaptic weight from Feedforward cells (for
% Pyramidal cells only)
% * wAHP = AHP synaptic weight used to calculate AHP activation curves.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str_aux = 1;
% Find parameter name
while str(str_aux) ~= ' '
    str_aux = str_aux + 1;
end

ParName = str(1:str_aux - 1); % name of the parameter
ParValue = str(str_aux + 1:end); % value of the parameter

switch lower(ParName)
    case 'tau'
        for ii = 1:ncells
            N{ii}.tau = str2num(ParValue); %time in ms
        end
    case 'tauahp'
        for ii = 1:ncells
            N{ii}.tauAHP = str2num(ParValue); %time in ms
        end
    case 'fthresh'
        for ii = 1:ncells
            N{ii}.FThresh = str2num(ParValue); %potential in volts
        end
    case 'vrest'
        for ii = 1:ncells
            N{ii}.Vrest = str2num(ParValue); %potential in volts
        end
    case 'vhyper'
        for ii = 1:ncells
            N{ii}.Vhyper = str2num(ParValue); %potential in volts
        end
    case 'eampa'
        for ii = 1:ncells
            N{ii}.EAMPA = str2num(ParValue); %potential in volts
        end
    case 'egaba'
        for ii = 1:ncells
            N{ii}.EGABA = str2num(ParValue); %potential in volts
        end
    case 'eahp'
        for ii = 1:ncells
            N{ii}.EAHP = str2num(ParValue); %potential in volts
        end
    case 'tauampa1'
        for ii = 1:ncells
            N{ii}.tauAMPA1 = str2num(ParValue); %time in ms
        end
    case 'tauampa2'
        for ii = 1:ncells
            N{ii}.tauAMPA2 = str2num(ParValue); %time in ms
        end
    case 'taugaba1'
        for ii = 1:ncells
            N{ii}.tauGABA1 = str2num(ParValue); %time in ms
        end
    case 'taugaba2'
        for ii = 1:ncells
            N{ii}.tauGABA2 = str2num(ParValue); %time in ms
        end
    case 'tauca1'
        for ii = 1:ncells
            N{ii}.tauCA1 = str2num(ParValue); %time in ms
        end
    case 'tauca2'
        for ii = 1:ncells
            N{ii}.tauCA2 = str2num(ParValue); %time in ms
        end
    case 'wampa'
        for ii = 1:ncells
            N{ii}.wAMPA = str2num(ParValue); % excitatory synaptic weight
        end
    case 'wampapy'
        for ii = 1:ncells
            N{ii}.wAMPAPY = str2num(ParValue); % excitatory synaptic weight
            % from associative connections (pyramidal cells only)
        end
    case 'wgaba'
        for ii = 1:ncells
            N{ii}.wGABA = str2num(ParValue); % inhibitory synaptic weight
        end
    case 'wgabaff'
        for ii = 1:ncells
            N{ii}.wGABAFF = str2num(ParValue); % inhibitory synaptic weight
            % from feedforward cells (pyramidal cells only
        end
    case 'wahp'
        for ii = 1:ncells
            N{ii}.wAHP = str2num(ParValue); % AHP amplitude
        end
        
        % New parameters must be added here...
        
    otherwise
        disp(['parameter ' ParName ' does not exist']);
        
end
end

function Mitral = ReadMitralData(param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function reads the mitral cell data and creates a mitral cell group
% The main function of this program is piriformmain.m
%
% Licurgo de Almeida
% 03/10/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mitfile = strcat(param.outputPath,param.Mitralsource);

load(Mitfile,'MitralMat');

Mitral = cell(size(MitralMat,1),1);

for ii = 1:length(Mitral)
    Mitral{ii}.S = MitralMat(ii,:);
    Mitral{ii}.label = 'Mitral';
end
end

function [Feedforward Pyramidal Feedback param] = NeuroActivity(Feedforward,Pyramidal,Feedback,Mitral,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates the neuronal activity
% The main function of this program is piriformmain.m
%
% Licurgo de Almeida
% 02/28/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create mitral activity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MitralMat = zeros(length(Mitral),length(Mitral{1}.S));
for ii = 1:length(Mitral)
    MitralMat(ii,:) = Mitral{ii}.S;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set connections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if param.flagPreSetConn  == true
    Savefile = strcat(param.outputPath,'ConnData');
    load(Savefile,'MatMitPyr','MatPyrFba','MatFbaPyr','MatFfoPyr','MatPyrPyr','MatMitFfo');
    Savefile = strcat(param.outputPath,'LearnData'); %pyr-pyr learned weights
    load(Savefile,'W','Wsource');
    
else
    typecon = 'other'; % the type of connection is different for different pairs of neurons
    MatMitPyr = SetConnections(param.nPyramidal,size(MitralMat,1),param.CChanceMitPyr,typecon);
    MatMitFfo = SetConnections(param.nPyramidal,size(MitralMat,1),param.CChanceMitFfo,typecon);
    MatPyrFba = SetConnections(param.nFeedback,param.nPyramidal,param.CChancePyrFba,typecon);
    MatFfoPyr = SetConnections(param.nPyramidal,param.nPyramidal,param.CChanceFfoPyr,typecon);
    typecon = 'FbaPyr';
    MatFbaPyr = MatPyrFba';
    % The 5th parameter in SetConnections is only used when typecon = 'FbaPyr'
    MatFbaPyr = SetConnections(param.nPyramidal,param.nFeedback,param.CChanceFbaPyr,typecon,MatFbaPyr);
    typecon = 'PyrPyr';
    MatPyrPyr = SetConnections(param.nPyramidal,param.nPyramidal,param.CChancePyrPyr,typecon);
    W = zeros(param.nPyramidal,param.nPyramidal);
    Wsource = param.Mitralsource(1:5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set weights matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wGABAPyr = zeros(param.nPyramidal,1);
wAMPAPyr = wGABAPyr;
wAMPAPYPyr = wGABAPyr;
wGABAFFPyr = wGABAPyr;
wAMPAFfo = wGABAPyr;
for ii = 1:param.nPyramidal
    wGABAPyr(ii) = Pyramidal{ii}.wGABA;
    wAMPAPyr(ii) = Pyramidal{ii}.wAMPA;
    wAMPAPYPyr(ii) = Pyramidal{ii}.wAMPAPY;
    wGABAFFPyr(ii) = Pyramidal{ii}.wGABAFF;
    wAMPAFfo(ii) = Feedforward{ii}.wAMPA;
end
wFbaPyr = SetWeights(MatFbaPyr,wGABAPyr,param);
wMitPyr = SetWeights(MatMitPyr,wAMPAPyr,param);
wPyrPyr = SetWeights(MatPyrPyr,wAMPAPYPyr,param);
wFfoPyr = SetWeights(MatFfoPyr,wGABAFFPyr,param);

wMitFfo = SetWeights(MatMitFfo,wAMPAFfo,param);

wAMPAFba = zeros(param.nFeedback,1);
for ii = 1:param.nFeedback
    wAMPAFba(ii) = Feedback{ii}.wAMPA;
end
wPyrFba = SetWeights(MatPyrFba,wAMPAFba,param);

flag_normalize = false; % flag_normalize is true when the learning rule has
% no unlearning and we try to store many memories.

if flag_normalize == false
    wPyrPyr = wPyrPyr .* W;
else
    X = reshape(W,1,size(W,1) * size(W,2));
    I = X == 0;
    X(I) = [];
    X = sort(X,'descend');
    if strcmp(Wsource,'AChOn')
        wPyrPyr = wPyrPyr .* (W / mean(X(1:round(length(X) * 0.04))));
    elseif strcmp(Wsource,'AChOf')
        wPyrPyr = wPyrPyr .* (W / (mean(X(1:round(length(X) * 0.04))) / 0.56)); %correction factor
    else
        wPyrPyr = wPyrPyr .* 0;
    end
end


if strcmp(param.Mitralsource(1:5),'AChOn')
    wPyrPyr = wPyrPyr .* 0; %if ACho is on, these connections are inactive
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Feedforward neurons parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stores the voltage of each Feedforward cell at a given time
Vffo = zeros(param.nPyramidal,round(param.tsim / param.dt));
% Binary matrix recording only the spikes
Sffo = Vffo;
Cmitffo = Vffo;

refracffo = 30; % Refractory period after firing (ms)
countrefracffo = zeros(param.nPyramidal,1);

Restffo = ExtractParameter(Feedforward,'Vrest',param);
Hyperffo = ExtractParameter(Feedforward,'Vhyper',param);
tauffo = ExtractParameter(Feedforward,'tau',param);
Threshffo = ExtractParameter(Feedforward,'FThresh',param);
tauAMPA1ffo = ExtractParameter(Feedforward,'tauAMPA1',param);
tauAMPA2ffo = ExtractParameter(Feedforward,'tauAMPA2',param);
EAMPAffo = ExtractParameter(Feedforward,'EAMPA',param);

% ...New parameters should be added here...

% Initialize Feedforward cells potentials
Vffo(:,1) = Restffo;
% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tAMPA0ffo = zeros(size(MitralMat,1),1) - 10000000;
Imitffo = zeros(param.nPyramidal,1); % Input coming from Mitral cells
maxgAMPAffo = getmaxg(param.dt,tauAMPA1ffo,tauAMPA2ffo); % Get max conductance
% amplitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Feedforward neurons parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Pyramidal cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stores the voltage of each Pyramidal cell at a given time
Vpyr = zeros(param.nPyramidal,round(param.tsim / param.dt));
% Binary matrix recording only the spikes
Spyr = Vpyr;
Cmitpyr = Vpyr;
Cpyrpyr = Vpyr;
Cfbapyr = Vpyr;
Cffopyr = Vpyr;

% Stores the AHP values over time
VAHPpyr = Vpyr;

refracpyr = 2; % Refractory period after firing (ms)
countrefracpyr = zeros(param.nPyramidal,1);

Restpyr = ExtractParameter(Pyramidal,'Vrest',param);
Hyperpyr = ExtractParameter(Pyramidal,'Vhyper',param);
taupyr = ExtractParameter(Pyramidal,'tau',param);
Threshpyr = ExtractParameter(Pyramidal,'FThresh',param);
tauAMPA1pyr = ExtractParameter(Pyramidal,'tauAMPA1',param);
tauAMPA2pyr = ExtractParameter(Pyramidal,'tauAMPA2',param);
EAMPApyr = ExtractParameter(Pyramidal,'EAMPA',param);
tauGABA1pyr = ExtractParameter(Pyramidal,'tauGABA1',param);
tauGABA2pyr = ExtractParameter(Pyramidal,'tauGABA2',param);
EGABApyr = ExtractParameter(Pyramidal,'EGABA',param);
tauAHPpyr = ExtractParameter(Pyramidal,'tauAHP',param);
tauCA1pyr = ExtractParameter(Pyramidal,'tauCA1',param);
tauCA2pyr = ExtractParameter(Pyramidal,'tauCA2',param);
EAHPpyr = ExtractParameter(Pyramidal,'EAHP',param);

% ...New parameters should be added here...

% Initialize Pyramidal cells potentials
Vpyr(:,1) = Restpyr;

% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tAMPA0pyr = zeros(size(MitralMat,1),1) - 10000000;
Imitpyr = zeros(param.nPyramidal,1); % Input coming from Mitral cells
maxgAMPApyr = getmaxg(param.dt,tauAMPA1pyr,tauAMPA2pyr); % Get max conductance
% amplitude

% AMPAPY time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tAMPAPY0pyr = zeros(param.nPyramidal,1) - 10000000;
Ipyrpyr = zeros(param.nPyramidal,1); % Input coming from Pyramidal cells

% GABA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. Each Pyramidal cell can be
% connected to a varied number of Feedback cells
tGABA0pyr = zeros(param.nFeedback,1) - 10000000;
Ifbapyr = zeros(param.nPyramidal,1); % Input coming from Feedback cells
maxgGABApyr = getmaxg(param.dt,tauGABA1pyr,tauGABA2pyr); % Get max conductance
% amplitude

% GABAFF time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. Each Pyramidal cell can be
% connected to one Feedforward cell
tGABAFF0pyr = zeros(param.nPyramidal,1) - 10000000;
Iffopyr = zeros(param.nPyramidal,1); % Input coming from Feedforward cells

% AHP time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. The AHP in Pyramidal cells
% depends on their own fire.
tAHP0pyr = zeros(param.nPyramidal,1) - 10000000;
maxgAHPpyr = getmaxg(param.dt,tauCA1pyr,tauCA2pyr); % Get max conductance
% amplitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Pyramidal cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Feedback cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stores the voltage of each Feedback cell at a given time
Vfba = zeros(param.nFeedback,round(param.tsim / param.dt));
% Binary matrix recording only the spikes
Sfba = Vfba;
Cpyrfba = Vfba;

refracfba = 2; % Refractory period after firing (ms)
countrefracfba = zeros(param.nFeedback,1);

Restfba = ExtractParameter(Feedback,'Vrest',param);
Hyperfba = ExtractParameter(Feedback,'Vhyper',param);
taufba = ExtractParameter(Feedback,'tau',param);
Threshfba = ExtractParameter(Feedback,'FThresh',param);
tauAMPA1fba = ExtractParameter(Feedback,'tauAMPA1',param);
tauAMPA2fba = ExtractParameter(Feedback,'tauAMPA2',param);
EAMPAfba = ExtractParameter(Feedback,'EAMPA',param);

% ...New parameters should be added here...

% Initialize Feedback cells potentials
Vfba(:,1) = Restfba;

% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tAMPA0fba = zeros(param.nPyramidal,1) - 10000000;
Ipyrfba = zeros(param.nFeedback,1); % Input coming from Pyramidal cells
maxgAMPAfba = getmaxg(param.dt,tauAMPA1fba,tauAMPA2fba); % Get max conductance
% amplitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Pyramidal cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin neuron simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start loop
for tt = 2:round(param.tsim / param.dt)
    t = tt * param.dt; % current time
    
    % Mitral cell activity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I = MitralMat(:,tt - 1) == 1;
    
    % Pyramidal cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tAMPA0pyr(I) = t; % new t0 for the Pyramidal AMPA current
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Feedforward cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tAMPA0ffo(I) = t; % new t0 for the Pyramidal AMPA current
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Feedforward Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get Mitral input to Feedforward cells
    Imitffo(:) = SetI(tauAMPA1ffo(1),tauAMPA2ffo(1),t,tAMPA0ffo,maxgAMPAffo(1),...
        wMitFfo,EAMPAffo(1),Vffo(:,tt - 1));

    % Feedforward cell potential
    Vffo(:,tt) = Vffo(:,tt - 1) + (param.dt ./ tauffo(:)) .* (Imitffo(:) - ...
        Vffo(:,tt - 1) + Restffo(:));
    
    Cmitffo(:,tt) = Imitffo(:);
    
    
    % If the neuron fired last cycle, neuron potential hyperpotentializes
    I = Vffo(:,tt - 1) == param.SpikeV;
    Vffo(I,tt) = Hyperffo(I);
    
    % Pyramidal cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tGABAFF0pyr(I) = t; % new t0 for the Pyramidal cell GABAFF current
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    countrefracffo(I) = refracffo / param.dt; % neurons that just fired get into the
    % refractory period
    I = countrefracffo > 0;
    countrefracffo(I) = countrefracffo(I) - 1;
    
    I = find(countrefracffo == 0); % if countrefracffo = 0 the neuron can fire
    % again
    if ~isempty(I)
        if param.flagnoiseffo == true
            aux_J = SpikeNoise(Restffo(I),Threshffo(I),param,Vffo(I,tt),'Feedforward');
            J = find(aux_J);
        else
            J = find(Vffo(I,tt) >= Threshffo(I));
        end
        if ~isempty(J)
            Vffo(I(J),tt) = param.SpikeV; % Action potential
            Sffo(I(J),tt) = 1; % Record spike time
        end
    end
    
    % Pyramidal Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get Mitral input to Pyramidal cells
    Imitpyr(:) = SetI(tauAMPA1pyr(1),tauAMPA2pyr(1),t,tAMPA0pyr,maxgAMPApyr(1),...
        wMitPyr,EAMPApyr(1),Vpyr(:,tt - 1));
    
    % Get Pyramidal input to Pyramidal cells
    Ipyrpyr(:) = SetI(tauAMPA1pyr(1),tauAMPA2pyr(1),t,tAMPAPY0pyr,maxgAMPApyr(1),...
        wPyrPyr,EAMPApyr(1),Vpyr(:,tt - 1));
    
    % Get Feedback input to Pyramidal cells
    Ifbapyr(:) = SetI(tauGABA1pyr(1),tauGABA2pyr(1),t,tGABA0pyr,maxgGABApyr(1),...
        wFbaPyr,EGABApyr(1),Vpyr(:,tt - 1));
    
    % Get Feedback input to Pyramidal cells
    Iffopyr(:) = SetI(tauGABA1pyr(1),tauGABA2pyr(1),t,tGABAFF0pyr,maxgGABApyr(1),...
        wFfoPyr,EGABApyr(1),Vpyr(:,tt - 1));
    
    if param.flagAHP == true
        VAHPpyr(:,tt) = VAHPpyr(:,tt - 1) + (param.dt / tauAHPpyr(1)) *...
            (EAHPpyr(1) * Spyr(:,tt -1) - VAHPpyr(:,tt - 1));
    end

    % Feedforward cell potential
    Vpyr(:,tt) = Vpyr(:,tt - 1) + (param.dt ./ taupyr(:)) .* ((Imitpyr(:)...
         + Ipyrpyr(:) + Ifbapyr(:) + Iffopyr(:)) - Vpyr(:,tt - 1)...
         + Restpyr(:) + VAHPpyr(:,tt));
%      Vpyr(:,tt) = Vpyr(:,tt - 1) + (param.dt ./ taupyr(:)) .* (30e-3 ...
%          - Vpyr(:,tt - 1) + Restpyr(:) + VAHPpyr(:,tt));
     
     Cmitpyr(:,tt) = Imitpyr(:);
     Cpyrpyr(:,tt) = Ipyrpyr(:);
     Cfbapyr(:,tt) = Ifbapyr(:);
     Cffopyr(:,tt) = Iffopyr(:);
    
    % If the neuron fired last cycle, neuron potential hyperpotentializes
    I = Vpyr(:,tt - 1) == param.SpikeV;
    Vpyr(I,tt) = Hyperpyr(I);
    
    % Pyramidal cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tAMPAPY0pyr(I) = t; % new t0 for the Pyramidal cell AMPA current
    tAHP0pyr(I) = t; % new t0 for the Pyramidal cell AHP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Feedback cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tAMPA0fba(I) = t; % new t0 for the Feedback cell AMPA current
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    countrefracpyr(I) = refracpyr / param.dt; % neurons that just fired get into the
    % refractory period
    I = countrefracpyr > 0;
    countrefracpyr(I) = countrefracpyr(I) - 1;
    
    I = find(countrefracpyr == 0); % if countrefracpyr = 0 the neuron can fire
    % again
    if ~isempty(I)
        if param.flagnoisepyr == true
            aux_J = SpikeNoise(Restpyr(I),Threshpyr(I),param,Vpyr(I,tt),'Pyramidal');
            J = find(aux_J);
        else
            J = find(Vpyr(I,tt) >= Threshpyr(I));
        end
        if ~isempty(J)
            Vpyr(I(J),tt) = param.SpikeV; % Action potential
            Spyr(I(J),tt) = 1; % Record spike time
        end
    end
    
    
    % Feedback Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get Pyramidal input to Feedback cells
    Ipyrfba(:) = SetI(tauAMPA1fba(1),tauAMPA2fba(1),t,tAMPA0fba,maxgAMPAfba(1),...
        wPyrFba,EAMPAfba(1),Vfba(:,tt - 1));

    % Feedforward cell potential
    Vfba(:,tt) = Vfba(:,tt - 1) + (param.dt ./ taufba(:)) .*...
        (Ipyrfba(:) - Vfba(:,tt - 1) + Restfba(:));
    
    Cpyrfba(:,tt) = Ipyrfba(:);
    
    
    % If the neuron fired last cycle, neuron potential hyperpotentializes
    I = Vfba(:,tt - 1) == param.SpikeV;
    Vfba(I,tt) = Hyperfba(I);
    
    % Pyramidal cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tGABA0pyr(I) = t; % new t0 for the Pyramidal cell GABA current
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    countrefracfba(I) = refracfba / param.dt; % neurons that just fired get into the
    % refractory period
    I = countrefracfba > 0;
    countrefracfba(I) = countrefracfba(I) - 1;
    
    I = find(countrefracfba == 0); % if countrefracfba = 0 the neuron can fire
    % again
    if ~isempty(I)
        if param.flagnoisefba == true
            aux_J = SpikeNoise(Restfba(I),Threshfba(I),param,Vfba(I,tt),'Feedback');
            J = find(aux_J);
        else
            J = find(Vfba(I,tt) >= Threshfba(I));
        end
        if ~isempty(J)
            Vfba(I(J),tt) = param.SpikeV; % Action potential
            Sfba(I(J),tt) = 1; % Record spike time
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Feedforward and Pyramidal cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:param.nPyramidal
    Feedforward{ii}.V = Vffo(ii,:); % Save neuronal activity
    Feedforward{ii}.S = Sffo(ii,:); % Save spike time
    Feedforward{ii}.Cmit = Cmitffo(ii,:);
    Pyramidal{ii}.V = Vpyr(ii,:); % Save neuronal activity
    Pyramidal{ii}.AHP = VAHPpyr(ii,:); % Save AHP activity
    Pyramidal{ii}.S = Spyr(ii,:); % Save spike time
    Pyramidal{ii}.Cmit = Cmitpyr(ii,:);
    Pyramidal{ii}.Cffo = Cffopyr(ii,:);
    Pyramidal{ii}.Cpyr = Cpyrpyr(ii,:);
    Pyramidal{ii}.Cfba = Cfbapyr(ii,:);
    Feedforward{ii}.MitCon = MatMitFfo(ii,:); % Save Mitral connections
    Pyramidal{ii}.MitCon = MatMitPyr(ii,:); % Save Mitral connections
    Pyramidal{ii}.FfoCon = MatFfoPyr(ii,:); % Save Feedforward connections
    Pyramidal{ii}.FbaCon = MatFbaPyr(ii,:); % Save Feedback connections
    Pyramidal{ii}.PyrCon = MatPyrPyr(ii,:); % Save Pyramidal connections
end

% Feedback cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:param.nFeedback
    Feedback{ii}.V = Vfba(ii,:); % Save neuronal activity
    Feedback{ii}.S = Sfba(ii,:); % Save spike time
    Feedback{ii}.Cpyr = Cpyrfba(ii,:);
    Feedback{ii}.PyrCon = MatPyrFba(ii,:); % Save Pyramidal connections
end

% Connections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Savefile = strcat(param.outputPath,'ConnData');
save(Savefile,'MatMitPyr','MatPyrFba','MatFbaPyr','MatFfoPyr','MatPyrPyr','MatMitFfo');
Savefile = strcat(param.outputPath,'LearnData'); %pyr-pyr learned weights
save(Savefile,'W','Wsource');
end

function Mat = SetConnections(cell1,cell2,cchance,typecon,Maux)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the connection between different neurons
% The main function of this program is piriformmain.m
%
% Licurgo de Almeida
% 03/01/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randtype = 'fixed';

switch randtype
    case 'poisson'
        if strcmpi(typecon,'pyrpyr')
            Mat = rand(cell1,cell2);
            auxcchance = (cchance * cell2) / (cell2 - 1);
            Mat = Mat <= auxcchance;
            for ii = 1:cell1
                Mat(ii,ii) = 0;
            end
        elseif strcmpi(typecon,'fbapyr')
            Mat = rand(cell1,cell2);
            Mat = Mat <= cchance;
            Mat = Mat + Maux;
            Mat = Mat > 0;
        else
            Mat = rand(cell1,cell2);
            Mat = Mat <= cchance;
        end
    case 'fixed'
        if strcmpi(typecon,'pyrpyr')
            if cchance == 1
                Mat = ones(cell1,cell2);
                for ii = 1:cell1
                    Mat(ii,ii) = 0;
                end
            else
                Mat = zeros(cell1,cell2);
                for ii = 1:cell1
                    convector = randperm(cell2);
                    convector = convector(1:(round(cell2 * cchance) + 1 ...
                    + round(0.5 * randn)));
                    I = find(convector == ii);
                    if ~isempty(I)
                        convector(I) = [];
                    else
                        convector(1) = [];
                    end
                    Mat(ii,convector) = 1;
                end
            end
        elseif strcmpi(typecon,'fbapyr')
            Mat = zeros(cell1,cell2);
            nconnections = round(cell2 * cchance);
            for ii = 1:cell1
                convector = randperm(cell2);
                I = find(Maux(ii,:));
                for jj = 1:length(I)
                    J = convector == I(jj);
                    convector(J) = [];
                end
                if nconnections - length(I) > 0
                    convector = convector(1:(nconnections - length(I)));
                else
                    convector = convector(1:nconnections);
                end
                Mat(ii,convector) = 1;
                Mat(ii,I) = 1;
            end
        else
            Mat = zeros(cell1,cell2);
            for ii = 1:cell1
                convector = randperm(cell2);
                convector = convector(1:round(cell2 * cchance));
                Mat(ii,convector) = 1;
            end
        end
end
end

function w = SetWeights(Mat,weights,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the synaptic weights between different neurons
% The main function of this program is piriformmain.m
%
% Licurgo de Almeida
% 03/02/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = zeros(size(Mat));
for ii = 1:length(weights)
    w(ii,:) = Mat(ii,:) * weights(ii);
end

if param.NoiseParam == true
    for ii = 1:length(weights)
        signvar = rand(1,size(Mat,2));
        I = signvar < 0.5;
        signvar(I) = -1;
        I = signvar >= 0.5;
        signvar(I) = 1;
        
        changeval = (rand(1,size(Mat,2)) * param.NoiseLevel * weights(ii)) .* Mat(ii,:);
        w(ii,:) = w(ii,:) + (signvar .* changeval);
    end
end
end

function parameter = ExtractParameter(neuron,ParName,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function extract the parameters from the neuronfile
% The main function of this program is piriformmain.m
%
% Licurgo de Almeida
% 03/01/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameter = zeros(length(neuron),1);

switch lower(ParName)
    case 'vrest'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.Vrest;
        end
    case 'vhyper'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.Vhyper;
        end
    case 'r'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.R;
        end
    case 'tau'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tau;
        end
    case 'tauahp'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauAHP;
        end
    case 'fthresh'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.FThresh;
        end
    case 'gmaxampa'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.gmaxAMPA;
        end
    case 'gmaxampapy'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.gmaxAMPAPY;
        end
    case 'gmaxgaba'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.gmaxGABA;
        end
    case 'gmaxgabaff'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.gmaxGABAFF;
        end
    case 'tauampa1'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauAMPA1;
        end
    case 'tauampa2'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauAMPA2;
        end
    case 'taugaba1'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauGABA1;
        end
    case 'taugaba2'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauGABA2;
        end
    case 'tauca1'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauCA1;
        end
    case 'tauca2'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauCA2;
        end
    case 'eampa'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.EAMPA;
        end
    case 'egaba'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.EGABA;
        end
    case 'eahp'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.EAHP;
        end
end

if param.NoiseParam == true
    for ii = 1:length(parameter)
        signvar = rand;
        if signvar < 0.5
            signivar = -1;
        else
            signivar = 1;
        end
        parameter(ii) = parameter(ii) + (signivar * parameter(ii) * rand* param.NoiseLevel);
    end
end
end

function maxg = getmaxg(dt,tau1,tau2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function finds the max amplitude for g, so we can use this value to
% normalize the curve
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 11/10/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmax = 100;

x = 0:dt:tmax;
maxg = zeros(length(tau1),1);

for ii = 1: length(tau1)
    if ii == 1
        y = ((tau1(ii) * tau2(ii)) / (tau1(ii) - tau2(ii))) * (exp(-(x)...
            / tau1(ii)) - exp(-(x) / tau2(ii)));
        maxg(ii) = max(y);
    else
        if tau1(ii) == tau1(ii -1) && tau2(ii) == tau2(ii -1)
            maxg(ii) = maxg(ii - 1);
        else
            y = ((tau1(ii) * tau2(ii)) / (tau1(ii) - tau2(ii)))...
                * (exp(-(x) / tau1(ii)) - exp(-(x) / tau2(ii)));
            maxg(ii) = max(y);
        end
    end
end
end

function Ic = SetI(tau1,tau2,t,t0,normg,W,E,V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the currents at each time step.
% This function is called from NeuroActivity.m
%
% Licurgo de Almeida
% 11/18/2010
%
% Ic: channel current
% gmax: maximum conductance
% tau1: channel's rising time
% tau2: channel's falling time
% t: step time
% t0: last time the pre synaptic neuron fired
% normg: normalizes the conductance curve between 0 and 1
% W: synaptic weights
% E: Nernst potential
% V: neuron's potential from last timestep
% Mcon: connection matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


g = (((tau1 * tau2) / (tau1 - tau2)) *...
    (exp(-(t - t0) / tau1) - exp(-(t - t0) / tau2))) / normg;

Ic = (W * g) .* (E - V);

end

function I = SpikeNoise(Rest,Thres,param,V,tnet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the spiking chance
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 03/03/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(tnet,'OSN')
    pfactor = 3;
elseif strcmp(tnet,'Pyramidal')
    pfactor = 10;
else
    pfactor = 5;
end

Vlimit = 2e-3; % Limits where firing chances varied between 0 and 1
Rest = Rest - Vlimit;
Thres = Thres + Vlimit;
y = (V - Rest) ./ (Thres -  Rest);
J = y < 0;
y(J) = 0;
J = y > 1;
y(J) = 1;
y = y.^pfactor;
I = rand(length(Rest),1) <= y * param.dt;
end

function SaveFiles(Feedforward,Pyramidal,Feedback,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function save the cells in different files if we don't want to
% return anything back in the piriformmain.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = clock;

Y = num2str(T(1));
Y = Y(3:4);
Mo = num2str(T(2),'%02.0f');
D = num2str(T(3),'%02.0f');
H = num2str(T(4),'%02.0f');
Mi = num2str(T(5),'%02.0f');

if strcmp(param.Mitralsource(1:5),'AChOf') == true
    ACh = 'F';
else
    ACh = 'N';
end

outputfolder = 'AChout/';

filedate = strcat(Y,Mo,D,H,Mi,ACh);

outputFile = strcat(outputfolder,'Feedforward',filedate,'.txt');

fid2 = fopen(outputFile,'w');

for ii = 1:length(Feedforward)
    fprintf(fid2,'%1i',Feedforward{ii}.S);
    if ii < length(Feedforward)
        fprintf(fid2,'\n');
    end
end
fclose(fid2);

outputFile = strcat(outputfolder,'Pyramidal',filedate,'.txt');

fid2 = fopen(outputFile,'w');

for ii = 1:length(Pyramidal)
    fprintf(fid2,'%1i',Pyramidal{ii}.S);
    if ii < length(Pyramidal)
        fprintf(fid2,'\n');
    end
end
fclose(fid2);

outputFile = strcat(outputfolder,'Feedback',filedate,'.txt');

fid2 = fopen(outputFile,'w');

for ii = 1:length(Feedback)
    fprintf(fid2,'%1i',Feedback{ii}.S);
    if ii < length(Feedback)
        fprintf(fid2,'\n');
    end
end
fclose(fid2);

end

function SavePyramidalFile(Pyramidal,MitFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function saves Pyramidal cells, so we can use it on
% the HDB program.
%
% Licurgo de Almeida
% 02/02/2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..

sfile = strcat('HDBModel/Data/',MitFile);
save(sfile,'Pyramidal');

cd PiriformModel
end