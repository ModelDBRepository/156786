%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% bulbmain.m is the main function of this program. It reads data from
% a input and call all other functions
%
% Licurgo de Almeida
% 10/29/2010
%
% Other useful functions:
%
% * To test the sparseness in the network, use TestSparseness.m
% * To calculate correlation use [OSNsource Grasource] =
% getcorrelation(OSNsource,Grasource,Matrices). This function works only if
% ConnType = spatial.
% * To plot activity maps WHEN CONNTYPE = SPATIAL use Matrices =
% PlotMaps(OSNsource,Grasource,Matrices,OSN,Mitral,fignum).
% * To plot the activity of some type of cell use
% PlotActivity(Cell,param,CellLabel,fignum).
% * To plot the voltage of some type of cell use
% PlotVoltage(Cell,param,CellLabel,fignum).
% * To plot specific neurons from a network, use Rasterplot.m
% * To plot the circular maps use
% PlotCircularMaps(OSN,Pglo,Mitral,Granule,fignum,param). This function
% works when ConnType = circular.
% * To plot non-spiking Pglo cells, use PlotPglo(Pglo).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starting program

function [OSN Pglo Glo Mitral Granule OSNsource Grasource param] = bulbmain(inputFile,odorant,test)
tic;

if strcmpi(inputFile(end - 2:end),'txt') % if we use a txt as input, the
    % program reads the txt file and create the variables, if we use a mat file,
    % the program loads the variables.
    
    % Set new rand seed.
    s = RandStream.create('mt19937ar','seed',sum(100*clock));
    RandStream.setDefaultStream(s);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Reading parameters from input file and creating neurons
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param.inputFile = inputFile;
    
    
    % Open the input file for reading
    fid1 = fopen(param.inputFile,'r');
    if fid1 == -1
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
    
    param.mitfilename = strcat(param.mitfilename,num2str(odorant,'%02.0f'),num2str(test,'%02.0f'));
    param.Odorant = odorant;
    
    % Create cells
    if param.flagpreset == true
        if strcmpi(param.ConnType,'randperm')
            param.ConnType = 'spatialfix';
        end
        OSNfile = strcat(param.outputPath,param.OSNsource);
        if strcmp(param.ConnType(1:7),'spatial')
            Grafile = strcat(param.outputPath,param.Grasource);
        else
            Grafile = []; % if network is circular, grafile is not used.
        end
        [param,OSNsource,Grasource,OSN,Mitral,Pglo,Glo,Granule] = CreateCells(param,OSNfile,Grafile);
    else
        OSN = cell(param.nMitral,1);
        Mitral = cell(param.nMitral,1);
        Pglo = cell(param.nMitral,1);
        Glo = cell(param.nMitral,1);
        for ii = 1:param.nMitral
            OSN{ii}.input = []; %no input for now
            Mitral{ii}.input = [];
            Pglo{ii}.input = [];
            Glo{ii}.input = [];
            Mitral{ii}.Acorrec = 1;
        end
        
        Granule = cell(param.nGranule,1);
        for ii = 1:param.nGranule
            Granule{ii}.input = [];
            Granule{ii}.Acorrec = 1;
        end
    end
    
    str = fgetl(fid1);
    
    % Get cell parameters.
    while ~strcmpi(str,'end')
        
        switch lower(str)
            case '' % It's possible to use blank lines to organize the
                % neuronal parameters
            case 'osn'
                celltype = 'osn';
            case 'periglomerular'
                celltype = 'periglomerular';
            case 'glomerulus'
                celltype = 'glomerulus';
            case 'mitral'
                celltype = 'mitral';
            case 'granule'
                celltype = 'granule';
            otherwise
                switch celltype
                    case 'osn'
                        OSN = SetNeuronParameters(OSN,param.nMitral,str);
                    case 'periglomerular'
                        Pglo = SetNeuronParameters(Pglo,param.nMitral,str);
                    case 'glomerulus'
                        Glo = SetNeuronParameters(Glo,param.nMitral,str);
                    case 'mitral'
                        Mitral = SetNeuronParameters(Mitral,param.nMitral,str);
                    case 'granule'
                        Granule = SetNeuronParameters(Granule,param.nGranule,str);
                end
                
        end
        
        str = fgetl(fid1);
    end
    
    fclose(fid1); % Close input file
    fname = inputFile(1:end - 3);
    fname = strcat(fname,'mat');
    save(fname,'OSN','Pglo','Glo','Mitral','Granule','OSNsource','Grasource','param');
    
elseif strcmpi(inputFile(end - 2:end),'mat') %if the input file is .mat
    % we have to
    load(inputFile,'OSN','Pglo','Glo','Mitral','Granule','OSNsource','Grasource','param');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set cell position
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(param.ConnType(1:7),'spatial')
    % Show in terminal parameters that are not used for this simulation
    disp('Parameter NGraCon is not being used');
    if param.flagpreset == true
        disp('Using pre set data from source file');
        param.MitralRadius = Mitral{1}.CellRadius;
        param.GranuleRadius = Granule{1}.CellRadius;
        % Inside the PreSetCellPosition we also define the area-input
        % correction factor. Whe the shape of the bulb is not a toroid, we
        % have to multiply the inputs received by neurons on the borders so
        % that the overall activity to each neurons will be similar.
        Mitral = PreSetCellPosition(Mitral,param,OSNsource);
        Granule = PreSetCellPosition(Granule,param,OSNsource);
    else
        Mitral = SetCellPosition(Mitral,param);
        Granule = SetCellPosition(Granule,param);
    end
elseif strcmp(param.ConnType,'randperm')
    % Show in terminal parameters that are not used for this simulation
    disp('Parameters ConnChance, BulbH and BulbW are not being used');
    if param.flagpreset == true
        disp('NSO source file is not being used!');
    end
elseif strcmp(param.ConnType,'circular')
    Mitral = PreSetCellPosition(Mitral,param,OSNsource);
    Granule = PreSetCellPosition(Granule,param,Grasource);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Neuronal activity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[OSN Pglo Glo Mitral Granule param] = NeuroActivity(OSN,Pglo,Glo,Mitral,Granule,param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reduce granular matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(param.ConnType(1:7),'spatial')
    Matrices = CompactGraMat(param,Grasource,Granule);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Saving data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if param.flagsavemitral == true
    SaveMitralFile(Mitral,param);
end

toc;
end

function param = SetNetworkParameters(param,str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the parameters for the different neurons
%
% Licurgo de Almeida
% 12/20/2010
% Information not related with the parameters of different neuros.
% dt: timestep (in ms)
% tsim: simulation time (in ms)
% tinit: stimulus begin (in ms)
% tfinal: stimulus end (in ms)
% nMitral: number of mitral cells (and olfactory sensory neurons (OSN))
% nGranule: number of granule cells
% ConnType = type of connection between cells. Can be 'randperm' cells are
% randomly connected; 'spatialvar' the closer/father to each other the cells are
% the higher chance of being connected; 'spatialfix' cells that are close
% to each other have a fixed chance of being connected; 'spatialfile' uses
% the granule cell file to define the chance of connection between two
% neurons; 'circle' works as a 1d plot to test contrast;
% ConnChanceMitGra = in case the chance of connection is fixed, tell the
% exact percentage of chances of Mitral - Granule connection
% ConnChanceGraGra = in case the chance of connection is fixed, tell the
% exact percentage of chances of Granule - Granule connection
% GraGracon = if true, granule cells connect to each other
% SpikingInput = if true, OSNs are integrate and fire neurons
% SpikingPglo = if true, Pglos are integrate and fire neurons, otherwise 
% they're a logistic function.
% BulbH = Bulb height (in distance units)
% BulbW = Bulb width (in distance units)
% NoiseOSN = true if we have noise in the OSN
% NoisePgl = true if we have noise in the Periglomerular cells
% NoiseMit = true if we have noise in the Mitral cells
% NoiseGra = true if we have noise in the Granule cells
% PreSet = true if we use pre set position for the cells
% SaveMitral = true if we want to save the mitral cell activity for the
% piriform cortex model
% MitFileName = Mitral cell file name
% mFactor = multiplicative factor. number of granule cells = number of
% mitral cells * mFactor
% OSNsource = source of the OSN inputs.
% Odorant = odortant number.
% Grasource = source of the granule cells (used only for comparison
% (add the name of new parameters here)
% loadOSNdata = true if we want to use pre-processed data.
% toroid = true if we want to use a toroid shape, if false, we use a
% cilinder shape
% weights = if true, the synatic weight of some connectios is variable
% (using learning)
% AHP = if true, mitral cells show AHP
% Respiration = if true, the input is modulated by an oscillation
% representing the respiration
% RespFreq = Respiratory frequency
% LoadConn = true if we want to use pre-processed Mit-Gra connections
% Iext = external current to OSNs
% SpikeV = Spike voltage
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
    case 'tinit'
        param.tinit = str2num(ParValue); %stimulus begin
    case 'tfinal'
        param.tfinal = str2num(ParValue); %stimulus end
    case 'nmitral'
        param.nMitral = str2num(ParValue); %number of Mitral cells (and OSN)
    case 'ngranule'
        param.nGranule = str2num(ParValue); %number of Granule cells
    case 'conntype'
        param.ConnType = ParValue; %type of connection
    case 'connchancemitgra'
        param.ConnChanceMitGra = str2num(ParValue); %chance of connection between mit and gra
    case 'connchancegragra'
        param.ConnChanceGraGra = str2num(ParValue); %chance of connection between gra and gra
    case 'gragracon'
        param.GraGracon = str2num(ParValue); %Gra-Gra connections
    case 'spikinginput'
        param.SpikingInput = str2num(ParValue); %flag spiking OSNs
    case 'spikingpglo'
        param.SpikingPglo = str2num(ParValue); %flag spiking Pglo cells
    case 'bulbh'
        param.BulbH = str2num(ParValue); %Bulb height
    case 'bulbw'
        param.BulbW = str2num(ParValue); %Bulb width
    case 'noiseosn'
        param.flagnoiseosn = str2num(ParValue); %flag noise OSN
    case 'noisepgl'
        param.flagnoisepgl = str2num(ParValue); %flag noise Periglomerular
    case 'noisemit'
        param.flagnoisemit = str2num(ParValue); %flag noise Mitral
    case 'noisegra'
        param.flagnoisegra = str2num(ParValue); %flag noise granule
    case 'preset'
        param.flagpreset = str2num(ParValue); %flag preset positions
    case 'mfactor'
        param.mFactor = str2num(ParValue); % multiplicative factor
    case 'osnsource'
        param.OSNsource = ParValue; %name of the file source for OSN information
    case 'odorant'
        param.Odorant = str2num(ParValue); %odorant number
    case 'grasource'
        param.Grasource = ParValue; %name of the file source for granule cells information
    case 'loadosndata'
        param.flagOSNdata = str2num(ParValue); %flag load preset OSN data
    case 'toroid'
        param.flagtoroid = str2num(ParValue); %flag toroid shape
    case 'weights'
        param.flagweights = str2num(ParValue); %flag use of synaptic weights
    case 'savemitral'
        param.flagsavemitral = str2num(ParValue); %flag save mitral file
    case 'mitfilename'
        param.mitfilename = ParValue; % Mitral file name.
    case 'ahp'
        param.flagAHP = str2num(ParValue); %flag presence of AHP
    case 'respiration'
        param.flagRespiration = str2num(ParValue); %flag respiratory modulation
    case 'respfreq'
        param.RespFreq = str2num(ParValue); %respiratory frequency
    case 'loadconn'
        param.flagLoadConn = str2num(ParValue); %load Mit-Gra conn
    case 'iext'
        param.Iext = str2num(ParValue); %External current to OSNs
    case 'spikev'
        param.SpikeV = str2num(ParValue); %Spike voltage
        
        % New parameters must be added here...
        
    otherwise
        disp(['parameter ' ParName ' does not exist']);
        
end

end

function [param,OSNsource,Grasource,OSN,Mitral,Pglo,Glo,Granule] = CreateCells(param,OSNfile,Grafile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates neurons based on the OSN file source (excel file).
%
% Licurgo de Almeida
% 12/21/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(param.ConnType(1:7),'spatial')
    
    mfactor = round(sqrt(param.mFactor)); %multiplicative factor of columns and rows, the number of
    %granule cell
    
    param.mFactor = mfactor^2; %mFactor is corrected to the real value
    
    warning('off','MATLAB:xlsread:ActiveX');
    OSNsource.inputs = xlsread(OSNfile); %read excel file with NSO information
    Grasource.outputs = xlsread(Grafile); %read excel file with Granule cells information
    warning('on','MATLAB:xlsread:ActiveX')
    OSNsource.inputs = OSNsource.inputs / max(max(OSNsource.inputs)); %normalize data
    [X,Y] = size(OSNsource.inputs);
    Grasource.outputs = Grasource.outputs / max(max(Grasource.outputs)); %normaliza data
    disp(['Number of Mitral cells changed to ' num2str(X * Y)]);
    param.nMitral = X * Y;
    OSN = cell(param.nMitral,1);
    Mitral = cell(param.nMitral,1);
    Pglo = cell(param.nMitral,1);
    Glo = cell(param.nMitral,1);
    OSNsource.X = zeros(param.nMitral,1);
    OSNsource.Y = OSNsource.X;
    count = 1;
    for ii = 1:Y
        for jj = 1:X
            OSN{count}.input = OSNsource.inputs(jj,ii);
            Mitral{count}.input = []; %no input for now
            Pglo{count}.input = [];
            Glo{count}.input = [];
            OSNsource.X(count) = jj;
            OSNsource.Y(count) = ii;
            count = count + 1;
        end
    end
    
    disp(['Number of Granule cells changed to ' num2str((X * mfactor) * (Y * mfactor))]);
    param.nGranule = (X * mfactor) * (Y * mfactor);
    Granule = cell(param.nGranule,1);
    Grasource.X = zeros(param.nGranule,1);
    Grasource.Y = Grasource.X;

    if strcmp(param.ConnType,'spatialfile')
        auxMat = imresize(Grasource.outputs,mfactor); %resize matrix
    end
    count = 1;
    for ii = 1:(Y * mfactor)
        for jj = 1:(X * mfactor)
            Granule{count}.input = [];
            Grasource.X(count) = jj;
            Grasource.Y(count) = ii;
            if strcmp(param.ConnType,'spatialfile')
                Granule{count}.Connchance = auxMat(jj,ii);
            end
            count = count + 1;
        end
    end
    
    % Preset also changes the size of the network
    disp(['Bulb high changed to ' num2str(X) ' and Buld width changed to ' num2str(Y)]);
    param.BulbH = X;
    param.BulbW = Y;
    
elseif strcmp(param.ConnType,'circular')
    aux_inputs = load(OSNfile); %read file with NSO information
    %OSNsource.inputs = aux_inputs.inputmatrand(param.Odorant,:);
    OSNsource.inputs = aux_inputs.inputmat(param.Odorant,:);
    %OSNsource.inputs = OSNsource.inputs / max(OSNsource.inputs); %normalize data
    OSN = cell(param.nMitral,1);
    Mitral = cell(param.nMitral,1);
    Pglo = cell(param.nMitral,1);
    Glo = cell(param.nMitral,1);
    OSNsource.X = zeros(param.nMitral,1);
    if length(OSNsource.inputs) ~= param.nMitral
        msgbox('Input file and number of cells are different!','ERROR');
        return;
    end
    for ii = 1:param.nMitral
        OSN{ii}.input = OSNsource.inputs(ii);
        OSN{ii}.label = 'OSN';
        Mitral{ii}.input = []; %no input for now
        Mitral{ii}.label = 'Mitral';
        Pglo{ii}.input = [];
        Pglo{ii}.label = 'Pglo';
        Glo{ii}.input = [];
        Glo{ii}.label = 'Glo';
        OSNsource.X(ii) = ii;
    end
    Granule = cell(param.nGranule,1);
    Grasource.X = zeros(param.nGranule,1);
    x = param.nMitral / param.nGranule:param.nMitral / param.nGranule:param.nMitral;
    for ii = 1:param.nGranule
        Granule{ii}.input = [];
        Granule{ii}.label = 'Granule';
        Grasource.X(ii) = x(ii);
    end
end
end

function N = SetNeuronParameters(N,ncells,str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the parameters for the different neurons
%
% Licurgo de Almeida
% 11/02/2010
% Neurons cah present the following set of parameters:
% * tau = charging time constant of the neuron (ms). ms (and nos seconds) is
% the basic time unit in this program.
% * R = Membrane resistance (in ohms)
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
% * tauAHP1 = AHP's rising tau.
% * tauAHP2 = AHP's falling tau.
% * gmaxAMPA = AMPA's max conductance
% * gmaxGABA = GABA's max conductance
% * gmaxGABAP = Periglomerula GABA's max conductance (for Mitral cells
% only)
% * gmaxAHP = AHP's max conductance
% * IACh = Addition current when ACh is ON in non-spiking cells.
% * nGracon = number of granule cells connected to mitral or granule cells
% * CellRadius = radius of mitral and granule cells' dendrites
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
    case 'r'
        for ii = 1:ncells
            N{ii}.R = str2num(ParValue); %resistance in ohms
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
    case 'noise'
        for ii = 1:ncells
            N{ii}.Noise = str2num(ParValue);
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
    case 'tauahp1'
        for ii = 1:ncells
            N{ii}.tauAHP1 = str2num(ParValue); %time in ms
        end
    case 'tauahp2'
        for ii = 1:ncells
            N{ii}.tauAHP2 = str2num(ParValue); %time in ms
        end
    case 'gmaxampa'
        for ii = 1:ncells
            N{ii}.gmaxAMPA = str2num(ParValue); %AMPA channel
            % conductance in siemens
        end
    case 'gmaxgaba'
        for ii = 1:ncells
            N{ii}.gmaxGABA = str2num(ParValue); %GABA channel
            % conductance in siemens
        end
    case 'gmaxahp'
        for ii = 1:ncells
            N{ii}.gmaxAHP = str2num(ParValue); %"AHP" channel (for the simulation)
            % conductance in siemens
        end
    case 'iach'
        for ii = 1:ncells
            N{ii}.IACh = str2num(ParValue); %Additional current when ACh is ON
            % in non-spiking neurons
        end
    case 'gmaxgabap'
        for ii = 1:ncells
            N{ii}.gmaxGABAP = str2num(ParValue); %GABA channel from PG cells
            % conductance in siemens
        end
    case 'ngracon'
        for ii = 1:ncells
            N{ii}.NGraCon = str2num(ParValue); %number of granule cells
            % connected to a mitral or granule cell
        end
    case 'cellradius'
        for ii = 1:ncells
            N{ii}.CellRadius = str2num(ParValue); % (in distance units) If 
            % param.conntype is 'spatial' the user must set the radius of
            % mitral and granule cells
        end
        
        % New parameters must be added here...
        
    otherwise
        disp(['parameter ' ParName ' does not exist']);
        
end
end

function Cell = SetCellPosition(Cell,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the position on different neurons when no preset is
% presented.
% The main function of this program is neurogenesismain.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncells = length(Cell);

auxposy = param.BulbH / floor(sqrt(ncells));

if sqrt(ncells) == floor(sqrt(ncells))
    auxposx = param.BulbW / sqrt(ncells);
else
    auxposx = param.BulbW / (sqrt(ncells) + 1);
end

countcell = 0;
x = 0;
y = 0;

while countcell < ncells
    countcell = countcell + 1;
    if y >= param.BulbH
        y = 0;
        x = x + auxposx;
    end
    Cell{countcell}.X = mean([x,(x + auxposx)]);
    Cell{countcell}.Y = mean([y,(y + auxposy)]);
    y = y + auxposy;
end
end

function Cell = PreSetCellPosition(Cell,param,source)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the position on different neurons based on presets from
% excel files.
% Licurgo de Almeida 12/22/2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(param.ConnType(1:7),'spatial')
    
    if length(Cell) == param.nMitral
        for ii = 1:param.nMitral
            Cell{ii}.X = source.X(ii);
            Cell{ii}.matX = source.X(ii); %position in the pcolor matrix
            Cell{ii}.Y = source.Y(ii);
            Cell{ii}.matY = source.Y(ii);
            if param.flagtoroid == true
                Cell{ii}.Acorrec = 1;
            else
                Cell{ii}.Acorrec = AreaCorrection(source.Y(ii),param.BulbW,param.MitralRadius + param.GranuleRadius);
            end
        end
    else
        dSize = round(sqrt(param.mFactor));
        x = (1 / dSize:1 / dSize:param.BulbH);
        y = (1 / dSize:1 / dSize:param.BulbW);
        countx = 1;
        county = 1;
        for ii = 1:param.nGranule
            Cell{ii}.X = x(countx);
            Cell{ii}.matX = countx;
            Cell{ii}.Y = y(county);
            Cell{ii}.matY = county;
            if param.flagtoroid == true
                Cell{ii}.Acorrec = 1;
            else
                Cell{ii}.Acorrec = AreaCorrection(y(county),param.BulbW,param.MitralRadius + param.GranuleRadius);
            end
            countx = countx + 1;
            if countx > length(x)
                countx = 1;
                county = county + 1;
            end
        end
    end
    
elseif strcmp(param.ConnType,'circular')
    
    for ii = 1:length(source.X)
        Cell{ii}.X = source.X(ii);
    end
    
end
end

function AC = AreaCorrection(pos,Wsize,radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function defines the area-input correction factor. When the shape of
% the bulb is not a toroid, we have to multiply the inputs received by 
% neurons on the borders so that the overall activity to each neurons will
% be similar.
% This function is called from PreSetCellPosition.m
% Licurgo de Almeida
% 01/26/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (pos + radius) <= Wsize && (pos - radius) >= 0
    AC = 1;
else
    x = min([pos,Wsize - pos]);
    PA = ((pi / 2) * (radius^2)) + (x * sqrt((radius^2) - (x^2))) +...
    ((radius^2) * atan(x / sqrt((radius^2) - (x^2))));
    AC = (pi * (radius^2)) / PA;

end
end

function [OSN Pglo Glo Mitral Granule param] = NeuroActivity(OSN,Pglo,Glo,Mitral,Granule,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates the neuronal activity
%
% Licurgo de Almeida
% 11/03/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create initial setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gPosition = CreateGraPosition(Granule,param);

if param.flagRespiration == true
    Respiration = CreateRespFreq(param);
else
    Respiration = ones(1,round(param.tsim / param.dt));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set OSN parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% OSN external input
Iextosn = SetExtInput(param,OSN);


% Stores the voltage of each OSN at a given time
Vosn = zeros(param.nMitral,round(param.tsim / param.dt));
% Binary matrix recording only the spikes
Sosn = Vosn;

if param.SpikingInput == false
    Pfireosn = Vosn; %Pfireosn is the probability of firing when we use
    % non-spiking cells
end

refracosn = 2; % Refractory period after firing (ms)

Restosn = zeros(param.nMitral,1); % resting potential
Threshosn = Restosn; % current firing threshold
Hyperosn = Restosn; % hyperpolarization potential
Rosn = Restosn; % membrane resistance
tauosn = Restosn; % tau neuron
countrefracosn = Restosn; % Refractory period counter. This is assumed to
% be the same for all OSN so we don't need to add it to the for below.
Noiseosn = Restosn; % The parameter noise gives que number of changes in
% neuronal threshold

% ...New parameters should be added here and inside the 'for'...

for ii = 1:param.nMitral
    Restosn(ii) = OSN{ii}.Vrest;
    Hyperosn(ii) = OSN{ii}.Vhyper;
    Rosn(ii) = OSN{ii}.R;
    tauosn(ii) = OSN{ii}.tau;
    Threshosn(ii) = OSN{ii}.FThresh;
    Noiseosn(ii) = OSN{ii}.Noise;
end

% Initialize OSN potentials
Vosn(:,1) = Restosn;


if param.flagOSNdata == true %read pre-processed data for OSN
    OSNfile = strcat(param.outputPath,'OSNdata');
    load(OSNfile,'auxVosn');
    if size(auxVosn,1) == size(Vosn,1) && size(auxVosn,2) == size(Vosn,2)
        Vosn = auxVosn;
        load(OSNfile,'auxSosn');
        Sosn = auxSosn;
    else
        param.flagOSNdata = false; %if the sizes are not correct we don't change
        % the variables
        disp('param.flagOSNdata was changed to false');
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of OSN parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set periglomerular cell parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MOSNConpgl = eye(param.nMitral); % matrix with OSN-Pglo connection

Wosnpgl = eye(param.nMitral,param.nMitral); % Excitatory synaptic weights
% between OSNs and Pglo cells

if param.flagweights == true
    Wosnpgl = Wosnpgl * 0.5; % If we have learning all synapses start with
    % an average weight
end

% Stores the voltage of each Pglo at a given time
Vpgl = zeros(param.nMitral,round(param.tsim / param.dt));
% Binary matrix recording only the spikes
Spgl = Vpgl;

if param.SpikingPglo == false
    Pfirepgl = Vpgl; %Pfirepgl is the probability of firing when we use
    % non-spiking cells
end

refracpgl = 2; % Refractory period after firing (ms)

Restpgl = zeros(param.nMitral,1); % resting potential
Threshpgl = Restpgl; % current firing threshold
Hyperpgl = Restpgl; % hyperpolarization potential
Rpgl = Restpgl; % membrane resistance
taupgl = Restpgl; % tau neuron
countrefracpgl = Restpgl; % Refractory period counter. This is assumed to
% be the same for all PG cells so we don't need to add it to the for below.
Noisepgl = Restpgl; % The parameter noise gives que number of changes in
% neuronal threshold
gmaxAMPApgl = Restpgl; % Max AMPA conductance
tauAMPA1pgl = Restpgl; % AMPA's rising tau.
tauAMPA2pgl = Restpgl; % AMPA's falling tau.
EAMPApgl = Restpgl; % AMPA's Nernst potential.
IAChpgl = Restpgl; % Additional current when ACh is ON in non-spiking neurons.

% ...New parameters should be added here and inside the 'for'...

for ii = 1:param.nMitral
    Restpgl(ii) = Pglo{ii}.Vrest;
    Hyperpgl(ii) = Pglo{ii}.Vhyper;
    Rpgl(ii) = Pglo{ii}.R;
    taupgl(ii) = Pglo{ii}.tau;
    Threshpgl(ii) = Pglo{ii}.FThresh;
    Noisepgl(ii) = Pglo{ii}.Noise;
    gmaxAMPApgl(ii) = Pglo{ii}.gmaxAMPA;
    tauAMPA1pgl(ii) = Pglo{ii}.tauAMPA1;
    tauAMPA2pgl(ii) = Pglo{ii}.tauAMPA2;
    EAMPApgl(ii) = Pglo{ii}.EAMPA;
    IAChpgl(ii) = Pglo{ii}.IACh;
end

% Initialize OSN potentials
Vpgl(:,1) = Restpgl;

% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tAMPA0pgl = zeros(param.nMitral,1) - 10000000;
Iosnpgl = zeros(param.nMitral,1); % Input coming from OSN
maxgAMPApgl = getmaxg(param.dt,tauAMPA1pgl,tauAMPA2pgl); % Get max conductance
% amplitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of periglomerular parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set glomerulus cell parameters and variables. These cells are actually a
% compartment of mitral cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MOSNConglo = eye(param.nMitral); % matrix with OSN-Glo connection

Wosnglo = eye(param.nMitral,param.nMitral); % Excitatory synaptic weights
% between OSNs and Glo cells

if param.flagweights == true
    Wosnglo = Wosnglo * 0.5; % If we have learning all synapses start with
    % an average weight
end

% Stores the voltage of each Glo at a given time
Vglo = zeros(param.nMitral,round(param.tsim / param.dt));

Pfireglo = Vglo; %Pfireglo is the probability of firing when we use
% non-spiking cells


Restglo = zeros(param.nMitral,1); % resting potential
Threshglo = Restglo; % current firing threshold
Rglo = Restglo; % membrane resistance
tauglo = Restglo; % tau neuron
Noiseglo = Restglo; % The parameter noise gives que number of changes in
% neuronal threshold
gmaxAMPAglo = Restglo; % Max AMPA conductance
tauAMPA1glo = Restglo; % AMPA's rising tau.
tauAMPA2glo = Restglo; % AMPA's falling tau.
EAMPAglo = Restglo; % AMPA's Nernst potential.
gmaxGABAglo = Restglo; % Max GABA conductance
tauGABA1glo = Restglo; % GABA's rising tau.
tauGABA2glo = Restglo; % GABA's falling tau.
EGABAglo = Restglo; % GABA's Nernst potential.
IAChglo = Restglo; % Additional current when ACh is ON in non-spiking neurons.

% ...New parameters should be added here and inside the 'for'...

for ii = 1:param.nMitral
    Restglo(ii) = Glo{ii}.Vrest;
    Rglo(ii) = Glo{ii}.R;
    tauglo(ii) = Glo{ii}.tau;
    Threshglo(ii) = Glo{ii}.FThresh;
    Noiseglo(ii) = Glo{ii}.Noise;
    gmaxAMPAglo(ii) = Glo{ii}.gmaxAMPA;
    tauAMPA1glo(ii) = Glo{ii}.tauAMPA1;
    tauAMPA2glo(ii) = Glo{ii}.tauAMPA2;
    EAMPAglo(ii) = Glo{ii}.EAMPA;
    gmaxGABAglo(ii) = Glo{ii}.gmaxGABA;
    tauGABA1glo(ii) = Glo{ii}.tauGABA1;
    tauGABA2glo(ii) = Glo{ii}.tauGABA2;
    EGABAglo(ii) = Glo{ii}.EGABA;
    IAChglo(ii) = Glo{ii}.IACh;
end

% Initialize potentials
Vglo(:,1) = Restglo;

% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tAMPA0glo = zeros(param.nMitral,1) - 10000000;
Iosnglo = zeros(param.nMitral,1); % Input coming from OSN
maxgAMPAglo = getmaxg(param.dt,tauAMPA1glo,tauAMPA2glo); % Get max conductance
% amplitude

% GABA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tGABA0glo = zeros(param.nMitral,1) - 10000000;
Ipglglo = zeros(param.nMitral,1); % Input coming from Pglo
maxgGABAglo = getmaxg(param.dt,tauGABA1glo,tauGABA2glo); % Get max conductance
% amplitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of glomerulus parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Mitral cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Stores the voltage of each Mitral cell at a given time
Vmit = zeros(param.nMitral,round(param.tsim / param.dt));
% Binary matrix recording only the spikes
Smit = Vmit;
AHPmit = Vmit;

refracmit = 2; % Refractory period after firing (ms)


Restmit = zeros(param.nMitral,1); % resting potential
Threshmit = Restmit; % current firing threshold
Hypermit = Restmit; % hyperpolarization potential
Rmit = Restmit; % membrane resistance
taumit = Restmit; % tau neuron
countrefracmit = Restmit; % Refractory period counter. This is assumed to
% be the same for all Mitral cells so we don't need to add it to the for below.
Noisemit = Restmit; % The parameter noise gives que number of changes in
% neuronal threshold
gmaxAMPAmit = Restmit; % Max AMPA conductance
tauAMPA1mit = Restmit; % AMPA's rising tau.
tauAMPA2mit = Restmit; % AMPA's falling tau.
EAMPAmit = Restmit; % AMPA's Nernst potential.
gmaxGABAmit = Restmit; % Max GABA conductance
tauGABA1mit = Restmit; % GABA's rising tau.
tauGABA2mit = Restmit; % GABA's falling tau.
EGABAmit = Restmit; % GABA's Nernst potential.
Acorrmit = Restmit; % Area correction of the inputs.
if param.flagAHP == true
    gmaxAHPmit = Restmit; % Max AHP conductance
    tauAHP1mit = Restmit; % AHP's rising tau.
    tauAHP2mit = Restmit; % AHP's falling tau.
    EAHPmit = Restmit; % AHP's Nernst potential.
end

Wglomit = eye(param.nMitral,param.nMitral); % Excitatory synaptic weights
% between Glos and mitral cells
if param.flagweights == true
    Wglomit = Wglomit * 0.5; % If we have learning all synapses start with
    % an average weight
end

Wpglmit = Wglomit; %for now, these matrices are the same

if param.flagLoadConn  == true
    Savefile = strcat(param.outputPath,'ConnData');
    load(Savefile,'MGraConmit');
    Savefile = strcat(param.outputPath,'LearnData'); %Mit-gra learned weights
    load(Savefile,'Wgramit');
else
    MGraConmit = zeros(param.nMitral,param.nGranule); % matrix with gra-mit
    % connection
    Wgramit = zeros(param.nMitral,param.nGranule); % Excitatory synaptic weights
    % between mitral cells and granule cells
end

% ...New parameters should be added here and inside the 'for'...

for ii = 1:param.nMitral
    Restmit(ii) = Mitral{ii}.Vrest;
    Hypermit(ii) = Mitral{ii}.Vhyper;
    Rmit(ii) = Mitral{ii}.R;
    taumit(ii) = Mitral{ii}.tau;
    gmaxAMPAmit(ii) = Mitral{ii}.gmaxAMPA;
    tauAMPA1mit(ii) = Mitral{ii}.tauAMPA1;
    tauAMPA2mit(ii) = Mitral{ii}.tauAMPA2;
    EAMPAmit(ii) = Mitral{ii}.EAMPA;
    gmaxGABAmit(ii) = Mitral{ii}.gmaxGABA;
    tauGABA1mit(ii) = Mitral{ii}.tauGABA1;
    tauGABA2mit(ii) = Mitral{ii}.tauGABA2;
    EGABAmit(ii) = Mitral{ii}.EGABA;
    if param.flagAHP == true
        gmaxAHPmit(ii) = Mitral{ii}.gmaxAHP;
        tauAHP1mit(ii) = Mitral{ii}.tauAHP1;
        tauAHP2mit(ii) = Mitral{ii}.tauAHP2;
        EAHPmit(ii) = Mitral{ii}.EAHP;
    end
    
    
    if strcmp(param.ConnType(1:7),'spatial')
        Acorrmit(ii) = Mitral{ii}.Acorrec;
    else
        Acorrmit(ii) = 1;
    end
    
    % The firing threshold is different for each neuron. After the neuron
    % fires, this threshold changes
    Threshmit(ii) = Mitral{ii}.FThresh;
    Noisemit(ii) = Mitral{ii}.Noise;
    
    if param.flagLoadConn  == false
        [MGraConmit,Wgramit] = SetConnection(MGraConmit,Wgramit,ii,Mitral{ii},param,gPosition,'M');
    end
    
    Mitral{ii}.Connections = MGraConmit(ii,:);
    Mitral{ii}.ConWeights = Wgramit(ii,:);
end

% Initialize Mitral cells potentials
Vmit(:,1) = Restmit;

% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tAMPA0mit = zeros(param.nMitral,1) - 10000000;

Iglomit = zeros(param.nMitral,1); % Input coming from Glo
maxgAMPAmit = getmaxg(param.dt,tauAMPA1mit,tauAMPA2mit); % Get max conductance
% amplitude

% GABA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. Each Mitral cell can be
% connected to a varied number of granule cells
tGABA0mit = zeros(param.nGranule,1) - 10000000;
Igramit = zeros(param.nMitral,1); % Input coming from Granule cells
maxgGABAmit = getmaxg(param.dt,tauGABA1mit,tauGABA2mit); % Get max conductance
% amplitude


% AHP time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0.
if param.flagAHP == true
    tAHP0mit = zeros(param.nMitral,1) - 10000000;
    maxgAHPmit = getmaxg(param.dt,tauAHP1mit,tauAHP2mit); % Get max conductance
end
Iahpmit = zeros(param.nMitral,1); % AHP current
% amplitude


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Mitral cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Granule cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stores the voltage of each Granule cell at a given time
Vgra = zeros(param.nGranule,round(param.tsim / param.dt));
% Binary matrix recording only the spikes
Sgra = Vgra;


refracgra = 2; % Refractory period after firing (ms)

Restgra = zeros(param.nGranule,1); % resting potential
Threshgra = Restgra; % current firing threshold
Hypergra = Restgra; % hyperpolarization potential
Rgra = Restgra; % membrane resistance
taugra = Restgra; % tau neuron
countrefracgra = Restgra; % Refractory period counter. This is assumed to
% be the same for all Granule cells so we don't need to add it to the for below.
Noisegra = Restgra; % The parameter noise gives que number of changes in
% neuronal threshold
gmaxAMPAgra = Restgra; % Max AMPA conductance
tauAMPA1gra = Restgra; % AMPA's rising tau.
tauAMPA2gra = Restgra; % AMPA's falling tau.
EAMPAgra = Restgra; % AMPA's Nernst potential.
gmaxGABAgra = Restgra; % Max GABA conductance
tauGABA1gra = Restgra; % GABA's rising tau.
tauGABA2gra = Restgra; % GABA's falling tau.
EGABAgra = Restgra; % GABA's Nernst potential.
Acorrgra = Restgra; % Area correction of the inputs.

if param.flagAHP == true
    gmaxAHPgra = Restgra; % Max AHP conductance
    tauAHP1gra = Restgra; % AHP's rising tau.
    tauAHP2gra = Restgra; % AHP's falling tau.
    EAHPgra = Restgra; % AHP's Nernst potential.
end


MMitCongra = MGraConmit'; % matrix with mit-gra connection
if param.flagLoadConn  == true
    Savefile = strcat(param.outputPath,'ConnData');
    load(Savefile,'MGraCongra');
    Savefile = strcat(param.outputPath,'LearnData'); %Mit-gra learned weights
    load(Savefile,'Wgragra');
else
    MGraCongra = zeros(param.nGranule); % matrix with gra-gra connection
    Wgragra = zeros(param.nGranule);
end

MahpCongra = eye(param.nGranule); % matrix with AHP connection
Wahpgra = MahpCongra; % AHP weights


% ...New parameters should be added here and inside the 'for'...

for ii = 1:param.nGranule
    Restgra(ii) = Granule{ii}.Vrest;
    Hypergra(ii) = Granule{ii}.Vhyper;
    Rgra(ii) = Granule{ii}.R;
    taugra(ii) = Granule{ii}.tau;
    gmaxAMPAgra(ii) = Granule{ii}.gmaxAMPA;
    tauAMPA1gra(ii) = Granule{ii}.tauAMPA1;
    tauAMPA2gra(ii) = Granule{ii}.tauAMPA2;
    EAMPAgra(ii) = Granule{ii}.EAMPA;
    gmaxGABAgra(ii) = Granule{ii}.gmaxGABA;
    tauGABA1gra(ii) = Granule{ii}.tauGABA1;
    tauGABA2gra(ii) = Granule{ii}.tauGABA2;
    EGABAgra(ii) = Granule{ii}.EGABA;
    if param.flagAHP == true
        gmaxAHPgra(ii) = Granule{ii}.gmaxAHP;
        tauAHP1gra(ii) = Granule{ii}.tauAHP1;
        tauAHP2gra(ii) = Granule{ii}.tauAHP2;
        EAHPgra(ii) = Granule{ii}.EAHP;
    end
    
    if strcmp(param.ConnType(1:7),'spatial')
        Acorrgra(ii) = Granule{ii}.Acorrec;
    else
        Acorrgra(ii) = 1;
    end
    
    % The firing threshold is different for each neuron. After the neuron
    % fires, this threshold changes
    Threshgra(ii) = Granule{ii}.FThresh;
    Noisegra(ii) = Granule{ii}.Noise;
    
    if param.GraGracon == true && param.flagLoadConn  == false
        [MGraCongra Wgragra]= SetConnection(MGraCongra,Wgragra,ii,Granule{ii},param,gPosition,'G');
    end
    Granule{ii}.Connections = MMitCongra(ii,:);
end

% Initialize Granule cells potentials
Vgra(:,1) = Restgra;

% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tAMPA0gra = zeros(param.nMitral,1) - 10000000;
Imitgra = zeros(param.nGranule,1); % Input coming from mitral cells
Wmitgra = ones(param.nGranule,param.nMitral); % Synaptic weights
maxgAMPAgra = getmaxg(param.dt,tauAMPA1gra,tauAMPA2gra); % Get max conductance
% amplitude

% GABA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. Each Granule cell can be
% connected to a varied number of other granule cells
tGABA0gra = zeros(param.nGranule,1) - 10000000;
Igragra = zeros(param.nGranule,1); % Input coming from other granule cells
maxgGABAgra = getmaxg(param.dt,tauGABA1gra,tauGABA2gra); % Get max conductance
% amplitude

% AHP time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0.
if param.flagAHP == true
    tAHP0gra = zeros(param.nGranule,1) - 10000000;
    maxgAHPgra = getmaxg(param.dt,tauAHP1gra,tauAHP2gra); % Get max conductance
end
Iahpgra = zeros(param.nGranule,1); % AHP current
% amplitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Granule cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin neuron simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start loop
for tt = 2:round(param.tsim / param.dt)
    t = tt * param.dt; % current time
    
    % OSN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if param.flagOSNdata == true
        I = Vosn(:,tt - 1) == param.SpikeV;
        
        % Mitral cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tAMPA0glo(I) = t; % new t0 for the Glomerulus AMPA current
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Periglomerular cell cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if param.SpikingPglo == true
            tAMPA0pgl(I) = t; % new t0 for the Pglo AMPA current
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        Vosn(:,tt) = Vosn(:,tt - 1) + (param.dt ./ tauosn(:)) .* (Rosn(:) .*...
            Respiration(tt) .* Iextosn(:,tt) - Vosn(:,tt - 1) + Restosn(:));
        
        if param.SpikingInput == true
            % If the neuron fired last cycle, neuron potential hyperpotentializes
            I = Vosn(:,tt - 1) == param.SpikeV;
            Vosn(I,tt) = Hyperosn(I);
            
            % Mitral cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tAMPA0mit(I) = t; % new t0 for the Mitral AMPA current
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Periglomerular cell cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tAMPA0pgl(I) = t; % new t0 for the Pglo AMPA current
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            countrefracosn(I) = refracosn / param.dt; % neurons that just fired
            % get into the refractory period
            I = countrefracosn > 0;
            countrefracosn(I) = countrefracosn(I) - 1;
            I = find(countrefracosn == 0); % if countrefracosn = 0 the neuron can fire
            % again
            if ~isempty(I)
                if param.flagnoiseosn == true
                    aux_J = SpikeNoise(Restosn(I),Threshosn(I),param,...
                        Vosn(I,tt),Noiseosn(I),'OSN');
                    J = find(aux_J);
                else
                    J = find(Vosn(I,tt) >= Threshosn(I));
                end
                if ~isempty(J)
                    Vosn(I(J),tt) = param.SpikeV; % Action potential
                    Sosn(I(J),tt) = 1; % Record spike time
                end
            end
        else
            Pfireosn(:,tt) = NSpikeP(Vosn(:,tt),Restosn,Threshosn,Noiseosn);
        end
    end
    
    % Periglomerular Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if param.SpikingInput == true
        % Get OSN input to Pglo cells
        Iosnpgl(:) = SetI(gmaxAMPApgl(1),tauAMPA1pgl(1),tauAMPA2pgl(1),t,tAMPA0pgl,maxgAMPApgl(1),...
            Wosnpgl,EAMPApgl(1),Vpgl(:,tt - 1),MOSNConpgl);
    else
        Iosnpgl(:) = SetInoSpike(gmaxAMPApgl,Pfireosn(:,tt),EAMPApgl,Vpgl(:,tt - 1));
    end
    
    % Periglomerular cell potential
    Vpgl(:,tt) = Vpgl(:,tt - 1) + (param.dt ./ taupgl(:)) .* (Rpgl(:) .*...
        Iosnpgl(:) + IAChpgl(:) - Vpgl(:,tt - 1) + Restpgl(:));
    
    if param.SpikingPglo == true
        % If the neuron fired last cycle, neuron potential hyperpotentializes
        I = Vpgl(:,tt - 1) == param.SpikeV;
        Vpgl(I,tt) = Hyperpgl(I);
        
        % Mitral cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tGABA0glo(I) = t; % new t0 for the Glo cell GABA current
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        countrefracpgl(I) = (refracpgl / param.dt);
        I = countrefracpgl > 0;
        countrefracpgl(I) = countrefracpgl(I) - 1;
        
        
        I = find(countrefracpgl == 0); % if countrefracmit = 0 the neuron can fire
        % again
        
        
        if ~isempty(I)
            if param.flagnoisepgl == true
                aux_J = SpikeNoise(Restpgl(I),Threshpgl(I),param,Vpgl(I,tt),Noisepgl(I),'Pglo');
                J = find(aux_J);
            else
                J = find(Vpgl(I,tt) >= Threshpgl(I));
            end
            if ~isempty(J)
                Vpgl(I(J),tt) = param.SpikeV; % Action potential
                Spgl(I(J),tt) = 1; % Record spike time
            end
        end
    else
        Pfirepgl(:,tt) = NSpikeP(Vpgl(:,tt),Restpgl,Threshpgl,Noisepgl);
    end
    
    % Glomeruli Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get OSN input to Glo cells
    if param.SpikingInput == true
        Iosnglo(:) = SetI(gmaxAMPAglo(1),tauAMPA1glo(1),tauAMPA2glo(1),t,tAMPA0glo,maxgAMPAglo(1),...
            Wosnglo,EAMPAglo(1),Vglo(:,tt - 1),MOSNConglo);
    else
        Iosnglo(:) = SetInoSpike(gmaxAMPAglo,Pfireosn(:,tt),EAMPAglo,Vglo(:,tt - 1));
    end
    
    % Get Periglomerular input to Glo cells
    if param.SpikingPglo == true
        Ipglglo(:) = SetI(gmaxGABAglo(1),tauGABA1glo(1),tauGABA2glo(1),t,tGABA0glo,maxgGABAglo(1),...
            Wpglglo,EGABAglo(1),Vglo(:,tt - 1),MOSNConglo);
    else
        Ipglglo(:) = SetInoSpike(gmaxGABAglo,Pfirepgl(:,tt),EGABAglo,Vglo(:,tt - 1));
    end
    
    % Glo cell potential
    Vglo(:,tt) = Vglo(:,tt - 1) + (param.dt ./ tauglo(:)) .* (Rglo(:) .*...
        (Iosnglo(:) + Ipglglo(:)) - Vglo(:,tt - 1) + Restglo(:));

    Pfireglo(:,tt) = NSpikeP(Vglo(:,tt),Restglo,Threshglo,Noiseglo);

    % Mitral Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get Glo input to Mitral cells. This current is applied directly in
    % Vmit, no channel conductance.
    Iglomit(:) = Pfireglo(:,tt) .* gmaxAMPAmit;
    
    % Get Granule input to Mitral cells
    Igramit(:) = SetI(gmaxGABAmit(1),tauGABA1mit(1),tauGABA2mit(1),t,tGABA0mit,maxgGABAmit(1),...
        Wgramit,EGABAmit(1),Vmit(:,tt - 1),MGraConmit);
    
    % Get AHP to Mitral cells
    if param.flagAHP == true
        Iahpmit(:) = SetI(gmaxAHPmit(1),tauAHP1mit(1),tauAHP2mit(1),t,tAHP0mit,maxgAHPmit(1),...
            Wpglmit,EAHPmit(1),Vmit(:,tt - 1),MOSNConglo);
    end
    
    % Mitral cell potential
    Vmit(:,tt) = Vmit(:,tt - 1) + (param.dt ./ taumit(:)) .* (Rmit(:) .*...
        (Iahpmit(:) + Iglomit(:) + Igramit(:) .* ...
        Acorrmit(:)) - Vmit(:,tt - 1) + Restmit(:));
    
    AHPmit(:,tt) = Iahpmit(:);
    
    % If the neuron fired last cycle, neuron potential hyperpotentializes
    I = Vmit(:,tt - 1) == param.SpikeV;
    Vmit(I,tt) = Hypermit(I);
    

    % Granule cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tAMPA0gra(I) = t; % new t0 for the Granule cell AMPA current
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if param.flagAHP == true
        tAHP0mit(I) = t; % new t0 for the Mitral cell AHP
    end
    
    countrefracmit(I) = refracmit / param.dt; % neurons that just fired get into the
    % refractory period
    I = countrefracmit > 0;
    countrefracmit(I) = countrefracmit(I) - 1;
    
    
    I = find(countrefracmit == 0); % if countrefracmit = 0 the neuron can fire
    % again
    
    if ~isempty(I)
        if param.flagnoisemit == true
            aux_J = SpikeNoise(Restmit(I),Threshmit(I),param,Vmit(I,tt),...
                Noisemit(I),'Mitral');
            J = find(aux_J);
        else
            J = find(Vmit(I,tt) >= Threshmit(I));
        end
        if ~isempty(J)
            Vmit(I(J),tt) = param.SpikeV; % Action potential
            Smit(I(J),tt) = 1; % Record spike time
        end
    end
    
    % Granule Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get Mitral input to Granule cells
    Imitgra(:) = SetI(gmaxAMPAgra(1),tauAMPA1gra(1),tauAMPA2gra(1),t,tAMPA0gra,maxgAMPAgra(1),...
        Wmitgra,EAMPAgra(1),Vgra(:,tt - 1),MMitCongra);
    
    % Get Granule input to Granule cells
    if param.GraGracon == true
        Igragra(:) = SetI(gmaxGABAgra(1),tauGABA1gra(1),tauGABA2gra(1),t,tGABA0gra,maxgGABAgra(1),...
            Wgragra,EGABAgra(1),Vgra(:,tt - 1),MGraCongra);
    end
    % Get AHP to Granule cells
    if param.flagAHP == true
        Iahpgra(:) = SetI(gmaxAHPgra(1),tauAHP1gra(1),tauAHP2gra(1),t,tAHP0gra,maxgAHPgra(1),...
            Wahpgra,EAHPgra(1),Vgra(:,tt - 1),MahpCongra);
    end
    
    % Granule cell potential
    Vgra(:,tt) = Vgra(:,tt - 1) + (param.dt ./ taugra(:)) .* (Rgra(:) .*...
        (Iahpgra(:) + Igragra(:) + Imitgra(:) .* Acorrgra(:))...
        - Vgra(:,tt - 1) + Restgra(:));
    
    % If the neuron fired last cycle, neuron potential hyperpotentializes
    I = Vgra(:,tt - 1) == param.SpikeV;
    Vgra(I,tt) = Hypergra(I);
    
    
    % Mitral cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tGABA0mit(I) = t; % new t0 for the Mitral cell GABA current
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tGABA0gra(I) = t; % new t0 for the Granule cell GABA current
    
    if param.flagAHP == true
        tAHP0gra(I) = t; % new t0 for the Granule cell AHP
    end
    
    countrefracgra(I) = refracgra / param.dt; % neurons that just fired get
    % into the refractory period
    I = countrefracgra > 0;
    countrefracgra(I) = countrefracgra(I) - 1;
    
    I = find(countrefracgra == 0); % if countrefracgra = 0 the neuron can fire
    % again
    
    if ~isempty(I)
        if param.flagnoisegra == true
            aux_J = SpikeNoise(Restgra(I),Threshgra(I),param,Vgra(I,tt),...
                Noisegra(I),'Granule');
            J = find(aux_J);
        else
            J = find(Vgra(I,tt) >= Threshgra(I));
        end
        if ~isempty(J)
            Vgra(I(J),tt) = param.SpikeV; % Action potential
            Sgra(I(J),tt) = 1; % Record spike time
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OSN and Mitral cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:param.nMitral
    OSN{ii}.V = Vosn(ii,:); % Save neuronal activity
    if param.SpikingInput == true
        OSN{ii}.S = Sosn(ii,:); % Save spike time
    else
        OSN{ii}.S = Pfireosn(ii,:);
    end
    Mitral{ii}.V = Vmit(ii,:); % Save neuronal activity
    Mitral{ii}.S = Smit(ii,:); % Save spike time
    Pglo{ii}.V = Vpgl(ii,:); % Save neuronal activity
    if param.SpikingPglo == true
        Pglo{ii}.S = Spgl(ii,:); % Save spike time
    else
        Pglo{ii}.S = Pfirepgl(ii,:);
    end
    Glo{ii}.V = Vglo(ii,:); % Save neuronal activity
    Glo{ii}.S = Pfireglo(ii,:);
end

% Granule cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:param.nGranule
    Granule{ii}.V = Vgra(ii,:); % Save neuronal activity
    Granule{ii}.S = Sgra(ii,:); % Save spike time
end

% Save OSN file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param.flagOSNdata == false && param.SpikingInput == true
    OSNfile = strcat(param.outputPath,'OSNdata');
    auxVosn = Vosn;
    auxSosn = Sosn;
    save(OSNfile,'auxVosn','auxSosn');
end

% Save connections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param.flagLoadConn  == false
    Savefile = strcat(param.outputPath,'ConnData');
    save(Savefile,'MGraConmit','MGraCongra');
    Savefile = strcat(param.outputPath,'LearnData'); %Mit-Gra learned weights
end
save(Savefile,'Wgramit','Wgragra');
end

function gPosition = CreateGraPosition(Cell,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create position vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(param.ConnType(1:7),'spatial')
    gPosition = zeros(param.nGranule,3);
    for ii = 1:param.nGranule
        gPosition(ii,1) = Cell{ii}.X;
        gPosition(ii,2) = Cell{ii}.Y;
        if strcmp(param.ConnType,'spatialfile');
            gPosition(ii,3) = Cell{ii}.ConnChanceMitGra;
        else
            gPosition(ii,3) = Cell{ii}.CellRadius;
        end
    end
elseif strcmp(param.ConnType,'circular')
    gPosition = zeros(param.nGranule,2);
    for ii = 1:param.nGranule
        gPosition(ii,1) = Cell{ii}.X;
        gPosition(ii,2) = Cell{ii}.CellRadius;
    end
else
    gPosition = []; % If ConnType is not spatial, the position is not necessary
end
end

function Iext = SetExtInput(param,Neuron)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the inputs to OSN neurons
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 11/02/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Iext = zeros(param.nMitral,round(param.tsim / param.dt));

for ii = 1:param.nMitral
    Iext(ii,round(param.tinit / param.dt) : round(param.tfinal / param.dt)) = Neuron{ii}.input;
end

Iext = Iext * param.Iext;

end

function R = CreateRespFreq(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates an artificial respiration to modulate the bulb
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 04/06/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = param.dt:param.dt:param.tsim;
time = time / 1000; %converting the time to seconds
R = -cos(2 * pi * param.RespFreq * time);
R = R + 1;
R = R / 2;
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

function Ic = SetInoSpike(gmax,P,E,V)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns the nonspike pglo input for a given excitation
%
% Licurgo de Almeida
% 08/04/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = P .* gmax;
    Ic = g .* (E - V);
end

function [Mat Weight] = SetConnection(Mat,Weight,celln,cell,param,gPosition,conntypecell)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates a matrix of connections between mitral-granule
% cells or granule-granule cells.
% This function is called from NeuroActivity.m
%
% Licurgo de Almeida
% 11/17/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The connection configuration changes depending on the method:
% randperm = spiacial configuration isn't taken in consideration, cells
% just connect randomly to each other.
% spatialvar = connection are determined based on the position of each
% cell, the farther the cells are, the higher the chances.
% spatialfix =  same as spatialvar, but the chances of connection are
% fixed.
% spatialfile = neurons are spatially organized, but the chance of
% connection is preset and given by the granule cell layer output matrix.
% In this example, gposition(:,3) stores the connection chance and not the
% radius.
% circular = neurons are organized around a ring (1D).

switch lower(param.ConnType)
    case 'randperm'
        if strcmp(conntypecell,'G')
            Mat = randpermG(Mat,celln,cell.NGraCon,param.nGranule);
            Weight = ones(size(Weight));
        elseif strcmp(conntypecell,'M')
            Mat = randpermM(Mat,celln,cell.NGraCon,param.nGranule);
            Weight = ones(size(Weight));
        end
        
    case 'spatialvar'
        [Mat Weight] = spatialconn(Mat,Weight,celln,cell,param,gPosition,conntypecell);
        
    case 'spatialfix'
        [Mat Weight] = spatialconn(Mat,Weight,celln,cell,param,gPosition,conntypecell);
        
    case 'spatialfile'
        [Mat Weight] = fileconn(Mat,Weight,celln,param,gPosition,conntypecell);
        
    case 'circular'
        [Mat Weight] = circularconn(Mat,Weight,celln,cell,param,gPosition,conntypecell);
        
    otherwise
        disp(['Unknown method: ' connmethod]);
end


end

function Mat = randpermG(Mat,celln,Nconn,Ngra)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is only used for connetions type G using method randperm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Ngra == 1
    disp('The network has only one granule cell. No G-G connetions was set');
else
    conns = zeros(Nconn,1);
    countconn = 1;
    while countconn <= Nconn
        conns(countconn) = ceil(rand * Ngra);
        if conns(countconn) ~= celln
            countconn = countconn + 1;
        end
    end
    for ii = 1:length(conns)
        Mat(celln,conns(ii)) = Mat(celln,conns(ii)) + 1;
    end
end

end

function Mat = randpermM(Mat,celln,Nconn,Ngra)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is only used for connetions type M using method randperm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Nconn <= Ngra % in this case, we can use randperm
    conns = randperm(Ngra);
    conns = conns(1:Nconn);
    Mat(celln,conns) = Mat(celln,conns) + 1;
else % if the number of connections is larger than the number of granule
    % cells, some cells are going to be connected twice
    conns = ceil(rand(1,Nconn) * Ngra);
    for ii = 1:Nconn
        Mat(celln,conns(ii)) = Mat(celln,conns(ii)) + 1;
    end
end

end

function [Mat Weight] = spatialconn(Mat,Weight,celln,cell,param,gPosition,conntypecell)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used for spatial connections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode = 'farther';

distcells = zeros(size(gPosition,2),1); % store the distance betwenn two
% cells
x = cell.X;
y = cell.Y;
if param.flagtoroid == true % here we use the complete toroid
    
    for ii = 1:size(gPosition,1)
        distcells(ii) = sqrt((x - gPosition(ii,1))^2 + (y - gPosition(ii,2))^2);
        if (x + cell.CellRadius > param.BulbH) && (y + cell.CellRadius > param.BulbW) % > >
            auxdist = sqrt(((x - param.BulbH) - gPosition(ii,1))^2 + ((y - param.BulbW) - gPosition(ii,2))^2);
            if auxdist < distcells(ii)
                distcells(ii) = auxdist;
            end
        end
        if (x + cell.CellRadius > param.BulbH) && (y - cell.CellRadius < 0) % > <
            auxdist = sqrt(((x - param.BulbH) - gPosition(ii,1))^2 + ((y + param.BulbW) - gPosition(ii,2))^2);
            if auxdist < distcells(ii)
                distcells(ii) = auxdist;
            end
        end
        if (x - cell.CellRadius < 0) && (y + cell.CellRadius > param.BulbW) % < >
            auxdist = sqrt(((x + param.BulbH) - gPosition(ii,1))^2 + ((y - param.BulbW) - gPosition(ii,2))^2);
            if auxdist < distcells(ii)
                distcells(ii) = auxdist;
            end
        end
        if (x - cell.CellRadius < 0) && (y - cell.CellRadius < 0) % < <
            auxdist = sqrt(((x + param.BulbH) - gPosition(ii,1))^2 + ((y + param.BulbW) - gPosition(ii,2))^2);
            if auxdist < distcells(ii)
                distcells(ii) = auxdist;
            end
        end
        if x + cell.CellRadius > param.BulbH
            auxdist = sqrt(((x - param.BulbH) - gPosition(ii,1))^2 + (y - gPosition(ii,2))^2);
            if auxdist < distcells(ii)
                distcells(ii) = auxdist;
            end
        end
        if x - cell.CellRadius < 0
            auxdist = sqrt(((x + param.BulbH) - gPosition(ii,1))^2 + (y - gPosition(ii,2))^2);
            if auxdist < distcells(ii)
                distcells(ii) = auxdist;
            end
        end
        if y + cell.CellRadius > param.BulbW
            auxdist = sqrt((x - gPosition(ii,1))^2 + ((y - param.BulbW) - gPosition(ii,2))^2);
            if auxdist < distcells(ii)
                distcells(ii) = auxdist;
            end
        end
        
        if y - cell.CellRadius < 0
            auxdist = sqrt((x - gPosition(ii,1))^2 + ((y + param.BulbW) - gPosition(ii,2))^2);
            if auxdist < distcells(ii)
                distcells(ii) = auxdist;
            end
        end
        
    end
else % if we're using a cilinder shape, only the high circular connections are important
    for ii = 1:size(gPosition,1)
        distcells(ii) = sqrt((x - gPosition(ii,1))^2 + (y - gPosition(ii,2))^2);
        if x + cell.CellRadius > param.BulbH
            auxdist = sqrt(((x - param.BulbH) - gPosition(ii,1))^2 + (y - gPosition(ii,2))^2);
            if auxdist < distcells(ii)
                distcells(ii) = auxdist;
            end
        end
        if x - cell.CellRadius < 0
            auxdist = sqrt(((x + param.BulbH) - gPosition(ii,1))^2 + (y - gPosition(ii,2))^2);
            if auxdist < distcells(ii)
                distcells(ii) = auxdist;
            end
        end
    end
end
    
switch conntypecell
    case 'M'
        ConnChance = param.ConnChanceMitGra;
    case 'G'
        ConnChance = param.ConnChanceGraGra;
end
if strcmpi(param.ConnType,'spatialfix')
    I = find(distcells(:) < cell.CellRadius + gPosition(:,3));
    if param.flagweights == false
        J = rand(length(I),1) > ConnChance;
        I(J) = [];
        Weight(celln,I) = 1;
    else
        for ii = 1:length(I)
            Weight(celln,I(ii)) = rand * ConnChance * 2;
        end
    end
elseif strcmpi(param.ConnType,'spatialvar')
    if param.flagweights == false
        if strcmpi(mode,'closer')
            % The closer the higher chance of connection
            chance = (1 - (distcells - gPosition(:,3)) ./ cell.CellRadius) * (ConnChance * 2);
        elseif strcmpi(mode,'farther')
            % The farther, the higher chance of connection
            chance = ((distcells - gPosition(:,3)) ./ cell.CellRadius);
        end
        I = chance > 1;
        chance(I) = 0;
        chance = chance * (ConnChance * 2);
        I = find(rand(length(distcells),1) < chance);
    else
        I = find(distcells(:) < cell.CellRadius + gPosition(:,3));
        if strcmpi(mode,'closer')
            for ii = 1:length(I)
                Weight(celln,I(ii)) = rand * (1 - (distcells(I(ii)) - gPosition(I(ii),3)) ./ cell.CellRadius) * (ConnChance * 2);
            end
        elseif strcmpi(mode,'farther')
            for ii = 1:length(I)
                Weight(celln,I(ii)) = rand * ((distcells(I(ii)) - gPosition(I(ii),3)) ./ cell.CellRadius) * (ConnChance * 2);
            end
        end
    end
        
end

Mat(celln,I) = 1;

if strcmp(conntypecell,'G') % In gran-gran connections, cells are not
    % connected to each other
    Mat(celln,celln) = 0;
end


end

function [Mat Weight] = fileconn(Mat,Weight,celln,param,gPosition,conntypecell)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is for preset connections base on the granule layer file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch conntypecell
    case 'M'
        ConnChance = param.ConnChanceMitGra;
    case 'G'
        ConnChance = param.ConnChanceGraGra;
end

chance = gPosition(:,3)' * ConnChance;
I = rand(1,length(chance)) < chance;
Mat(celln,I) = 1;
Weight(celln,I) = 1;
end

function [Mat Weight] = circularconn(Mat,Weight,celln,cell,param,gPosition,conntypecell)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is for preset circular 1D networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_UseCellRadius = false; % if true, the distance between cells and the radius
% of their dendrited are taken in consideration; if false, the cells are
% don't have any spacial limitation

switch conntypecell
    case 'M'
        ConnChance = param.ConnChanceMitGra;
    case 'G'
        ConnChance = param.ConnChanceGraGra;
end

if flag_UseCellRadius == true
    x = cell.X;
    distcells = abs(x - gPosition(:,1)); % store the distance betwenn two
    % cells
    I = find(distcells > (param.nMitral / 2));
    distcells(I) = abs(distcells(I) - param.nMitral);
    totaldist = gPosition(:,2) + cell.CellRadius;
    I = find(distcells <= totaldist);
    J = randperm(length(I));
    I = (I(J(1:ceil(length(I) * ConnChance))));
    if param.flagweights == false
        Weight(celln,I) = 1;
    else
        for ii = 1:length(I)
            Weight(celln,I(ii)) = rand * ConnChance * 2;
        end
    end
else
    I = randperm(param.nGranule);
    nconn = round(param.nGranule * ConnChance);
    if nconn > 0
        I = I(1:nconn);
        if strcmp(conntypecell,'G')
            J = find(I == celln);
            if ~isempty(J)
                I(J) = [];
            else
                I(1) = [];
            end
            Mat(celln,celln) = 1;
        end
    end
end

Mat(celln,I) = 1;

end

function I = SpikeNoise(Rest,Thres,param,V,noise,tnet)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the spiking chance
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 12/20/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(tnet,'OSN')
    Vlimit = 0e-3; % Limits where firing chances varied between 0 and 1
    pfactor = 3;
    chance = 1;
elseif strcmp(tnet,'Pglo')
    Vlimit = 0e-3; % Limits where firing chances varied between 0 and 1
    pfactor = 1;
    chance = 1;
elseif strcmp(tnet,'Granule')
    Vlimit = 0e-3; % Limits where firing chances varied between 0 and 1
    pfactor = 3;
    chance = 1;
elseif strcmp(tnet,'Mitral')
    Vlimit = 0.1e-3; % Limits where firing chances varied between 0 and 1
    pfactor = 2;
    chance = 1;
end

Rest = Rest - Vlimit;
%Thres = Thres + Vlimit;
y = (V - Rest + noise) ./ (Thres -  Rest + noise);
J = y <= 0;
y(J) = 0;
J = y > 1;
y(J) = 1;
y = y.^pfactor;
I = rand(length(Rest),1) <= y * (param.dt * chance);
    
end

function P = NSpikeP(V,R,T,noise)
P = (V - R + noise) ./ (T - R + noise);
J = P <= 0;
P(J) = 0;
J = P > 1;
P(J) = 1;

end

function Ic = SetI(gmax,tau1,tau2,t,t0,normg,W,E,V,Mcon)

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

g = gmax * (((tau1 * tau2) / (tau1 - tau2)) *...
    (exp(-(t - t0) / tau1) - exp(-(t - t0) / tau2))) / normg;

Ic = (Mcon * g) .* (E - V);

% Ic = zeros(length(tau1),1);

% if param.flagweights == false
%     %Mcon = Mcon';
%     for ii = 1:length(tau1)
%         I = Mcon(ii,:) > 0;
%         
%         g = gmax(ii) * ((tau1(ii) * tau2(ii)) / (tau1(ii) - tau2(ii))) *...
%             (exp(-(t - t0(I)) / tau1(ii)) - exp(-(t - t0(I)) / tau2(ii))) / normg(ii);
% 
%         
%         auxIc = Mcon(ii,I)' .* g .* (E(ii) - V(ii));
%         Ic(ii) = sum(auxIc);
%     end
%     
% else
%     
%     for ii = 1:length(tau1)
%         I = Mcon(ii,:) > 0;
%         
%         g = gmax(ii) * ((tau1(ii) * tau2(ii)) / (tau1(ii) - tau2(ii))) *...
%             (exp(-(t - t0(I)) / tau1(ii)) - exp(-(t - t0(I)) / tau2(ii))) / normg(ii);
%         
%         auxIc = Mcon(ii,I)' .* W(ii,I)' .* g .* (E(ii) - V(ii));
%         
%         Ic(ii) = sum(auxIc);
%     end
% end
end

function Matrices = CompactGraMat(param,Grasource,Granule)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function reduces the size of the granular activity matrix.
%
% Licurgo de Almeida
% 01/07/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mfactor = sqrt(param.mFactor); %this number is fixed in CreateCells.m and
% necessarily a integer

x = max(Grasource.X);
y = max(Grasource.Y);
Gramato = zeros(x,y); %Original granular matrix
for ii = 1:length(Granule)
    Gramato(Grasource.X(ii),Grasource.Y(ii)) = sum(Granule{ii}.S);
end
Matrices.OriginalGraMat = Gramato;

[A B] = size(Grasource.outputs);
mat = zeros(A,B);

for ii = 0:A - 1
    for jj = 0:B - 1
        C = mean(mean(Gramato((ii * mfactor) + 1:1:(ii * mfactor) + mfactor,(jj * mfactor) + 1:1:(jj * mfactor) + mfactor)));
        mat(ii + 1,jj + 1) = C;
    end
end

Matrices.Gramat = mat;
end

function SaveMitralFile(Mitral,param)

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

MitralMat = zeros(length(Mitral),length(Mitral{1}.S));

for ii = 1:length(Mitral)
    MitralMat(ii,:) = Mitral{ii}.S;
end

cd ..

sfile = strcat('PiriformModel/',param.outputPath,param.mitfilename);
save(sfile,'MitralMat');

cd BulbModel
end