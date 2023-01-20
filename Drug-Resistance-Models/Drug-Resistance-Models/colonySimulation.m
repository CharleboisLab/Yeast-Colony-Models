function [phenotypeAves,genotypeAves,appearanceTimes,establishmentTimes,fixationTimes,cellGenotypeNos,cellNumbers,cellCountTable,time] = colonySimulation(diffusionRates,endConditions,mutationSettings,latticeSize,lattices,DRUGS_TYPE,displaySettings)
% COLONYSIMULATION models the growth of a colony and the development of 
% genetic and non-genetic drug resistance under the influence of both
% nutrient concentrations and drugs.
%
%   INPUT
%
%       diffusionRates is a structure array containing the following fields:
%            -- nutrientDiff is the number of steps that a nutrient packet takes during its 
%               random walk, correlates to the diffusion coefficient of the media.
%            -- drugDiff is the number of steps that drug packet takes during its
%               random walk in SNG mode, correlates to the diffusion coefficient
%               of the drugs through the media.
%
%       endConditions is a structure array containing the following fields:
%            -- 'finalCellNo' is the number  of cells that colony should reach in 
%                order fo the program to stop.
%            -- 'finalTime' is the number of timesteps until the program stops.
%            -- 'endCondition' is a str which determines whether the end condition will
%                be based on timesteps of final cell count: 
%               'timestep' = timestep, 'cell count' = cell count
%
%       mutationSettings is a structure array containing the following fields:
%            -- 'mutationProb' is a float between 0.0 and 1.0 giving the probability
%                of a cell mutating during division
%            -- 'mutationType' is a str which sets the mutation model:
%               'SNG' = cells can be susceptible, non-genetically drug
%                resistant, or genetically drug resistant
%               'megaPlate' = cells are assigned a genotype which sets the level
%                of resistance and a phenotype number which is added to the
%                genotype number to either increase or decrease the resistance
%            -- 'switchUpProb' is a float between 0.0 and 1.0 giving the probability
%                of a cell's phenotype increasing in resistance: in SNG mode it is the
%                probability of a cell becoming non-genetically resistant, 
%                in megaPlate mode it is the probability of a cell's phenotype 
%                number increasing
%            -- 'switchDownProb' is a float between 0.0 and 1.0 giving the probability
%                of a cell's phenotype decreasing in resistance: in SNG mode it is the
%                probability of a cell losing non-genetic resistance, 
%                in megaPlate mode it is the probability of a cell's phenotype
%                number decreasing
%            -- 'inheritPhenotype' is an int. 0 if the daughter cell does not
%                inherit the phenotype of the mother, 1 if the daughter cell does
%                inherit the phenotype of the mother.
%            -- 'NDivisionProb' is a float between 0.0 and 1.0 giving the
%                probability of a non-genetically resistant cell dividing in the
%                presence of a drug in the 'SNG' mode.
%
%        latticeSize is a structure array containing the following fields:
%            -- 'x' is the x dimension size
%            -- 'y' is the y dimension size
%
%        lattices is a structure array containing the following fields:
%            -- 'cells' is a matrix of structure arrays with the following field:
%                      + 'cell' is 0 if no cell is present, a Cell object if
%                         a cell is present
%            -- 'nutrients' is a matrix containing the number of nutrient packets
%                at each point on the lattice
%            -- 'states' is a matrix of ones and zeros, where a 1 states that a
%                cell is a present and 0 states that no cell is present
%            -- 'drugs' is a matrix containing the number of drug packets at each
%                point on the lattice
%
%       DRUGS_TYPE is a str setting what type of drug is being applied:
%                  'static' = static drug, 'cidal' = cidal drug, 'off' = no drug
%
%        displaySettings is a structure array with the following fields:
%            -- 'displayType' is a str which determines how the colony is visualised:
%                 'show' = display the colony updating as the simulation runs,
%                 'save' = save an image of the colony with each update,
%                 'none' = do not visualise the colony.
%            -- 'folderName' is a str giving the folder the images will be saved
%                to if 'displayType' == 'save'
%            -- 'map' is a three-column matrix of values between 0.0 and 1.0
%                creating a custom colormap for visualising the colony and
%                different cell genotypes/phenotypes and drug concentrations
%
%   OUTPUT
%
%       phenotypeAves is an array with the average phenotype+genotype values of all
%           the cells at various timesteps in the first row and the timesteps
%           when averages were sampled in the second row for 'megaPlate' mode
%           Returns an array of zeroes for 'SNG' mode.
%
%       genotypeAves is an array with the average genotype values of all the cells
%           at various timesteps in the first row and the timesteps when averages
%           were sampled in the second row for 'megaPlate' mode. 
%           Returns an array of zeroes for 'SNG' mode.
%
%       appearanceTimes is an array of the times at which each genotype first
%           appeared for 'megaPlate' mode, a single value for the time of first
%           appearance of the mutant for 'SNG' mode.
%
%       establishmentTimes is an array of the times at which each genotype
%           established (makes up 5% of the population) for 'megaPlate' mode, a
%           single value for the establishment time of the mutant for 'SNG' mode.
%
%       fixationTimes is an array of the times at which each genotype fixated
%           (makes up 95% of the population) for 'megaPlate' mode, a single value for
%           the fixation time of the mutant for 'SNG' mode.
%
%       cellGenotypeNos is an array of the final number of cells of each genotype
%           for 'megaPlate' mode. Returns array of zeros for 'SNG'.
%       
%       cellNumbers is a structrue array containing the total number of cells
%           and for 'SNG' mode, the number of cells of each SNG type.
%           Contains the following fields:
%               -- 'cellsno' is the total number of living cells
%               -- 'susceptibles' is the number of susceptible cells
%               -- 'nongenet' is the number of non-genetically drug
%                   resistant cells
%               -- 'genet' is the number of genetically drug resistant cells
%
%       cellCountTable is a table containing the number of each type of
%           cell in SNG mode at various times through the simulation. The
%           first column contains the timesteps when the cell counts were
%           sampled, the second column contains susceptible cell counts,
%           the third column contains non-genetically resistant cell
%           counts, the fourth column contains genetically resistant cell
%           counts. In 'megaPlate' mode just returns 0.
%
%       time is the number of timesteps the simulation ran for
%
%   See also CELL

%% Initialisation

% Real time equivalent to timestep in minutes
stepConv = 0.0129;

% Image magnification
magnification = 4;

% Get end condition info
endCondition = endConditions.('endCondition');
finalCellNo = endConditions.('finalCellNo');
finalTime = endConditions.('finalTime');

% Get diffusion info
nutrientDiff = diffusionRates.('nutrientDiff');
if strcmp(mutationSettings.('mutationType'),'SNG')
    drugDiff = diffusionRates.('drugDiff');
else
    drugDiff = 0;
end

% Get display info
displayType = displaySettings.('displayType');
map = displaySettings.('map');
folderName = displaySettings.('folderName');

% Sets the range of cells which will be considered under the rules as the
% cells immediate surrounding the start cell.
[iMin,iMax,jMin,jMax] = boundaryCoors(lattices.('states'),0,0,0,0,0,0,false);
[iMinD,iMaxD,jMinD,jMaxD] = boundaryCoors(lattices.('drugs'),0,0,0,0,0,0,false);

colouredStateLattice = ones(latticeSize.('y'),latticeSize.('x'));
if strcmp(mutationSettings.('mutationType'),'megaPlate')
    colouredStateLattice = assignDrugColours(colouredStateLattice,lattices,DRUGS_TYPE,iMinD,iMaxD,jMinD,jMaxD);
end

[phenotypeImage,cellPhenotypeNos] = createColouredLattice(map,colouredStateLattice,lattices,DRUGS_TYPE,mutationSettings.('mutationType'),true,iMinD,iMaxD,jMinD,jMaxD,iMin,iMax,jMin,jMax,magnification);
[genotypeImage,cellGenotypeNos] = createColouredLattice(map,colouredStateLattice,lattices,DRUGS_TYPE,mutationSettings.('mutationType'),false,iMinD,iMaxD,jMinD,jMaxD,iMin,iMax,jMin,jMax,magnification);

% Display conditions

if strcmp(displayType,'show')
    % Visualise initial states.

    disp('Press any key...')
    
    if strcmp(mutationSettings.('mutationType'),'megaPlate')
        imshowpair(phenotypeImage,genotypeImage,'montage')
    elseif strcmp(mutationSettings.('mutationType'),'SNG')
        imshow(genotypeImage)
    end
    tt = 2;

    pause
elseif strcmp(displayType,'save')
    % Save initial states.
    if strcmp(mutationSettings.('mutationType'),'megaPlate')
        frame = imfuse(phenotypeImage,genotypeImage,'montage');
        imFileName = strcat(folderName,'/frame1.jpg');
        imwrite(frame,imFileName,'jpg')
    elseif strcmp(mutationSettings.('mutationType'),'SNG')
        mkdir(folderName)
        imFileName = strcat(folderName,'/frame1.jpeg');
        imwrite(genotypeImage,imFileName,'jpg')
    end
    tt = 2;
end

time = 0; 

%% Run
cellNumbers = countCells(lattices,iMin,jMin,iMax,jMax);   % number of cells alive on the lattice
prevtime = 0;
appearanceTimeFound = false;
establishmentTimeFound = false;
fixationTimeFound = false;
genotypeAves = zeros(2,finalTime);
phenotypeAves = zeros(2,finalTime);

if strcmp(mutationSettings.('mutationType'),'megaPlate')
    appearanceTimes = zeros(1,15);
    establishmentTimes = zeros(1,15);
    fixationTimes = zeros(1,15);
    totalGenotype = 0;
    totalPhenotype = 0;
    for k = 1:length(cellGenotypeNos)
        totalPhenotype = totalPhenotype + k*cellPhenotypeNos(k);
        totalGenotype = totalGenotype + k*cellGenotypeNos(k);
    end
    aveGenotype = totalGenotype/cellNumbers.('cellsno');
    avePhenotype = totalPhenotype/cellNumbers.('cellsno');
    genotypeAves(1,1) = aveGenotype;
    genotypeAves(2,1) = time;
    phenotypeAves(1,1) = avePhenotype;
    phenotypeAves(2,1) = time;
elseif strcmp(mutationSettings.('mutationType'),'SNG')
    cellCountS = zeros(finalTime,1);
    cellCountN = zeros(finalTime,1);
    cellCountG = zeros(finalTime,1);
    timestep = zeros(finalTime,1);
    cellCountS(1) = cellNumbers.('susceptibles');
    cellCountN(1) = cellNumbers.('nongenet');
    cellCountG(1) = cellNumbers.('genet');
    timestep(1) = time;
    appearanceTimes = 0;
    establishmentTimes = 0;
    fixationTimes = 0;
end

t = 2;

whileCondition = setEndCondition(finalCellNo,finalTime,endCondition,cellNumbers.('cellsno'),time);

% Run program until the end condition is met
while whileCondition  
    % Apply nutrient movement rule
    containsNutrient = false;
    while ~containsNutrient
        i_r = randi([1 (latticeSize.('y'))],1);
        j_r = randi([1 (latticeSize.('x'))],1);
            if lattices.('nutrients')(i_r,j_r) > 0
                containsNutrient = true;
                [i_n,j_n] = randomWalk(i_r,j_r,nutrientDiff,latticeSize.('x'),latticeSize.('y'));
                lattices.('nutrients')(i_n,j_n) = lattices.('nutrients')(i_n,j_n) + 1;
                lattices.('nutrients')(i_r,j_r) = lattices.('nutrients')(i_r,j_r) - 1;
            end
    end
   
    % Apply drug movement rule
    if ~strcmp(DRUGS_TYPE,'off')
        if drugDiff > 0
            if iMinD ~= -1
                containsDrugs = false;
                while ~containsDrugs
                i_r = randi([iMinD iMaxD],1);
                j_r = randi([jMinD jMaxD],1);
                    if lattices.('drugs')(i_r,j_r) > 0
                        containsDrugs = true;
                        [i_n,j_n] = randomWalk(i_r,j_r,drugDiff,latticeSize.('x'),latticeSize.('y'));
                        lattices.('drugs')(i_n,j_n) = lattices.('drugs')(i_n,j_n) + 1;
                        lattices.('drugs')(i_r,j_r) = lattices.('drugs')(i_r,j_r) - 1;
                        if lattices.('states')(i_n,j_n) == 1
                            [death, lattices,cellNumbers] = drugResponse(lattices,i_n,j_n,DRUGS_TYPE,cellNumbers);
                            if death
                                stateLattice = lattices.('states');
                                [iMin,iMax,jMin,jMax] = boundaryCoors(stateLattice,iMin,iMax,jMin,jMax,latticeSize.('x'),latticeSize.('y'),true);
                            end
                        end
                        if (i_r > iMaxD) || (i_r < iMinD) || (j_r > jMaxD) || (j_r < jMinD)
                            [iMinD,iMaxD,jMinD,jMaxD] = boundaryCoors(lattices.('drugs'),iMinD,iMaxD,jMinD,jMaxD,latticeSize.('x'),latticeSize.('y'),true);
                        end
                    end
                end
            end
        end
    end

    aliveCellChosen = false;
    while ~aliveCellChosen
        % Randomly find a lattice site containing a cell, then apply
        % division rules.
        i_r = randi([iMin iMax],1);
        j_r = randi([jMin jMax],1);

        if lattices.('states')(i_r,j_r) == 1
            [lattices,cellNumbers] = updateCell(lattices,i_r,j_r,DRUGS_TYPE,cellNumbers);
            if emptyNeighbour(lattices.('states'),[i_r,j_r])
                if ~isa(lattices.('cells')(i_r,j_r).('cell'),'double')
                    % Otherwise, attempt to bud
                    [divided,lattices,cellNumbers] = lattices.('cells')(i_r,j_r).('cell').divide(lattices,[i_r,j_r],cellNumbers);
                    if divided
                        stateLattice = lattices.('states');
                        [iMin,iMax,jMin,jMax] = boundaryCoors(stateLattice,iMin,iMax,jMin,jMax,latticeSize.('x'),latticeSize.('y'),true);
                    end
                end
            end
            aliveCellChosen = true;
        end
    end
    
    time = time + 1;
    whileCondition = setEndCondition(finalCellNo,finalTime,endCondition,cellNumbers.('cellsno'),time);
    
    % Check for first appearances, establishment, and fixation
    
    if strcmp(mutationSettings.('mutationType'),'SNG')
        if ~appearanceTimeFound
            if cellNumbers.('genet') > 0
                if strcmp(displayType,'show')
                    disp('first appearance time:')
                    disp('timesteps:')
                    disp(time)
                    disp('time in minutes:')
                    disp(time * stepConv)
                end
                appearanceTimeFound = true;
                appearanceTimes = time;
            end
        end
        if ~establishmentTimeFound
            if (cellNumbers.('genet')/cellNumbers.('cellsno')) >= 0.05
                if strcmp(displayType,'show')
                    disp('establishment time:')
                    disp('timesteps:')
                    disp(time)
                    disp('time in minutes:')
                    disp(time * stepConv)
                end
                establishmentTimeFound = true;
                establishmentTimes = time;
            end
        end
        if ~fixationTimeFound
            if (cellNumbers.('genet')/cellNumbers.('cellsno')) >= 0.95
                if strcmp(displayType,'show')
                    disp('fixation time:')
                    disp('timesteps:')
                    disp(time)
                    disp('time in minutes:')
                    disp(time * stepConv)
                end
                fixationTimeFound = true;
                fixationTimes = time;
            end
        end
    elseif strcmp(mutationSettings.('mutationType'),'megaPlate')
        for k = 2:15
            if appearanceTimes(k) == 0
                if cellGenotypeNos(k) > 0
                    string = strcat('first appearance time G',int2str(k),':');
                    if strcmp(displayType,'show')
                        disp(string)
                        disp('timesteps:')
                        disp(time)
                        disp('time in minutes:')
                        disp(time * stepConv)
                    end
                    appearanceTimes(k) = time;
                end
            end
            if establishmentTimes(k) == 0
                if (cellGenotypeNos(k)/cellNumbers.('cellsno')) >= 0.05
                    string = strcat('establishment time G',int2str(k),':');
                    if strcmp(displayType,'show')
                        disp(string)
                        disp('timesteps:')
                        disp(time)
                        disp('time in minutes:')
                        disp(time * stepConv)
                    end
                    establishmentTimes(k) = time;
                end
            end
            if fixationTimes(k) == 0
                if (cellGenotypeNos(k)/cellNumbers.('cellsno')) >= 0.95
                    string = strcat('fixation time G',int2str(k),':');
                    if strcmp(displayType,'show')
                        disp(string)
                        disp('timesteps:')
                        disp(time)
                        disp('time in minutes:')
                        disp(time * stepConv)
                    end
                    fixationTimes(k) = time;
                end
            end
        end
    end

    % Display or save lattice.
    if stepConv * (time - prevtime) > 1
        
        [phenotypeImage,cellPhenotypeNos] = createColouredLattice(map,colouredStateLattice,lattices,DRUGS_TYPE,mutationSettings.('mutationType'),true,iMinD,iMaxD,jMinD,jMaxD,iMin,iMax,jMin,jMax,magnification);
        [genotypeImage,cellGenotypeNos] = createColouredLattice(map,colouredStateLattice,lattices,DRUGS_TYPE,mutationSettings.('mutationType'),false,iMinD,iMaxD,jMinD,jMaxD,iMin,iMax,jMin,jMax,magnification);
        
        if strcmp(displayType,'show')

            if strcmp(mutationSettings.('mutationType'),'megaPlate')
                imshowpair(phenotypeImage,genotypeImage,'montage')
            elseif strcmp(mutationSettings.('mutationType'),'SNG')
                imshow(genotypeImage)
            end
            
            tt = tt + 1;
            
        elseif strcmp(displayType,'save')

            imFileName = strcat(folderName,'/frame',int2str(tt),'.jpg');
            if strcmp(mutationSettings.('mutationType'),'megaPlate')
                frame = imfuse(phenotypeImage,genotypeImage,'montage');
                imwrite(frame,imFileName,'jpg')
            elseif strcmp(mutationSettings.('mutationType'),'SNG')
                imwrite(genotypeImage,imFileName,'jpg')
            end
            
            tt = tt + 1;
        end
        
        if strcmp(mutationSettings.('mutationType'),'megaPlate')
            totalGenotype = 0;
            totalPhenotype = 0;
            for k = 1:length(cellGenotypeNos)
                totalPhenotype = totalPhenotype + k*cellPhenotypeNos(k);
                totalGenotype = totalGenotype + k*cellGenotypeNos(k);
            end
            aveGenotype = totalGenotype/cellNumbers.('cellsno');
            avePhenotype = totalPhenotype/cellNumbers.('cellsno');
            genotypeAves(1,t) = aveGenotype;
            genotypeAves(2,t) = time;
            phenotypeAves(1,t) = avePhenotype;
            phenotypeAves(2,t) = time;
        elseif strcmp(mutationSettings.('mutationType'),'SNG')
            cellCountS(t) = cellNumbers.('susceptibles');
            cellCountN(t) = cellNumbers.('nongenet');
            cellCountG(t) = cellNumbers.('genet');
            timestep(t) = time;
        end
        t = t + 1;
        
        prevtime = time;     
    end
end 

% Display or save final lattice.

[phenotypeImage,cellPhenotypeNos] = createColouredLattice(map,colouredStateLattice,lattices,DRUGS_TYPE,mutationSettings.('mutationType'),true,iMinD,iMaxD,jMinD,jMaxD,iMin,iMax,jMin,jMax,magnification);
[genotypeImage,cellGenotypeNos] = createColouredLattice(map,colouredStateLattice,lattices,DRUGS_TYPE,mutationSettings.('mutationType'),false,iMinD,iMaxD,jMinD,jMaxD,iMin,iMax,jMin,jMax,magnification);

if strcmp(mutationSettings.('mutationType'),'SNG')
    cellCountS(t) = cellNumbers.('susceptibles');
    cellCountN(t) = cellNumbers.('nongenet');
    cellCountG(t) = cellNumbers.('genet');
    timestep(t) = time;
    cellCountS = cellCountS(1:t);
    cellCountN = cellCountN(1:t);
    cellCountG = cellCountG(1:t);
    timestep = timestep(1:t);
    cellCountTable = table(timestep,cellCountS,cellCountN,cellCountG);
elseif strcmp(mutationSettings.('mutationType'),'megaPlate')
    totalGenotype = 0;
    totalPhenotype = 0;
    for k = 1:length(cellGenotypeNos)
        totalPhenotype = totalPhenotype + k*cellPhenotypeNos(k);
        totalGenotype = totalGenotype + k*cellGenotypeNos(k);
    end
    aveGenotype = totalGenotype/cellNumbers.('cellsno');
    avePhenotype = totalPhenotype/cellNumbers.('cellsno');
    genotypeAves(1,t) = aveGenotype;
    genotypeAves(2,t) = time;
    phenotypeAves(1,t) = avePhenotype;
    phenotypeAves(2,t) = time;
    genotypeAves = genotypeAves(:,1:t);
    phenotypeAves = phenotypeAves(:,1:t);
    cellCountTable = 0;   % Just to give program something to output
end

if strcmp(displayType,'show')
    if strcmp(mutationSettings.('mutationType'),'megaPlate')
        imshowpair(phenotypeImage,genotypeImage,'montage')
    elseif strcmp(mutationSettings.('mutationType'),'SNG')
        imshow(genotypeImage)
    end 
    disp('All done!')
elseif strcmp(displayType,'save')
    imFileName = strcat(folderName,'/frame',int2str(tt),'.jpg');
    if strcmp(mutationSettings.('mutationType'),'megaPlate')
        frame = imfuse(phenotypeImage,genotypeImage,'montage');
        imwrite(frame,imFileName,'jpg')
    elseif strcmp(mutationSettings.('mutationType'),'SNG')
        imwrite(genotypeImage,imFileName,'jpg')
        
        writetable(cellCountTable,[folderName,'/cellCounts.csv'])
    end
end
end

%% Create Genotype Coloured Lattice

function [image,cellResistanceNos] = createColouredLattice(map,colouredStateLattice,lattices,drugType,mutationType,phenotype,iMinD,iMaxD,jMinD,jMaxD,iMin,iMax,jMin,jMax,magnification)
    % CREATECOLOUREDLATTICE Updates colouredStateLattice based on the cell
    % lattice and for 'megaPlate' mode gives a count of the number of cells of
    % each genotype or phenotype level up to 15.
    %
    %   If mutationType = 'SNG':
    %       1 = no cell or drug packets are present
    %       2 = drug packet(s) is present, no cell is present
    %       3 = susceptible cell is present
    %       4 = genetically resistant cell is present
    %       5 = non-genetically resistant cell is present
    %   If mutationType = 'megaPlate':
    %       States:
    %           1-7 = no cell is present, assigned at start by assignDrugColours
    %           8 = cell with resistance level 0 or 1
    %           9 = cell with resistance level 2
    %           10 = cell with resistance level 3
    %           11 = cell with resistance level 4
    %           12 = cell with resistance level 5
    %           13 = cell with resistance level 6
    %           14 = cell with resistance level 7+
    %       If phenotype = true, resistance = genotype level + phenotype level
    %       If phenotype = false, resistance = genotype level
    %
    %   INPUT
    %       map is a three-column matrix of 5 RGB triplets for 'SNG' and 14 RGB triplets for 'megaPlate'
    %       colouredStateLattice is the lattice of colour states to be updated by createColouredLattice
    %       lattices is a structure array containing the following fields:
    %            -- 'cells' is a matrix of structure arrays with the following field:
    %                      + 'cell' is 0 if no cell is present, a Cell object if
    %                         a cell is present
    %            -- 'nutrients' is a matrix containing the number of nutrient packets
    %                at each point on the lattice
    %            -- 'states' is a matrix of ones and zeros, where a 1 states that a
    %                cell is a present and 0 states that no cell is present
    %            -- 'drugs' is a matrix containing the number of drug packets at each
    %                point on the lattice
    %       drugType is a str setting what type of drug is being applied:
    %           'static' = static drug, 'cidal' = cidal drug, 'off' = no drug
    %       mutationType is a str setting the mutation model: 'SNG' or 'megaPlate'
    %       phenotype is a bool which is true if phenotype resistance levels
    %         are to be counted and displayed and false if genotype levels are to
    %         be counted and displayed for 'megaPlate' mode.
    %       iMinD is the smallest row number of the region within which the
    %         function will update states according to drug concentration in
    %         'SNG' mode
    %       iMaxD is the largest row number of the region within which the
    %         function will update states according to drug concentration in
    %         'SNG' mode
    %       jMinD is the smallest column number of the region within which the
    %         function will update states according to drug concentration in
    %         'SNG' mode
    %       jMaxD is the largest column number of the region within which the
    %         function will update states according to drug concentration in
    %         'SNG' mode
    %       iMin is the smallest row number of the region within which the
    %         function will update states according to cell presence and resistance
    %       iMax is the largest row number of the region within which the
    %         function will update states according to cell presence and resistance
    %       jMin is the smallest column number of the region within which the
    %         function will update states according to drug concentration
    %       jMax is the largest column number of the region within which the
    %         function will update states according to drug concentration
    %       magnification is the magnification factor for imresize
    %
    %   OUTPUT
    %       image is an RGB image representing the current state of the lattice
    %       cellResistanceNos is a list containing the number of cells with
    %         resistance level 0 or 1 at index 1 and resistance level i at
    %         index i for i > 1, where resistance level is determined by the
    %         genotype level with the phenotype modifier if 'phenotype' is true
    %         and determined solely by genotype level if 'phenotype' is false.
    %         Applies to 'megaPlate' mode.

    cellResistanceNos = zeros(1,15);
    if strcmp(mutationType,'megaPlate')
        % assign states based on resistnce levels
        for i = iMin:iMax
            for j = jMin:jMax
                if (~isa(lattices.('cells')(i,j).('cell'),'double')) && lattices.('states')(i,j) == 1
                    if phenotype
                        resistance = lattices.('cells')(i,j).('cell').genotype + lattices.('cells')(i,j).('cell').phenotype;
                        % update cell number count
                        if resistance < 1
                            cellResistanceNos(1) = cellResistanceNos(1) + 1;
                        else
                            cellResistanceNos(resistance) = cellResistanceNos(resistance) + 1;
                        end
                    else
                        resistance = lattices.('cells')(i,j).('cell').genotype;
                        % update cell number count
                        cellResistanceNos(resistance) = cellResistanceNos(resistance) + 1;
                    end
                    if resistance == 0
                        colouredStateLattice(i,j) = 8;
                    elseif resistance < 7
                        colouredStateLattice(i,j) = resistance + 7;
                    else
                        colouredStateLattice(i,j) = 14;
                    end
                else
                    % return lattice point state to backgroun drug colour if no cell is present anymore
                    colouredStateLattice = assignDrugColours(colouredStateLattice,lattices,drugType,i,i,j,j);
                end
            end
        end
    elseif strcmp(mutationType,'SNG')
        % update drug colours
        if ~strcmp(drugType,'off')
            if iMinD ~= -1
                for i = iMinD:iMaxD
                    for j = jMinD:jMaxD
                        if lattices.('drugs')(i,j) > 0
                            colouredStateLattice(i,j) = 2;
                        else
                            colouredStateLattice(i,j) = 1;
                        end
                    end
                end
            end
        end
        % update cell colours
        for i = iMin:iMax
            for j = jMin:jMax
                if (~isa(lattices.('cells')(i,j).('cell'),'double')) && lattices.('states')(i,j) == 1
                    if lattices.('cells')(i,j).('cell').mutated
                        colouredStateLattice(i,j) = 4;
                    elseif lattices.('cells')(i,j).('cell').NGresistant
                        colouredStateLattice(i,j) = 5;
                    else
                        colouredStateLattice(i,j) = 3;
                    end
                end
            end
        end
    end
    % create cell image
    [newLat,newmap] = imresize(colouredStateLattice,map,magnification,'Colormap','optimized');
    image = ind2rgb(newLat,newmap);
end

function [colouredStateLattice] = assignDrugColours(colouredStateLattice,lattices,drugType,iMinD,iMaxD,jMinD,jMaxD)
    % ASSIGNDRUGCOLOURS Assign states based on drug concentration for 'megaPlate' mode.
    %
    %   INPUT
    %       colouredStateLattice is the lattice of colour states to be updated by assignDrugColours
    %       lattices is a structure array containing the following fields:
    %            -- 'cells' is a matrix of structure arrays with the following field:
    %                      + 'cell' is 0 if no cell is present, a Cell object if
    %                         a cell is present
    %            -- 'nutrients' is a matrix containing the number of nutrient packets
    %                at each point on the lattice
    %            -- 'states' is a matrix of ones and zeros, where a 1 states that a
    %                cell is a present and 0 states that no cell is present
    %            -- 'drugs' is a matrix containing the number of drug packets at each
    %                point on the lattice
    %       drugType is a str setting what type of drug is being applied:
    %           'static' = static drug, 'cidal' = cidal drug, 'off' = no drug
    %       iMinD is the smallest row number of the region within which the
    %         function will update states according to drug concentration
    %       iMaxD is the largest row number of the region within which the
    %         function will update states according to drug concentration 
    %       jMinD is the smallest column number of the region within which the
    %         function will update states according to drug concentration
    %       jMaxD is the largest column number of the region within which the
    %         function will update states according to drug concentration
    %   
    %   OUTPUT
    %       colouredStateLattice is the updated colouredStateLattice
    
    if ~strcmp(drugType,'off')
        if iMinD ~= -1
            for i = iMinD:iMaxD
                for j = jMinD:jMaxD
                    if lattices.('drugs')(i,j) > 99999
                        colouredStateLattice(i,j) = 7;
                    elseif lattices.('drugs')(i,j) > 9999
                        colouredStateLattice(i,j) = 6;
                    elseif lattices.('drugs')(i,j) > 999
                        colouredStateLattice(i,j) = 5;
                    elseif lattices.('drugs')(i,j) > 99
                        colouredStateLattice(i,j) = 4;
                    elseif lattices.('drugs')(i,j) > 9
                        colouredStateLattice(i,j) = 3;
                    elseif lattices.('drugs')(i,j) > 0
                        colouredStateLattice(i,j) = 2;
                    else
                        colouredStateLattice(i,j) = 1;
                    end
                end
            end
        end
    end
end

%% Update Yeast Cell Function

function [lattices, cellNumbers] = updateCell(lattices,i,j,drugType,cellNumbers)
    % UPDATECELL Updates the cell at (i,j), including updating its divisibility, its phenotype, its nutrient consumption and drug
    % response, and updates cell counts accordingly.
    %
    %   INPUT
    %       lattices is a structure array containing the following fields:
    %            -- 'cells' is a matrix of structure arrays with the following field:
    %                      + 'cell' is 0 if no cell is present, a Cell object if
    %                         a cell is present
    %            -- 'nutrients' is a matrix containing the number of nutrient packets
    %                at each point on the lattice
    %            -- 'states' is a matrix of ones and zeros, where a 1 states that a
    %                cell is a present and 0 states that no cell is present
    %            -- 'drugs' is a matrix containing the number of drug packets at each
    %                point on the lattice
    %       i is the row number of the location of the cell to be updated
    %       j is the column number of the location of the cell to be updated
    %       drugType is a str setting the type of drug being applied:
    %           'static' = static drug, 'cidal' = cidal drug, 'off' = no drug
    %       cellNumbers is a structrue array containing the total number of cells
    %           and for 'SNG' mode, the number of cells of each SNG type.
    %           Contains the following fields:
    %               -- 'cellsno' is the total number of living cells
    %               -- 'susceptibles' is the number of susceptible cells
    %               -- 'nongenet' is the number of non-genetically drug
    %                   resistant cells
    %               -- 'genet' is the number of genetically drug resistant cells
    %   
    %   OUTPUT
    %       lattices is the updated version of the struct that was input
    %       cellNumbers is the updated version of the struc that was input
    
    if isa(lattices.('cells')(i,j).('cell'),'Cell')
        
        lattices.('states')(i,j) = 1;
        lattices.('cells')(i,j).('cell') = lattices.('cells')(i,j).('cell').updateDivisibility();
        
        oldState = lattices.('cells')(i,j).('cell').NGresistant;
        lattices.('cells')(i,j).('cell') = lattices.('cells')(i,j).('cell').updatePhenotype();
        newState = lattices.('cells')(i,j).('cell').NGresistant;
        if ~lattices.('cells')(i,j).('cell').mutated
            if (oldState == false) && (newState == true)
                cellNumbers.('susceptibles') = cellNumbers.('susceptibles') - 1;
                cellNumbers.('nongenet') = cellNumbers.('nongenet') + 1;
            elseif (oldState == true) && (newState == false)
                cellNumbers.('nongenet') = cellNumbers.('nongenet') - 1;
                cellNumbers.('susceptibles') = cellNumbers.('susceptibles') + 1;
            end
        end
        
        if lattices.('nutrients')(i,j) > 0
            lattices.('cells')(i,j).('cell') = lattices.('cells')(i,j).('cell').consumeNutrient();
            lattices.('nutrients')(i,j) = lattices.('nutrients')(i,j) - 1;
        end
        
        [~,lattices, cellNumbers] = drugResponse(lattices,i,j,drugType,cellNumbers);
    end
end

%% Respond to drug

function [death, lattices, cellNumbers] = drugResponse(lattices,i,j,drugType,cellNumbers)
    % DRUGRESPONSE Updates a cell in response to a drug and update cell counts accordingly
    %
    %   INPUT
    %       lattices is a structure array containing the following fields:
    %            -- 'cells' is a matrix of structure arrays with the following field:
    %                      + 'cell' is 0 if no cell is present, a Cell object if
    %                         a cell is present
    %            -- 'nutrients' is a matrix containing the number of nutrient packets
    %                at each point on the lattice
    %            -- 'states' is a matrix of ones and zeros, where a 1 states that a
    %                cell is a present and 0 states that no cell is present
    %            -- 'drugs' is a matrix containing the number of drug packets at each
    %                point on the lattice
    %       i is the row number of the location of the cell to be updated
    %       j is the column number of the location of the cell to be updated
    %       drugType is a str setting the type of drug being applied:
    %           'static' = static drug, 'cidal' = cidal drug, 'off' = no drug
    %       cellNumbers is a structrue array containing the total number of cells
    %           and for 'SNG' mode, the number of cells of each SNG type.
    %           Contains the following fields:
    %               -- 'cellsno' is the total number of living cells
    %               -- 'susceptibles' is the number of susceptible cells
    %               -- 'nongenet' is the number of non-genetically drug
    %                   resistant cells
    %               -- 'genet' is the number of genetically drug resistant cells
    %   
    %   OUTPUT
    %       death is a bool which is false if the cell did not die and true if it did
    %       lattices is the updated version of the struct that was input
    %       cellNumbers is the updated version of the struc that was input
    
    death = false;
    
    NGresistant = lattices.('cells')(i,j).('cell').NGresistant;
    mutated = lattices.('cells')(i,j).('cell').mutated;
    if ~ strcmp(drugType,'off')
        lattices.('cells')(i,j).('cell') = lattices.('cells')(i,j).('cell').respondToDrug(lattices.('drugs')(i,j),drugType);
    end
    
    if lattices.('cells')(i,j).('cell').state == 0
        death = true;
        lattices.('cells')(i,j).('cell') = 0;
        lattices.('states')(i,j) = 0;
        cellNumbers.('cellsno') = cellNumbers.('cellsno') - 1;
        if mutated
            cellNumbers.('genet') = cellNumbers.('genet') - 1;
        elseif NGresistant
            cellNumbers.('nongenet') = cellNumbers.('nongenet') - 1;
        else
            cellNumbers.('susceptibles') = cellNumbers.('susceptibles') - 1;
        end
    end
end

%% Find Boundary Coordinates Function

function [iMin,iMax,jMin,jMax] = boundaryCoors(lattice,previMin,previMax,prevjMin,prevjMax,xlatticeSize,ylatticeSize,limit)
    % BOUNDARYCOORS Gives the minimum and maximum rows containing a cell/drug pack and likewise for columns.
    %
    %   INPUT
    %       lattice is a matrix of ones and zeros, where a 1 states that a 
    %         cell/drug packet is a present and 0 states that no cell/drug
    %       previMin is the previous value of iMin
    %       previMax is the previous value of iMax
    %       prevjMin is the previous value of jMin
    %       prevjMax is the previous value of jMax
    %       xlatticeSize is the x dimension size of the lattice
    %       ylatticeSize is the y dimension size of the lattice
    %       limit is a bool:
    %           if true, the search area for the minimum and maximum row
    %             and column values is limited to within the previously set
    %             region + one extra row/column in each direction
    %           if false, the search area is the entire lattice
    
    if limit    
        % increase region by one in every direction
        if previMax + 1 > ylatticeSize
            previMax = ylatticeSize;
        else
            previMax = previMax + 1;
        end

        if previMin - 1 < 1
            previMin = 1;
        else
            previMin = previMin - 1;
        end

        if prevjMax + 1 > xlatticeSize
            prevjMax = xlatticeSize;
        else
            prevjMax = prevjMax + 1;
        end

        if prevjMin - 1 < 1
            prevjMin = 1;
        else
            prevjMin = prevjMin - 1;
        end
        
        % create matrix of the region set by previous min/max rows/columns
        limitedLattice = lattice(previMin:previMax,prevjMin:prevjMax);
        
        % find the min/max rows/columns of the restricted matrix and
        % determine the corresponding row/column number of that point in
        % the original lattice
        jMin = find(any(limitedLattice),1,'first') + (prevjMin - 1);
        jMax = find(any(limitedLattice),1,'last') + (prevjMin - 1);

        nRows = any(limitedLattice,2);
        iMin = find(nRows,1,'first') + (previMin - 1);
        iMax = find(nRows,1,'last') + (previMin - 1);
    else
        % find the min/max rows/columns of the lattice
        jMin = find(any(lattice),1,'first');
        jMax = find(any(lattice),1,'last');

        nRows = any(lattice,2);
        iMin = find(nRows,1,'first');
        iMax = find(nRows,1,'last'); 
    end
    
    % return -1 if the state is 0 everywhere, ie no cell/drug packet is present
    if isempty(jMin)
        jMin = -1;
    end
    if isempty(jMax)
        jMax = -1;
    end
    if isempty(iMin)
        iMin = -1;
    end
    if isempty(iMax)
        iMax = -1;
    end
end

%% Nutrient Movement Rule

function [i_n,j_n] = randomWalk(i,j,nSteps,xlatticeSize,ylatticeSize)
    % RANDOMWALK Performs a random walk of of a nutrient or drug packet
    %
    %   INPUT
    %       i is the current row where the packet is present
    %       j is the current column where the packet is present
    %       nSteps is the length of the random walk
    %       xlatticeSize is the x dimension of the lattice
    %       ylaticeSize is the y dimension of the lattice
    %
    %   OUTPUT
    %       i_n is the row where the packet is present after its random walk
    %       j_n is the column where the packet is present after its random walk
    
    step = 0;
    i_n = i;
    j_n = j;
    while step < nSteps
        A = sign(randn);
        i_n = i_n + A;
        if i_n < 1
            i_n = i_n + 2;
        elseif i_n > ylatticeSize
            i_n = i_n - 2;
        end
        B = sign(randn);
        j_n = j_n + B;
        if j_n < 1
            j_n = j_n + 2;
        elseif j_n > xlatticeSize
            j_n = j_n - 2;
        end
        step = step + 1;
    end
    
end
                
%% Count Number of Cells

function cellNumbers = countCells(lattices,iMin,jMin,iMax,jMax)
    % COUNTCELLS Counts the number of cells on the lattice and the number
    % of cells of each SNG state for 'SNG' mode
    %
    %   INPUT
    %       lattices is a structure array containing the following fields:
    %            -- 'cells' is a matrix of structure arrays with the following field:
    %                      + 'cell' is 0 if no cell is present, a Cell object if
    %                         a cell is present
    %            -- 'nutrients' is a matrix containing the number of nutrient packets
    %                at each point on the lattice
    %            -- 'states' is a matrix of ones and zeros, where a 1 states that a
    %                cell is a present and 0 states that no cell is present
    %            -- 'drugs' is a matrix containing the number of drug packets at each
    %                point on the lattice
    %       iMin is the smallest row number of the region within which the
    %         function will count cells
    %       iMax is the largest row number of the region within which the
    %         function will count cells
    %       jMin is the smallest column number of the region within which the
    %         function will count cells
    %       jMax is the largest column number of the region within which the
    %         function will count cells       
    %
    %   OUTPUT
    %       cellNumbers is a structrue array containing the total number of cells
    %           and for 'SNG' mode, the number of cells of each SNG type.
    %           Contains the following fields:
    %               -- 'cellsno' is the total number of living cells
    %               -- 'susceptibles' is the number of susceptible cells
    %               -- 'nongenet' is the number of non-genetically drug
    %                   resistant cells
    %               -- 'genet' is the number of genetically drug resistant cells

    cellsno = 0;
    susceptibles = 0;
    nongenet = 0;
    genet = 0;
    
    latticeSizes = size(lattices.('states'));
    xlatticeSize = latticeSizes(2);
    ylatticeSize = latticeSizes(1);
    if iMin - 1 < 1
        iMin = 2;
    end
    if iMax + 1 > ylatticeSize
        iMax = ylatticeSize - 1;
    end
    
    if jMin - 1 < 1
        jMin = 2;
    end
    if jMax + 1 > xlatticeSize
        jMax = xlatticeSize - 1;
    end

    for i = iMin-1:iMax+1
       for j = jMin-1:jMax+1
           if lattices.('states')(i,j) == 1
               cellsno = cellsno + 1;
               if lattices.('cells')(i,j).('cell').mutated
                   genet = genet + 1;
               elseif lattices.('cells')(i,j).('cell').NGresistant
                   nongenet = nongenet + 1;
               else
                   susceptibles = susceptibles + 1;
               end
           end
       end
    end
    cellNumbers = struct('cellsno',cellsno,'susceptibles',susceptibles,'nongenet',nongenet,'genet',genet);
end

%% Timestep or Cell Count for While Condition

function whileCondition = setEndCondition(FINAL_CELL_COUNT,FINAL_TIMESTEP,TIME_OR_COUNT,cellsno,time)
    % SETENDCONDITION Returns true or false based on whether the end condition has been met
    %
    %   INPUT
    %       FINAL_CELL_COUNT is the final cell count
    %       FINAL_TIMESTEP is the final timestep
    %       TIME_OR_COUNT is a str setting the type of end condition:
    %           'time' if the end condition is based on final time step
    %           'count' if the end condition is based on final cell count
    %           'both' if the program will end as soon as it reaches either
    %             the final time step or the final cell count, whichever comes first
    %       cellsno is the current number of living cells on the lattice
    %       time is the current time step
    %
    %   OUTPUT
    %       whileCondition is a bool indicating whether the end condition
    %       has been met
    
    if strcmp(TIME_OR_COUNT,'time')
        whileCondition = time<FINAL_TIMESTEP;
    elseif strcmp(TIME_OR_COUNT,'count')
        whileCondition = cellsno<FINAL_CELL_COUNT;
    elseif strcmp(TIME_OR_COUNT,'both')
        whileCondition = (time<FINAL_TIMESTEP) && (cellsno<FINAL_CELL_COUNT);
    end
end

%% Count Live Neighbours
function emptyNeighbour = emptyNeighbour(stateLattice,coordinates)
    % EMPTYNEIGHBOUR Returns whether there are any neighbouring sites around the given coordinates that are empty
    %
    %   INPUT
    %       stateLattice is a matrix of ones and zeros indicating whether a cell is present
    %       coordinates is 1x2 array containing the row then column of the coordinates 
    %   
    %   OUTPUT
    %       emptyNeighbour is true if a neighbouring site is present
    
    emptyNeighbour = false;
    latticeSizes = size(stateLattice);
    xlatticeSize = latticeSizes(2);
    ylatticeSize = latticeSizes(1);
    i = coordinates(1);
    j = coordinates(2);
    for I = -1:1
        for J = -1:1
            if I ~= 0 || J ~= 0
                if i+I <= ylatticeSize && i+I > 0 && j+J <= xlatticeSize && j+J > 0
                    if stateLattice(i+I,j+J) == 0
                        emptyNeighbour = true;
                    end
                end
            end
        end
    end
end