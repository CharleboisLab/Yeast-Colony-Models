function runSimulation
% RUNSIMULATION Runs colonySimulation.m once and gives measurements

%% Parameters

SEED_TYPE = 'leftSide';   % centre, leftSide, bothSides, cross, square, checkerboard

% END CONDITION SETTINGS
TIME_OR_COUNT = 'both';   % string to set end condition. 'count' to stop at FINAL_CELL_COUNT cells, 'time' to stop after FINAL_TIMESTEP timesteps, 'both' to stop after whichever it reaches first.
FINAL_CELL_COUNT = 562000;   % number of cells that colony should reach in order fo the program to stop.
FINAL_TIMESTEP = 1000000;   % number of timesteps the program will go to before stopping.
endConditions = struct('endCondition',TIME_OR_COUNT,'finalCellNo',FINAL_CELL_COUNT,'finalTime',FINAL_TIMESTEP);

% MUTATION SETTINGS
MUTATION_PROB = 0.00018;   % set probability of mutation occuring during during division.
SWITCH_UP_PROB = 0.0625;   % set probability of a cell's phenotype increasing in resistance: 
                           % in SNG mode it is the probability of a cell becoming non-genetically resistant, 
                           % in megaPlate mode it is the probability of a cell's phenotype number increasing
SWITCH_DOWN_PROB = 0.0035;   % set probability of a cell's phenotype decreasing in resistance: 
                             % in SNG mode it is the probability of a cell losing non-genetic resistance, 
                             % in megaPlate mode it is the probability of a cell's phenotype number decreasing
MUTATION_TYPE = 'megaPlate';   % SNG or megaPlate
INHERIT_PHENOTYPE = 1;   % 0 for false, 1 for true. For SNG mode, true means a non-genetically resistant cell's 
                         % daughter will also be non-genetically resistant. In megaPlate mode, true means a 
                         % daughter cell will inherit its mother cell's phenotype modifier 
N_DIVISION_PROB = 0.75;   % set probability of a non-genetically resistant cell dividing or surviving in the 
                         % presence of a static or cidal drug respectively, in SNG mode
mutationSettings = struct('mutationProb',MUTATION_PROB,'switchUpProb',SWITCH_UP_PROB,'switchDownProb',SWITCH_DOWN_PROB,'mutationType',MUTATION_TYPE,'inheritPhenotype',INHERIT_PHENOTYPE,'NDivisionProb',N_DIVISION_PROB);

% DISPLAY SETTINGS
DISPLAY_TYPE = 'show';   % show, save, or none

% DIFFUSION SETTINGS
NUTRIENT_STEPS = 2;   % number of steps that a nutrient packet takes during its random walk, correlates to the diffusion coefficient of the media.
DRUG_STEPS = 10;   % number of steps that a drug particle takes during its random walk in SNG mode, 
                   % correlates to the diffusion coefficient of the drugs through the media.
diffusionRates = struct('nutrientDiff',NUTRIENT_STEPS,'drugDiff',DRUG_STEPS);

% LATTICE SIZE SETTINGS
xlatticeSize = 100;   % x dimension size (number of columns)
ylatticeSize = 200;   % y dimension size (number of rows)
latticeSize = struct('x',xlatticeSize,'y',ylatticeSize);

% NUTRIENT SETTINGS
NUTRS_FOR_DIVISION = 1;   % number of nutrients necessary for a cell to divide.
startNutrs = 20;   % number of initial nutrient packets at each lattice site.

INIT_NUTRS_LAT = zeros(ylatticeSize,xlatticeSize) + startNutrs;   % create initial nutrient lattice

% DRUG SETTINGS
DRUGS_TYPE = 'cidal';   % static, cidal, or off

if strcmp(mutationSettings.('mutationType'),'SNG')
    startDrugs = 10;   % number of initial drug packets at each lattice site in SNG mode.
    INIT_DRUG_LAT = zeros(ylatticeSize,xlatticeSize) + startDrugs;   % create initial drug lattice for SNG mode
elseif strcmp(mutationSettings.('mutationType'),'megaPlate')
    % create initial drug lattice for megaPlate mode
    % in this setup the lattice is divided into five equal bands, where the first
    % band contains no drug, the second has 10 packets at each lattice
    % point, and from there each band increases by a factor of 10
    INIT_DRUG_LAT = zeros(ylatticeSize,xlatticeSize);

    bandWidth = round(xlatticeSize/5);
    for k = 1:4
        drugConc = (10^(k));
        for i = (k*bandWidth):((k+1)*bandWidth)
            for j = 1:ylatticeSize
                INIT_DRUG_LAT(j,i) = drugConc;
            end
        end
    end
end

lattice = repmat(struct('cell',0,'angle',0),ylatticeSize,xlatticeSize);
nutrientLattice = INIT_NUTRS_LAT;
drugLattice = INIT_DRUG_LAT;
stateLattice = zeros(ylatticeSize,xlatticeSize);
lattices = struct('cells',lattice,'nutrients',nutrientLattice,'drugs',drugLattice,'states',stateLattice);

% SEED LATTICE WITH CELLS

% leftSide - cells are located along every point of the leftmost column
% all cells are susceptible/genotype = 1/phenotype = 0
if strcmp(SEED_TYPE,'leftSide')
    for k = 1:ylatticeSize
        seedCell = Cell(NUTRS_FOR_DIVISION,mutationSettings);
        lattices.('cells')(k,1).('cell') = seedCell;
        lattices.('states')(k,1) = 1;
    end
% bothSides - cells are located along every point of the left and rightmost columns
% all cells are susceptible/genotype = 1/phenotype = 0
elseif strcmp(SEED_TYPE,'bothSides') 
    for k = 1:yLatticeSize
        seedCell = Cell(NUTRS_FOR_DIVISION,mutationSettings);
        lattices.('cells')(k,1).('cell') = seedCell;
        lattices.('states')(k,1) = 1;
        
        lattices.('cells')(k,xlatticeSize).('cell') = seedCell;
        lattices.('states')(k,xlatticeSize) = 1;
    end
% centre - vertical line of five cells in the centre column
% all cells are susceptible/genotype = 1/phenotype = 0
elseif strcmp(SEED_TYPE,'centre')
    xcentre = xlatticeSize/2;
    ycentre = ylatticeSize/2;
    for k = -2:2
        seedCell = Cell(NUTRS_FOR_DIVISION,mutationSettings);
        lattices.('cells')(ycentre+k,xcentre).('cell') = seedCell;
        lattices.('states')(ycentre+k,xcentre) = 1;
    end
% cross - vertical line of ten cells in the centre column and horizontal
% line of ten cells in the centre row
% ~10% of cells are non-genetically resistant, rest are susceptible,
% all cells have genotype = 1/phenotype = 0
elseif strcmp(SEED_TYPE,'cross')
    xcentre = xlatticeSize/2;
    ycentre = ylatticeSize/2;
    for k = -5:5
        seedCell1 = Cell(NUTRS_FOR_DIVISION,mutationSettings);
        seedCell2 = Cell(NUTRS_FOR_DIVISION,mutationSettings);
        
        p1 = rand();
        if p1 < 0.1
            seedCell1.NGresistant = true;
        end
        p2 = rand();
        if p2 < 0.1
            seedCell2.NGresistant = true;
        end
        
        lattices.('cells')(ycentre+k,xcentre).('cell') = seedCell1;
        lattices.('states')(ycentre+k,xcentre) = 1;
        lattices.('cells')(ycentre,xcentre+k).('cell') = seedCell2;
        lattices.('states')(ycentre,xcentre+k) = 1;
    end
% square - 20x20 square of cells in the centre of the lattice
% ~20% are non-genetically resistant, rest are susceptible,
% all cells have genotype = 1/phenotype = 0
elseif strcmp(SEED_TYPE,'square')
    xcentre = xlatticeSize/2;
    ycentre = ylatticeSize/2;
    for k = -10:10
        for j = -10:10
            seedCell = Cell(NUTRS_FOR_DIVISION,mutationSettings);
            
            p = rand();
            if p < 0.2
                seedCell.NGresistant = true;
            end
            
            lattices.('cells')(ycentre+k,xcentre+j).('cell') = seedCell;
            lattices.('states')(ycentre+k,xcentre+j) = 1;
        end
    end
% checkerboard - 40x40 checkerboard of cells in the centre of the lattice
% ~20% of cells are non-genetically resistant, rest are susceptible,
% all cells have genotype = 1/phenotype = 0
elseif strcmp(SEED_TYPE,'checkerboard')
    xcentre = xlatticeSize/2;
    ycentre = ylatticeSize/2;
    for k = -20:20
        for j = -20:20
            if mod((k+j),2) == 0
                seedCell = Cell(NUTRS_FOR_DIVISION,mutationSettings);

                p = rand();
                if p < 0.2
                    seedCell.NGresistant = true;
                end

                lattices.('cells')(ycentre+k,xcentre+j).('cell') = seedCell;
                lattices.('states')(ycentre+k,xcentre+j) = 1;
            end
        end
    end
end

% COLOUR MAP
if strcmp(MUTATION_TYPE,'SNG')
    map = [0 0 0   % background (BLACK)
       0.5 0.5 0.5   % drugs (GREY)
       1 1 1   % susceptible (WHITE)
       1 0.5 0   % non-genetically resistant (ORANGE)
       0 0 1];   % genetically resistant (BLUE)
elseif strcmp(MUTATION_TYPE,'megaPlate')
    map = [0 0 0   % background (BLACK)
       0.5 0.5 0.5   % drugs 1 (GREY)
       0.45 0.45 0.45   % drugs 2 (DARKER GREY)
       0.4 0.4 0.4   % drugs 3 (DARKER DARKER GREY)
       0.35 0.35 0.35   % drugs 4 (DARKER DARKER DARKER GREY)
       0.3 0.3 0.3   % drugs 5 (DARKER DARKER DARKER DARKER GREY)
       0.25 0.25 0.25   % drugs 6 (DARKEST GREY)
       1 1 1   % WT (WHITE)
       1 0 0   % mutant 1 (RED)
       1 0.5 0   % mutant 2 (ORANGE)
       1 1 0   % mutant 3 (YELLOW)
       0 0.5 0   % mutant 4 (GREEN)
       0 0 1   % mutant 5 (BLUE)
       0.5 0 0.5];   % mutants 6+ (PURPLE)
end
FOLDER_NAME = 'images';
displaySettings = struct('displayType',DISPLAY_TYPE,'map',map,'folderName',FOLDER_NAME);

%% Run
disp('Mutation Settings:')
disp(mutationSettings)

disp('Drug Type:')
disp(DRUGS_TYPE)

disp('Lattice Size:')
disp(num2str([xlatticeSize ylatticeSize]))

[phenotypeAves,genotypeAves,appearanceTimes,establishmentTimes,fixationTimes,cellGenotypeNos,cellNumbers,cellCountTable,time] = colonySimulation(diffusionRates,endConditions,mutationSettings,latticeSize,lattices,DRUGS_TYPE,displaySettings);

%% Display Results

disp('time (timesteps)')
disp(time)

if strcmp(MUTATION_TYPE,'megaPlate')
    % Display cells counts and mutation time benchmarks
    
    disp('final cell count:')
    disp(cellNumbers.('cellsno'))
    disp('final genotype counts:')
    cellArray = num2cell(cellGenotypeNos);
    table = cell2table(cellArray,'VariableNames',{'WT','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12','G13','G14','G15'});
    disp(table)
    disp('first appearance times:')
    cellArray = num2cell(appearanceTimes);
    table = cell2table(cellArray,'VariableNames',{'WT','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12','G13','G14','G15'});
    disp(table)
    disp('establishment times:')
    cellArray = num2cell(establishmentTimes);
    table = cell2table(cellArray,'VariableNames',{'WT','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12','G13','G14','G15'});
    disp(table)
    disp('fixation times:')
    cellArray = num2cell(fixationTimes);
    table = cell2table(cellArray,'VariableNames',{'WT','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12','G13','G14','G15'});
    disp(table)
    
    % Display genotype and phenotype averages over time
    
    figure
    phenotypes = phenotypeAves(1,:);
    genotypes = genotypeAves(1,:);
    t = phenotypeAves(2,:);
    plot(t,phenotypes,t,genotypes)
    xlabel('time (timesteps)')
    ylabel('mean genotype and phenotype levels')
    legend({'mean phenotype level','mean genotype level'},'Location','southwest')

elseif strcmp(MUTATION_TYPE,'SNG')
    disp('final cell counts')
    disp(cellNumbers)
    
    cellCounts = table2array(cellCountTable);
    timestep = cellCounts(:,1);
    cellCountS = cellCounts(:,2);
    cellCountN = cellCounts(:,3);
    cellCountG = cellCounts(:,4);
    
    figure
    plot(timestep,cellCountN,timestep,cellCountG,timestep,cellCountS,'LineWidth',2)
    xlabel('time (timesteps)')
    ylabel('cell counts')
    legend({'non-genetic resistant','genetic resistant','susceptibles'},'Location','southwest')
end


end


