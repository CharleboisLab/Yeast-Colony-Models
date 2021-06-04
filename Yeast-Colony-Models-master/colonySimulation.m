function [cellsno, time, stateLattice] = colonySimulation(nSteps,FINAL_CELL_COUNT,FINAL_TIMESTEP,TIME_OR_COUNT,START_NUTRS,NUTRS_FOR_BUDDING,AXIAL_FRAC,MAGNETIC_FIELD,MF_STRENGTH,MIN_ANGLE,MAX_ANGLE,UNIPOLAR_ON,MUTATION_ON,MUTATION_PROB,DISPLAY_IMAGE)
% colonySimulation models the growth of a yeast colony under the
% influence of both nutrient concentrations and a static magnetic field.

p = inputParser;
% % parameters
% nSteps is the number of steps that a nutrient packet takes during its 
% random walk, correlates to the diffusion coefficient of the media.
addRequired(p,'nSteps',@(x)validateattributes(x,{'numeric'},{'nonnegative'}))

% FINAL_CELL_COUNT is the number  of cells that colony should reach in 
% order fo the program to stop.
addRequired(p,'FINAL_CELL_COUNT',@(x)validateattributes(x,{'numeric'},{'nonnegative'}))

% FINAL_TIMESTEP is the number of timesteps until the program stops.
addRequired(p,'FINAL_TIMESTEP',@(x)validateattributes(x,{'numeric'},{'nonnegative'}))

% TIME_OR_COUNT is a str which determines whether the end condition will
% be based on timesteps of final cell count: 
% 'time' = timestep, 'count' = cell count
addRequired(p,'TIME_OR_COUNT',@(x)validateattributes(x,{'string','char'},{'nonempty'}))

% START_NUTRS is the number of nutrient packets at each lattice site at the
% beginning.
addRequired(p,'START_NUTRS',@(x)validateattributes(x,{'numeric'},{'>',0}))

% NUTRS_FOR_BUDDING is the number of nutrients necessary for a cell to bud.
addRequired(p,'NUTRS_FOR_BUDDING',@(x)validateattributes(x,{'numeric'},{'>',0}))

% AXIAL_FRAC is the overall fraction of cells budding normally which bud 
% axially: 0 = average diploid colony, 0.6 = average haploid colony.
addRequired(p,'AXIAL_FRAC',@(x)validateattributes(x,{'numeric'},{'>=',0,'<=',1}))

% MAGNETIC_FIELD is the vector for the direction of the magnetic field, 
% according to x-y coordinates, not row-column, if there is no magnetic field set to [0 0].
addRequired(p,'MAGNETIC_FIELD',@(x)validateattributes(x,{'numeric'},{'size',[1,2]}))

% MF_STRENGTH is the probability the magnetic field bias will be applied, 
% used to control the strength of the magnetic field.
addRequired(p,'MF_STRENGTH',@(x)validateattributes(x,{'numeric'},{'>=',0}))

% MIN_ANGLE is the the minimum angle from the magnetic field of the range 
% of angles the MF biases budding towards.
addRequired(p,'MIN_ANGLE',@(x)validateattributes(x,{'numeric'},{'>=',0,'<=',180}))

% MAX_ANGLE is the maximum angle from the magnetic field of the range of 
% angles the MF biases budding towards.
addRequired(p,'MAX_ANGLE',@(x)validateattributes(x,{'numeric'},{'>=',0,'<=',180}))

% UNIPOLAR_ON is a bool indicating whether the colony will bud in a
% filamentous pattern in low nutrient conditions.
addRequired(p,'UNIPOLAR_ON',@(x)validateattributes(x,{'logical'},{'nonempty'}))

% MUTATION_ON is a bool indicating whether the cells might mutate.
addRequired(p,'MUTATION_ON',@(x)validateattributes(x,{'logical'},{'nonempty'}))

% MUTATION_PROB is the probability of a mutation occuring during budding.
addRequired(p,'MUTATION_PROB',@(x)validateattributes(x,{'numeric'},{'>=',0,'<=',1}))

% DISPLAY_COLONY is a bool indicating whether the colony should be
% visualised
addRequired(p,'DISPLAY_IMAGE',@(x)validateattributes(x,{'logical'},{'nonempty'}))

parse(p,nSteps,FINAL_CELL_COUNT,FINAL_TIMESTEP,TIME_OR_COUNT,START_NUTRS,NUTRS_FOR_BUDDING,AXIAL_FRAC,MAGNETIC_FIELD,MF_STRENGTH,MIN_ANGLE,MAX_ANGLE,UNIPOLAR_ON,MUTATION_ON,MUTATION_PROB,DISPLAY_IMAGE)

% Assertions
list = [0 1 -1];
val1chk = ismember(MAGNETIC_FIELD(1),list);
val2chk = ismember(MAGNETIC_FIELD(2),list);
assert((val1chk && val2chk), 'MAGNETIC_FIELD must be one of the following: [0 0], [0 1], [0 -1], [1 0], [1 1], [1 -1], [-1 0], [-1 1], [-1 -1]')
assert(MIN_ANGLE < MAX_ANGLE, 'MIN_ANGLE must be less than MAX_ANGLE')
assert((strcmp(TIME_OR_COUNT,'time') || strcmp(TIME_OR_COUNT,'count')),"TIME_OR_COUNT must be set to 'time' or 'count' strings")
if strcmp(TIME_OR_COUNT,'time')
    assert(FINAL_TIMESTEP > 0,'FINAL_TIMESTEP should be greater than 0')
else
    assert(FINAL_CELL_COUNT > 0,'FINAL_CELL_COUNT should be greater than0')
end

%% Initialisation

deltaUnipolarFrac = 1/(NUTRS_FOR_BUDDING*8);   % change in unipolarFrac for each nutrient packet available to the cell.

if strcmp(TIME_OR_COUNT,'count')
    % Determine the lattice size necessary for the set FINAL_CELL_COUNT
    % Determine the largest possible diagonal according to colony with
    % area = 2*FINAL_CELL_COUNT and elongation 2
    maxDiag = 2*(sqrt((2*FINAL_CELL_COUNT)/(0.5*pi))+250);
end

if strcmp(TIME_OR_COUNT,'time')
    % Determine the lattice size necessary for the set FINAL_TIMESTEP
    growthRate = 10000/320000;   % to match the lattice size of final cell count = 10 000 when final timestep = 320 000
    % Determine the largest possible diagonal according to colony with area
    % 2*growthRate*FINAL_TIMESTEP and elongation 2
    maxDiag = 2*(sqrt((2*growthRate*FINAL_TIMESTEP)/(0.5*pi))+250);
end

% Get the length of the lattice size using the diagonal and Pythagoras
latticeEdge = sqrt(((maxDiag^2)/2));
% Round the length to the closest even integer
latticeSize = ceil(latticeEdge);
if mod(latticeSize,2)==1
    latticeSize = latticeSize+1;
end
if latticeSize < 500
    latticeSize = 500;
end

% Make lattice of CellWithNutrients objects.
% See CellWithNutrients.m file for details of what information is saved.

% Build lattice.
cellLattice = [CellWithNutrients(0,[0 0],START_NUTRS,deltaUnipolarFrac)];
cellLattice = repelem(cellLattice,latticeSize,latticeSize);

% stateLattice is a lattice of ones and zeros indicating only state.
stateLattice = zeros(latticeSize);

% Set properties of seed cells.
centre = [latticeSize/2 latticeSize/2];

% The bud scar (position 2 in CellWithNutrients class constructor) is chosen
% randomly and the location of the seed cell can both be varied, if desired.
initialBudScars = [0 1;0 -1;1 1;1 0;1 -1;-1 1;-1 0;-1 -1];

if MUTATION_ON
    mutatedStateLattice = ones(latticeSize);
end

for k = -2:2
    i = randi([1,8]);
    initialScar = [initialBudScars(i,1),initialBudScars(i,2)];
    seedCell = CellWithNutrients(1,initialScar,START_NUTRS,deltaUnipolarFrac);
    cellLattice(centre(1)+k,centre(2)) = seedCell;
    stateLattice(centre(1)+k,centre(2)) = 1;
    if MUTATION_ON
        mutatedStateLattice(centre(1)+k,centre(2)) = 2;
    end
end

% Sets the range of cells which will be considered under the rules as the
% cells immediate surrounding the start cell.
[iMin,iMax,jMin,jMax] = boundaryCoors(stateLattice);

% Display conditions
if DISPLAY_IMAGE == true
    disp('Diffusion Steps:')
    disp(num2str(nSteps))

    disp('Initial Concentration:')
    disp(num2str(START_NUTRS))

    disp('Nutrients Needed to Bud:')
    disp(num2str(NUTRS_FOR_BUDDING))

    disp('Magnetic Field:')
    disp(num2str(MAGNETIC_FIELD))

    disp('Magnetic Field Strength:')
    disp(num2str(MF_STRENGTH))

    disp('Axial Frac:')
    disp(num2str(AXIAL_FRAC))

    % Visualise initial states.

    disp('Press any key...')
    imshow(stateLattice);
    pause
end

time=0; 

%% Run
cellsno = countCells(stateLattice,iMin,jMin,iMax,jMax);   % number of cells alive on the lattice
prevno = 0;

whileCondition = endCondition(FINAL_CELL_COUNT,FINAL_TIMESTEP,TIME_OR_COUNT,cellsno,time);

% Run program until the end condition is met
while whileCondition  
    % Apply nutrient movement rule
    containsNutrient = false;
    while ~containsNutrient
    i_r = randi([3 (latticeSize-3)],1);
    j_r = randi([3 (latticeSize-3)],1);
        if cellLattice(i_r,j_r).nutrientConcentration > 0
            containsNutrient = true;
            [i_n,j_n] = randomWalk(i_r,j_r,nSteps);
            cellLattice(i_r,j_r).nutrientConcentration = cellLattice(i_r,j_r).nutrientConcentration - 1;
            cellLattice(i_n,j_n).nutrientConcentration = cellLattice(i_n,j_n).nutrientConcentration + 1;
        end
    end

    aliveCellChosen = false;
    while ~aliveCellChosen
        % Randomly find a lattice site containing a cell, then apply
        % budding rules.
        i_r = randi([iMin iMax],1);
        j_r = randi([jMin jMax],1);

        if stateLattice(i_r,j_r) == 1
           % cell consumes a nutrient if one is present
           cellLattice(i_r,j_r) = cellLattice(i_r,j_r).updateNutrients();
           % update unipolar fraction according nutrient concentration
           cellLattice(i_r,j_r) = cellLattice(i_r,j_r).updateUnipolarFrac();
           % check if the cell is able to bud
           cellLattice(i_r,j_r) = cellLattice(i_r,j_r).updateBuddability(NUTRS_FOR_BUDDING);
           
           if cellLattice(i_r,j_r).buddable
                % Apply budding rules.

                % [i_d j_d] are the coordinates of the new daughter cell 

                % [i_mc j_mc] are the coordinates of the cell that gets moved
                % to make way for the new daughter cell or for the
                % extension of the daughter cell in the case of
                % unipolar budding.
                % ([i_mc j_mc] = [-1 -1] if not applicable)

                [didBud,unipolarBud,i_d,j_d,i_mc,j_mc,motherBudScar,daughterBudScar,extraBudScar] = buddingRules(MAGNETIC_FIELD,MF_STRENGTH,MIN_ANGLE,MAX_ANGLE,cellLattice,stateLattice,AXIAL_FRAC,UNIPOLAR_ON,i_r,j_r);
                
                if didBud
                    cellLattice(i_d,j_d).budded = true;
                    if i_mc ~= -1 && cellLattice(i_mc,j_mc).state
                        cellLattice(i_d,j_d).budded = true;
                    end
                    % mutate or inherit mutations
                    if MUTATION_ON
                        if ~unipolarBud
                            % if the new daughter cell is pushing a mutated
                            % cell out of the way, keep it mutated at its
                            % new location
                            if i_mc ~= -1 && cellLattice(i_d,j_d).state == 1 && cellLattice(i_d,j_d).mutated == true
                                cellLattice(i_mc,j_mc).mutated = true;
                            end
                        end
                        
                        % inherit mutation
                        if cellLattice(i_r,j_r).mutated == true
                            cellLattice(i_d,j_d).mutated = true;
                            if unipolarBud && i_mc ~= -1
                                cellLattice(i_mc,j_mc).mutated = true;
                            end
                        end
                        
                        % mutate
                        p = rand();
                        if p<MUTATION_PROB
                            cellLattice(i_d,j_d).mutated = true;
                            if unipolarBud && i_mc ~= -1  % if the daughter cell is elongated (composed of two lattice points) and has mutated, the tip of the daughter cell is mutated
                                cellLattice(i_mc,j_mc).mutated = true;
                            end
                        end
                    end
                    
                    % update mother cell
                    cellLattice(i_r,j_r).budScar = motherBudScar;
                    cellLattice(i_r,j_r).consumedNutrients = 0;

                    % update daughter cell
                    cellLattice(i_d,j_d).state = 1;
                    stateLattice(i_d,j_d) = 1;
                    cellLattice(i_d,j_d).budScar = daughterBudScar;

                    % update extra cell
                    if i_mc ~= -1   % if [i_mc j_mc] is a coordinate
                        cellLattice(i_mc,j_mc).state = 1;
                        stateLattice(i_mc,j_mc) = 1;
                        cellLattice(i_mc,j_mc).budScar = extraBudScar;
                    end
                cellsno = cellsno + 1; 
                [iMin,iMax,jMin,jMax] = boundaryCoors(stateLattice);
                end
           end
           aliveCellChosen = true;
        end
    end
    time = time + 1;
    
    whileCondition = endCondition(FINAL_CELL_COUNT,FINAL_TIMESTEP,TIME_OR_COUNT,cellsno,time);

    % Display stateLattice.
    if DISPLAY_IMAGE == true
        if cellsno - prevno > 100 || cellsno - prevno > prevno  
            if MUTATION_ON
                for i = iMin:iMax
                    for j = jMin:jMax
                        if cellLattice(i,j).state == 1
                            if cellLattice(i,j).mutated == true
                                mutatedStateLattice(i,j) = 3;
                            else
                                mutatedStateLattice(i,j) = 2;
                            end
                        else
                            mutatedStateLattice(i,j) = 1;
                        end
                    end
                end
                map = [0 0 0
                       1 1 1
                       1 0 0];
                image = ind2rgb(mutatedStateLattice,map);
                imshow(image);  
            else
                imshow(stateLattice); 
            end           
            drawnow;
            prevno = cellsno;
        end
    end 
end 
if DISPLAY_IMAGE == true
    if MUTATION_ON
        for i = iMin:iMax
            for j = jMin:jMax
                if cellLattice(i,j).state == 1
                    if cellLattice(i,j).mutated == true
                        mutatedStateLattice(i,j) = 3;
                    else
                        mutatedStateLattice(i,j) = 2;
                    end
                else
                    mutatedStateLattice(i,j) = 1;
                end
            end
        end
        map = [0 0 0
               1 1 1
               1 0 0];
        image = ind2rgb(mutatedStateLattice,map);
        imshow(image);  
    else
        imshow(stateLattice); 
    end 
    drawnow;
    disp('All done!')
end

end

%% Model Rules Function

function [didBud,unipolarBud,i_d,j_d,i_mc,j_mc,motherBudScar,daughterBudScar,extraBudScar] = buddingRules(magneticField,MFStrength,minAngle,maxAngle,cellLattice,stateLattice,axialFrac,unipolar_on,i,j)
% buddingRules describes the rules of budding the model.
% Gives coordinates of the new daughter cell and moved cell, if applicable,
% and bud scar of mother and daughter cells.

    aliveCellsno = aliveCells(stateLattice,i,j);
    
    neighbours = [i-1 j+1; i j+1; i+1 j+1; i-1 j; i+1 j; i-1 j-1; i j-1; i+1 j-1];
    
    i_d = 0;
    j_d = 0;
    i_mc = -1;
    j_mc = -1;
    motherBudScar = cellLattice(i,j).budScar;
    daughterBudScar = [0 0];
    extraBudScar = [0 0];
    didBud = false;
    unipolarBud = false;
    
    if aliveCellsno<8
        budScar = cellLattice(i,j).budScar;
        
        emptyCells = emptyCellArray(stateLattice,i,j,aliveCellsno);

        angleArray = createAngleArray(neighbours,budScar,i,j);

        emptyAngleArray = createAngleArray(emptyCells,budScar,i,j);
        
        buddable = true;

        % Determine what budding pattern the cell will follow.
        p1 = rand();
        if p1 < cellLattice(i,j).unipolarFrac || (unipolar_on == false)
            if cellLattice(i,j).budded == false
                p2 = rand();
                if p2<axialFrac
                    [didBud,i_d,j_d,i_mc,j_mc] = axialBuddingRule(magneticField,MFStrength,minAngle,maxAngle,stateLattice,emptyAngleArray,angleArray,neighbours,emptyCells,i,j);
                else
                    [didBud,i_d,j_d,i_mc,j_mc] = polarBuddingRule(magneticField,MFStrength,minAngle,maxAngle,stateLattice,emptyAngleArray,angleArray,neighbours,emptyCells,i,j);
                end
            else
                p3 = rand();
                if p3 < (1+axialFrac)/2
                    [didBud,i_d,j_d,i_mc,j_mc] = axialBuddingRule(magneticField,MFStrength,minAngle,maxAngle,stateLattice,emptyAngleArray,angleArray,neighbours,emptyCells,i,j);
                else
                    [didBud,i_d,j_d,i_mc,j_mc] = polarBuddingRule(magneticField,MFStrength,minAngle,maxAngle,stateLattice,emptyAngleArray,angleArray,neighbours,emptyCells,i,j);
                end
            end
         else
            [unipolarBud,i_d,j_d,i_mc,j_mc] = unipolarBuddingRule(magneticField,MFStrength,minAngle,maxAngle,stateLattice,emptyAngleArray,angleArray,emptyCells,i,j);
            if unipolarBud
                didBud = true;
            else
                buddable = false;
            end
        end

         % Set the new bud scars for both mother and daughter.
         if buddable
             if ~unipolarBud
                daughterBudScar = [i-i_d,j-j_d];
                motherBudScar = [i_d-i,j_d-j];
                if i_mc ~= -1
                    extraBudScar = cellLattice(i_d,j_d).budScar;
                end
             else
                extraBudScar = [i_d-i_mc,j_d-j_mc];
             end
         end
    end
end

%% Axial Budding Function

function [didBud,i_d,j_d,i_mc,j_mc] = axialBuddingRule(magneticField,MFStrength,minAngle,maxAngle,stateLattice,emptyAngleArray,angleArray,neighbours,emptyCells,i,j)

% axialBuddingRule describes how a cell will bud if it is budding axially.

% Returns: didBud - logical value indicating whether the cell has bud
%          i_d,j_d - coordinates of daughter cell
%          i_mc,j_mc - coordinates of extra cell, -1,-1 if there is none
    
    didBud = true;

    % Determine the minimum angle, which should be 45.
    nearbyAngle = min(angleArray(angleArray>44));
    i_mc = -1;
    j_mc = -1;
    % When at least one location 45 degrees away from the previous bud site
    % is empty, the cell will bud there, no cell will be moved out of the
    % way.
    if ismember(nearbyAngle,emptyAngleArray)
        % Determine the indices for every instance of the minimum angle.
        nearbyAngleInds = find(emptyAngleArray == nearbyAngle);
        
        % Add EMF bias condition.
        p = rand();
        if p < MFStrength
            indices = biasTowardEMF(nearbyAngleInds,emptyCells,magneticField,minAngle,maxAngle,i,j);
            cellList = emptyCells;
            if isempty(indices)
                if MFStrength <= 1
                    indices = nearbyAngleInds;
                else
                    nearbyAngleInds = find(angleArray == nearbyAngle);
                    indices = biasTowardEMF(nearbyAngleInds,neighbours,magneticField,minAngle,maxAngle,i,j);
                    cellList = neighbours;
                end
            end
        else
            indices = nearbyAngleInds;
            cellList = emptyCells;
        end
        
        % Randomly select one of these instances as the location for the
        % new daughter cell.
        randInd = randi(length(indices),1);
        randEmptyCellInd = indices(randInd);
        i_d = cellList(randEmptyCellInd,1);
        j_d = cellList(randEmptyCellInd,2);
    % If both sites 45 degrees away from the pevious bud site are full,
    % the cell will bud into one of the sites, pushing the existing cell out
    % of the way. If there are no empty site surround the existing cell,
    % nothing will happen.
    else
        % Determine the indices for every instance of the minimum angle.
        nearbyAngleInds = find(angleArray == nearbyAngle);
        
        % Add EMF bias condition.
        p = rand();
        if p < MFStrength
            indices = biasTowardEMF(nearbyAngleInds,neighbours,magneticField,minAngle,maxAngle,i,j);
        else
            indices = nearbyAngleInds;
        end
        
        % Randomly select one of these instances as the location for the
        % new daughter cell.
        randInd = randi(length(indices),1);
        randNeighboursInd = indices(randInd);
        i_d = neighbours(randNeighboursInd,1);
        j_d = neighbours(randNeighboursInd,2);
    end
    if stateLattice(i_d,j_d) == 1
        % Randomly select an empty site surround the preexisting cell to
        % move into.
        aliveCellsno = aliveCells(stateLattice,i_d,j_d);
        otherEmptyCells = emptyCellArray(stateLattice,i_d,j_d,aliveCellsno);
        
        if ~isempty(otherEmptyCells)
            [numRows,~] = size(otherEmptyCells);
            Inds = [1:numRows];

            % Add EMF bias condition.
            p = rand();
            if p < MFStrength
                indices = biasTowardEMF(Inds,otherEmptyCells,magneticField,minAngle,maxAngle,i,j);
                if isempty(indices)
                    indices = Inds;
                end
            else
                indices = Inds;
            end
           
            randInd = randi(length(indices),1);
            randNeighboursInd = indices(randInd);
            i_mc = otherEmptyCells(randNeighboursInd,1);
            j_mc = otherEmptyCells(randNeighboursInd,2);
        else
            didBud = false;
        end
    end
end

%% Polar Budding Function

function [didBud,i_d,j_d,i_mc,j_mc] = polarBuddingRule(magneticField,MFStrength,minAngle,maxAngle,stateLattice,emptyAngleArray,angleArray,neighbours,emptyCells,i,j)

% polarBuddingRule describes how a cell will bud when it follows a
% polar budding pattern.

% Returns: didBud - logical value indicating whether the cell has bud
%          i_d,j_d - coordinates of daughter cell
%          i_mc,j_mc - coordinates of extra cell, -1,-1 if there is none

    i_mc = -1;
    j_mc = -1;
    
    didBud = true;
    
    oppAngle = max(angleArray);
    secondAngle = max(angleArray(angleArray < oppAngle));
    % oppAngle should be 180 degrees.
    % secondAngle should be 135 degrees.
    
    % If at least one site which is either 135 or 180 degrees away from the
    % previous bud site is empty, the cell will bud into one of those sites
    % and no existing cells will be affected.
    if ismember(oppAngle,emptyAngleArray) || ismember(secondAngle,emptyAngleArray)
       % determine the indices for every instance of oppAngle and
       % secondAngle in emptyAngleArray
        maxAngleInds = find(emptyAngleArray == oppAngle | emptyAngleArray == secondAngle);
        
        % Add EMF bias condition
        p = rand();
        if p < MFStrength
            indices = biasTowardEMF(maxAngleInds,emptyCells,magneticField,minAngle,maxAngle,i,j);
            cellList = emptyCells;
            if isempty(indices)
                if MFStrength <= 1
                    indices = maxAngleInds;
                else
                    maxAngleInds = find(angleArray == oppAngle | angleArray == secondAngle);
                    indices = biasTowardEMF(maxAngleInds,neighbours,magneticField,minAngle,maxAngle,i,j);
                    cellList = neighbours;
                end
            end
        else
            indices = maxAngleInds;
            cellList = emptyCells;
        end
        
        % Randomly select one of these instances as the location for the
        % new daughter cell.
        randInd = randi(length(indices),1);
        randEmptyCellsInd = indices(randInd);
        i_d = cellList(randEmptyCellsInd,1);
        j_d = cellLIst(randEmptyCellsInd,2);
    % If all three sites which are 135 or 180 degrees away from the
    % previous bud site are full, a random site will be selected and the
    % existing cell will be moved into a random surrounding empty site. If
    % no such empty sites exist, nothing happens.
    else
        % determine the indices for every instance of oppAngle and
        % secondAngle in angleArray
        maxAngleInds = find(angleArray == oppAngle | secondAngle);
        
        % Add EMF bias condition.
        p = rand();
        if p < MFStrength
            indices = biasTowardEMF(maxAngleInds,neighbours,magneticField,minAngle,maxAngle,i,j);
        else
            indices = maxAngleInds;
        end
        
        % randomly select one of these instances as the location for the
        % new daughter cell
        randInd = randi(length(indices),1);
        randNeighboursInd = (indices(randInd));
        i_d = neighbours(randNeighboursInd,1);
        j_d = neighbours(randNeighboursInd,2);
    end
    
    if stateLattice(i_d,j_d) == 1
        % Randomly select an empty site for the existing cell to move into.
        aliveCellsno = aliveCells(stateLattice,i_d,j_d);
        otherEmptyCells = emptyCellArray(stateLattice,i_d,j_d,aliveCellsno);
        if ~isempty(otherEmptyCells)
            [numRows,~] = size(otherEmptyCells);
            Inds = [1:numRows];
            
            % Add EMF bias condition.
            p = rand();
            if p < MFStrength
                indices = biasTowardEMF(Inds,otherEmptyCells,magneticField,minAngle,maxAngle,i,j);
                if isempty(indices)
                    indices = Inds;
                end
            else
                indices = Inds;
            end
           
            randInd = randi(length(indices),1);
            randNeighboursInd = indices(randInd);
            i_mc = otherEmptyCells(randNeighboursInd,1);
            j_mc = otherEmptyCells(randNeighboursInd,2);
        else
            didBud = false;
        end
    end
end

%% Unipolar Budding Function

function [unipolarBud,i_d,j_d,i_mc,j_mc] = unipolarBuddingRule(magneticField,MFStrength,minAngle,maxAngle,stateLattice,emptyAngleArray,angleArray,emptyCells,i,j)

% unipolarBuddingRule describes how a cell will bud when it follows a
% unipolar budding pattern.

% Returns: didBud - logical value indicating whether the cell has bud
%          i_d,j_d - coordinates of back end of daughter cell
%          i_mc,j_mc - coordinates of front end of daughter cell

    i_mc = -1;
    j_mc = -1;
    i_d = -1;
    j_d = -1;
    unipolarBud = false;
    oppAngle = max(angleArray);
    secondAngle = max(angleArray(angleArray < oppAngle));
    % oppAngle should be 180 degrees.
    % secondAngle should be 135 degrees.

    % If at least one site which is either 135 or 180 degrees away from the
    % previous bud site is empty, the cell will bud into one of those sites
    % and no existing cells will be affected.
    if ismember(oppAngle,emptyAngleArray) || ismember(secondAngle,emptyAngleArray)
        possibleSite = false;

       % determine the indices for every instance of oppAngle and
       % secondAngle in emptyAngleArray
        maxAngleInds = find(emptyAngleArray == oppAngle | secondAngle);

        % Add EMF bias condition.
        p = rand();
        if p < MFStrength
            indices = biasTowardEMF(maxAngleInds,emptyCells,magneticField,minAngle,maxAngle,i,j);
        else
            indices = maxAngleInds;
        end
        while ~possibleSite && ~isempty(indices)
            % randomly select one of these instances as the location for the
            % new daughter cell
            randInd = randi(length(indices),1);
            randEmptyCellsInd = indices(randInd);
            i_d = emptyCells(randEmptyCellsInd,1);
            j_d = emptyCells(randEmptyCellsInd,2);

            % The cell will bud into two cells in a row, since cells undergoing
            % unipolar budding are elongated.
            I = i_d-i;
            J = j_d-j;
                i_mc = i_d + I;
                j_mc = j_d + J;
                if stateLattice(i_mc,j_mc)==0
                    possibleSite = true;
                else
                    indices = indices(indices~=randEmptyCellsInd);
                end
        end 
        if possibleSite
            unipolarBud = true;
        end
    end
end

%% Count Alive Neighbours Function

function aliveCellsno = aliveCells(stateLattice,i,j)

    % aliveCells counts the number of alive cells surround a cell at (i,j).

    aliveCellsno = 0;
    aliveCellsno = aliveCellsno + stateLattice(i-1,j-1) + stateLattice(i-1,j) + stateLattice(i-1,j+1);
    aliveCellsno = aliveCellsno + stateLattice(i,j-1) + stateLattice(i,j+1);
    aliveCellsno = aliveCellsno + stateLattice(i+1,j-1) + stateLattice(i+1,j) + stateLattice(i+1,j+1);

end

%% Empty Cells Function

function emptyCells = emptyCellArray(stateLattice,i,j,aliveCellsno)

    % emptyCellArray gives an array containing the indices of every empty
    % cell surrounding a cell at (i,j) with dimensions emptySpots x 2.
    
    ind = 1; 
    emptySpots = 8 - aliveCellsno;
    emptyCells = zeros(emptySpots,2);
    for I = [-1 0 1]
        for J = [-1 0 1]
            if I == 0 && J == 0
            else 
                if stateLattice(i+I, j+J) == 0
                    emptyCells(ind,1) = i+I;
                    emptyCells(ind,2) = j+J;
                    ind = ind + 1;
                end
            end
        end
    end      
end

%% Angle Array Function

function angleArray = createAngleArray(cellArray,vector,i,j)

    % createAngleArray creates an array containing the angle between each
    % lattice site in cellArray and the lattice site (i,j).
    
    [numRows,~] = size(cellArray);

    % create array containing angle between the vector of the
    % mother's bud scar and the vector from the mother cell to each
    % lattice site in cellArray
    angleArray = zeros(numRows,1);
    
    for k = 1:numRows
        i_2 = cellArray(k,1) - i;
        j_2 = cellArray(k,2) - j;
        dotProd = vector*[i_2;j_2];
        norm1 = norm(vector);
        norm2 = norm([i_2 j_2]);
        norms = norm1*norm2;
        normDot = dotProd/norms;
        angle = acosd(normDot);
        angleArray(k) = angle;
    end
    
end


%% Angle Array for EMF Bias Function

function angleArray = createEMFAngleArray(cellArray,vector,i,j)

    % createAngleArray creates an array containing the angle between each
    % lattice site in cellArray and the lattice site (i,j).
    
    [numRows,~] = size(cellArray);

    % create array containing angle between the vector of the
    % mother's bud scar and the vector from the mother cell to each
    % lattice site in cellArray
    angleArray = zeros(numRows,1);
    for k = 1:numRows
        if (cellArray(k,1) + cellArray(k,2)) ~= 0
            i_2 = cellArray(k,1) - i;
            j_2 = cellArray(k,2) - j;
            norm2 = norm([i_2 j_2]);
            dotProd = vector*[i_2;j_2];
            norm1 = norm(vector);
            norms = norm1*norm2;
            normDot = dotProd/norms;
            angle = acosd(normDot);
            angleArray(k) = angle;
            
        else
            angleArray(k) = 0;
        end
    end
    
end


%% Find Boundary Coordinates Function

function [iMin,iMax,jMin,jMax] = boundaryCoors(stateLattice)

    % boundaryCoors gives the minimum and maximum rows containing a cell
    % and likewise for columns.

    jMin = find(any(stateLattice),1,'first')-5;
    jMax = find(any(stateLattice),1,'last')+5;
  
    nRows = any(stateLattice,2);
    iMin = find(nRows,1,'first')-5;
    iMax = find(nRows,1,'last')+5;
end

%% Nutrient Movement Rule

function [i_n,j_n] = randomWalk(i,j,nSteps)

    % nutrientRandomWalk performs a random walk of nSteps of a nutrient packet at (i,j)
    
    step = 0;
    i_n = i;
    j_n = j;
    while step < nSteps
        A = sign(randn);
        i_n = i + A;
        j_n = j + A;
        step = step + 1;
    end
    
end
                
%% Count Number of Cells

function nCells = countCells(stateLattice,iMin,jMin,iMax,jMax)

    % countCells counts the number of cells on the lattice
    
    nCells = 0;
    for i = iMin-1:iMax+1
       for j = jMin-1:jMax+1
           nCells = nCells + stateLattice(i,j);
       end
    end
    
end

%% Bias Toward EMF

function indices = biasTowardEMF(currentInds,cellArray,magneticField,minAngle,maxAngle,i,j)

    % Change the list of the indices for the possible locations a cell
    % might bud into to bias growth according to the direction of the
    % magnetic field.
    
    indices = currentInds;
    
    if ~isequal(magneticField,[0 0])
        
        % Create array containing coordinates of the lattice sites
        % indicated by the list of indices. Lattice sites which are not
        % referenced 
        [numRows,~] = size(cellArray);
        possibleCells = zeros(numRows,2);
        y = 1;
        for x = 1:numRows
            if currentInds(y) == x
                possibleCells(x,1) = cellArray(x,1);
                possibleCells(x,2) = cellArray(x,2);
                if y<length(currentInds)
                    y = y+1;
                end
            else
               possibleCells(x,1) = 0;
               possibleCells(x,2) = 0;
            end
        end
        % Create an array of the angles of each 
        possibleEMFAngleArray = createEMFAngleArray(possibleCells,magneticField,i,j);
        
        for x = 1:length(possibleEMFAngleArray)
            if ~((possibleEMFAngleArray(x) > minAngle)&&(possibleEMFAngleArray(x) < maxAngle))
                possibleEMFAngleArray(x) = 0;
            end
        end
        indices = find(possibleEMFAngleArray ~= 0);
    end
end


%% Timestep or Cell Count for While Condition

function whileCondition = endCondition(FINAL_CELL_COUNT,FINAL_TIMESTEP,TIME_OR_COUNT,cellsno,time)
    if strcmp(TIME_OR_COUNT,'time')
        whileCondition = time<FINAL_TIMESTEP;
    elseif strcmp(TIME_OR_COUNT,'count')
        whileCondition = cellsno<FINAL_CELL_COUNT;
    end
end
