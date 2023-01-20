function runSimulation

% Runs the simulation once and gives measurements

%% Parameters

SEED_TYPE = 'centre';   % centre, leftSide, bothSides

% END CONDITION SETTINGS
TIME_OR_COUNT = 'time';   % string to set end condition. 'count' to stop at FINAL_CELL_COUNT cells, 'time' to stop after FINAL_TIMESTEP timesteps, 'both' to stop after whichever it reaches first.
FINAL_CELL_COUNT = 10000;   % number of cells that colony should reach in order fo the program to stop.
FINAL_TIMESTEP = 6977;   % number of timesteps the program will go to before stopping.
endConditions = struct('endCondition',TIME_OR_COUNT,'finalCellNo',FINAL_CELL_COUNT,'finalTime',FINAL_TIMESTEP);

% BUDDING PATTERN INFO
AXIAL_FRAC = 0.6;   % overall fraction of cells budding normally which bud axially: 0 = average diploid colony, 0.6 = average haploid colony.
UNIPOLAR_ON = false;   % bool set whether the colony will switch to filamentous growth in low nutrient conditions.

% MUTATION SETTINGS
% MUTATION_PROB = 0.0004104;
MUTATION_PROB = 0.0008208;
GENOTYPE_NUM = 1;

% DISPLAY SETTINGS
DISPLAY_TYPE = 'show';   % show, save, or none

% DIFFUSION SETTINGS
NUTRIENT_STEPS = 10;   % number of steps that a nutrient packet takes during its random walk, correlates to the diffusion coefficient of the media. 10 correlates to 1.5% agar weight

% LATTICE SIZE SETTINGS (auto set according to end conditions)
xylatticeSize = getLatticeSize(FINAL_CELL_COUNT,FINAL_TIMESTEP,TIME_OR_COUNT);
latticeSize = struct('x',xylatticeSize,'y',xylatticeSize);

% NUTRIENT SETTINGS
NUTRS_FOR_BUDDING = 1;   % number of nutrients necessary for a cell to bud.
startNutrs = 20;
AGAR_WEIGHT = 0.3;

INIT_NUTRS_LAT = zeros(latticeSize.('y'),latticeSize.('x')) + startNutrs;

% lattice = repmat(struct('cell',0,'angle',0,'emfAngle',0),latticeSize.('y'),latticeSize.('x'));
lattice = repmat(struct('cell',0,'angle',0),latticeSize.('y'),latticeSize.('x'));
nutrientLattice = INIT_NUTRS_LAT;
stateLattice = zeros(latticeSize.('y'),latticeSize.('x'));
mutationLattice = zeros(latticeSize.('y'),latticeSize.('x'));
lattices = struct('cells',lattice,'nutrients',nutrientLattice,'states',stateLattice,'mutants',mutationLattice);

% Set properties of seed cells.
initialBudScars = [0 1;0 -1;1 1;1 0;1 -1;-1 1;-1 0;-1 -1];

if strcmp(SEED_TYPE,'leftSide')
    for k = 1:latticeSize.('y')
        i1 = randi([1,8]);
        initialScar1 = [initialBudScars(i1,1) initialBudScars(i1,2)];

        seedCell1 = YeastCell(initialScar1,AXIAL_FRAC,NUTRS_FOR_BUDDING,GENOTYPE_NUM,MUTATION_PROB);
        lattices.('cells')(k,1).('cell') = seedCell1;
        lattices.('states')(k,1) = 1;
        lattices.('mutants')(k,1) = GENOTYPE_NUM;
    end
elseif strcmp(SEED_TYPE,'bothSides') 
    for k = 1:latticeSize.('y')
        i1 = randi([1,8]);
        initialScar1 = [initialBudScars(i1,1) initialBudScars(i1,2)];

        i2 = randi([1,8]);
        initialScar2 = [initialBudScars(i2,1) initialBudScars(i2,2)];

        seedCell1 = YeastCell(initialScar1,AXIAL_FRAC,NUTRS_FOR_BUDDING,GENOTYPE_NUM,MUTATION_PROB);
        lattices.('cells')(k,1).('cell') = seedCell1;
        lattices.('states')(k,1) = 1;
        lattices.('mutants')(k,1) = GENOTYPE_NUM;
    
        seedCell2 = YeastCell(initialScar2,AXIAL_FRAC,NUTRS_FOR_BUDDING,GENOTYPE_NUM,MUTATION_PROB);
        lattices.('cells')(k,latticeSize.('x')).('cell') = seedCell2;
        lattices.('states')(k,latticeSize.('x')) = 1;
        lattices.('mutants')(k,latticeSize.('x')) = GENOTYPE_NUM;
    end
elseif strcmp(SEED_TYPE,'centre')
    xcentre = latticeSize.('x')/2;
    ycentre = latticeSize.('y')/2;
    for k = -2:2
        i = randi([1,8]);
        initialScar = [initialBudScars(i,1) initialBudScars(i,2)];

        seedCell1 = YeastCell(initialScar,AXIAL_FRAC,NUTRS_FOR_BUDDING,GENOTYPE_NUM,MUTATION_PROB);
        lattices.('cells')(ycentre+k,xcentre).('cell') = seedCell1;
        lattices.('states')(ycentre+k,xcentre) = 1;
        lattices.('mutants')(ycentre+k,xcentre) = GENOTYPE_NUM;
    end
end

folderName = 'images';
displaySettings = struct('displayType',DISPLAY_TYPE,'folderName',folderName);

% MAGNETIC FIELD SETTINGS
MF_STRENGTH = 0; % probability the magnetic field bias will be applied, used to control the strength of the magnetic field. 1 for full strength field, 0 for no field.
MAGNETIC_FIELD = [0 0];   % vector for the direction of the magnetic field, according to x-y coordinates, not row-column, if there is no magnetic field set to [0 0].
MIN_ANGLE = 30;   % the minimum angle from the magnetic field of the range of angles the MF biases budding towards.
MAX_ANGLE = 150;   % the maximum angle from the magnetic field of the range of angles the MF biases budding towards.
magnet = MagneticField(MAGNETIC_FIELD,MF_STRENGTH,MIN_ANGLE,MAX_ANGLE);

%% Run
disp('Axial Frac:')
disp(num2str(AXIAL_FRAC))

disp('Lattice Size:')
disp(num2str([latticeSize.('x') latticeSize.('y')]))

[cellCount,mutantCount,subColonyCount,cellsnoList,time,lattices] = colonySimulation(NUTRIENT_STEPS,endConditions,latticeSize,lattices,AGAR_WEIGHT,magnet,UNIPOLAR_ON,displaySettings);
cellsno = cellsnoList(1);
wtCellsno = cellsnoList(2);
mutantCellsno = cellsnoList(3);

%% Display Results

disp('time')
disp(time)

disp('cell count')
disp(cellsno)

disp('WT cell count')
disp(wtCellsno)

disp('mutant cell count')
disp(mutantCellsno)

disp('mutant subcolony count')
disp(subColonyCount(end))

disp('mutant cell/subcolony counts')
disp([mutantCount subColonyCount])

% perimetre = findPerimetre(lattices.('states'),false);
% disp('perimetre')
% disp(perimetre)
% 
% area = bwarea(lattices.('states'));
% disp('area')
% disp(area)
% 
% convexHull = bwconvhull(lattices.('states'));
% 
% % Convexity gives a measurement of how irregular the boundary. Maximum
% % value of 1, the lower the number the less convex the colony is,
% % therefore the more irregular the boundary is.
% convexPerimetre = findPerimetre(convexHull,false);
% disp('convexity')
% disp(convexPerimetre/perimetre)
% 
% % Solidity gives a measurement of the density of the colony, both in
% % terms of how many holes the colony contains as well as how irregular
% % the boundary is. Maximum value of 1, the smaller the value the less
% % solid the colony is.
% convexArea = bwarea(convexHull);
% disp('solidity')
% disp(area/convexArea)
% 
% % Compactness is the ratio of the area of the colony with the area of a
% % circle with the same perimeter. Gives a measurement of much the
% % colony varies from a circle, both in the overall shape and in the
% % irregularity of the boundary.
% compactness = (4*pi*area)/(perimetre^2);
% disp('compactness')
% disp(compactness)
% 
% % Roundness is the ration of the area of the colony to the area of a
% % circle with the same convex perimeter. Gives a measurement of how
% % circular the colony is with less sensitivity to the irregularities in
% % the foundary, focusing more on the overall shape of the colony.
% % Maximum value of 1, the lower the value, the less circular the colony
% % is.
% roundness = (4*pi*area)/(convexPerimetre^2);
% disp('roundness')
% disp(roundness)
% 
% % Elongation is the ratio of the width and height of the bounding box
% % of the object. Maximum value of 1. As the value decreases, the colony
% % is more elongated.
% % The BoundingBox property given by regionprops does not
% % rotate according to the orientation of the colony, so the property
% % 'Orientation' is used to adjust for that.
% measurements = regionprops(lattices.('states'), 'Orientation');
% % Compute angle from measurements.Orientation
% angle = measurements.Orientation;
% % Rotate image into upright position.
% uprightImage = imrotate(lattices.('states'), angle);
% % Find the bounding box of the rotated image.
% boundingBox = regionprops(uprightImage,'BoundingBox');
% for k = 1 : length(boundingBox)
%      BB = boundingBox(k).BoundingBox;
% end
% % Calculate elongation.
% widthBoundingBox = min([BB(3) BB(4)]);
% lengthBoundingBox = max([BB(3) BB(4)]);
% disp('elongation')
% disp(lengthBoundingBox/widthBoundingBox)
% 
% % Holes measures the number of holes in the colony, giving a
% % measurement of how dense the colony is.
% eul = bweuler(lattices.('states'),8);
% disp('number of holes')
% disp(1-eul)
% 
% % The mean normalised radial length is the mean of all lengths between the
% % centroid of each pixel along the boundary. These radial lengths are
% % normalised by dividing by the maximum radial length. The mean and
% % standard deviation can be used to measure boundary fluctuations.
% [~,~,meanNormRadialLength,normRadialLengthStandDev] = findRadialLengths(lattices.('states'));
% disp('mean normalised radial length')
% disp(meanNormRadialLength)
% disp('standard deviation of normalised radial lengths')
% disp(normRadialLengthStandDev)

function latticeSize = getLatticeSize(finalCellCount,finalTimestep,timeOrCount)
    if strcmp(timeOrCount,'count')
        % Determine the lattice size necessary for the set FINAL_CELL_COUNT
        % Determine the largest possible diagonal according to colony with
        % area = 2*finalCellCount and elongation 2
        maxDiag = 2*(sqrt((2*finalCellCount)/(0.5*pi))+50);
    end

    if strcmp(timeOrCount,'time')
        % Determine the lattice size necessary for the set FINAL_TIMESTEP
        growthRate = 20000/320000;   % to match the lattice size of final cell count = 20 000 when final timestep = 320 000
        % Determine the largest possible diagonal according to colony with area
        % 2*growthRate*finalTimestep and elongation 2
        maxDiag = 2*(sqrt((2*growthRate*finalTimestep)/(0.5*pi))+50);
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
end

end


