function masterCode

% Runs the simulation once and gives measurements

%% Parameters

nSteps = 10;   % number of steps that a nutrient packet takes during its random walk, correlates to the diffusion coefficient of the media.
START_NUTRS = 2;   % number of nutrient packets at each lattice site at the beginning.
NUTRS_FOR_BUDDING = 1;   % number of nutrients necessary for a cell to bud.

AXIAL_FRAC = 0.6;   % overall fraction of cells budding normally which bud axially: 0 = average diploid colony, 0.6 = average haploid colony.
UNIPOLAR_ON = false;   % logical value sets whether the colony will switch to filamentous growth in low nutrient conditions.

MF_STRENGTH = 1; % probability the magnetic field bias will be applied, used to control the strength of the magnetic field. 1 for strong field, 0 for no field, >1 for full strength field.
MAGNETIC_FIELD = [1 1];   % vector for the direction of the magnetic field, according to x-y coordinates, not row-column, if there is no magnetic field set to [0 0].
MIN_ANGLE = 30;   % the minimum angle from the magnetic field of the range of angles the MF biases budding towards.
MAX_ANGLE = 150;   % the maximum angle from the magnetic field of the range of angles the MF biases budding towards.

TIME_OR_COUNT = 'time';   % string or char to set end condition. 'count' to stop at FINAL_CELL_COUNT cells, 'time' to stop after FINAL_TIMESTEP timesteps.
FINAL_CELL_COUNT = 10000;   % number of cells that colony should reach in order fo the program to stop.
FINAL_TIMESTEP = 320000;   % number of timesteps the program will go to before stopping.

MUTATION_ON = true;   % logical value sets whether mutations will occur.
MUTATION_PROB = 0;   % set probability of mutation occuring during during budding.

DISPLAY_IMAGE = true;   % logical value sets whether the colony will be displayed graphically as it grows throughout the simulation.

%% Run
[cellno, time, stateLattice] = colonySimulation(nSteps,FINAL_CELL_COUNT,FINAL_TIMESTEP,TIME_OR_COUNT,START_NUTRS,NUTRS_FOR_BUDDING,AXIAL_FRAC,MAGNETIC_FIELD,MF_STRENGTH,MIN_ANGLE,MAX_ANGLE,UNIPOLAR_ON,MUTATION_ON,MUTATION_PROB,DISPLAY_IMAGE);

%% Display Results

disp('time')
disp(time)

disp('final cell count')
disp(cellno)

perimetre = findPerimetre(stateLattice,false);
disp('perimetre')
disp(perimetre)

area = bwarea(stateLattice);
disp('area')
disp(area)

convexHull = bwconvhull(stateLattice);

% Convexity gives a measurement of how irregular the boundary. Maximum
% value of 1, the lower the number the less convex the colony is,
% therefore the more irregular the boundary is.
convexPerimetre = findPerimetre(convexHull,false);
disp('convexity')
disp(convexPerimetre/perimetre)

% Solidity gives a measurement of the density of the colony, both in
% terms of how many holes the colony contains as well as how irregular
% the boundary is. Maximum value of 1, the smaller the value the less
% solid the colony is.
convexArea = bwarea(convexHull);
disp('solidity')
disp(area/convexArea)

% Compactness is the ratio of the area of the colony with the area of a
% circle with the same perimeter. Gives a measurement of much the
% colony varies from a circle, both in the overall shape and in the
% irregularity of the boundary.
compactness = (4*pi*area)/(perimetre^2);
disp('compactness')
disp(compactness)

% Roundness is the ration of the area of the colony to the area of a
% circle with the same convex perimeter. Gives a measurement of how
% circular the colony is with less sensitivity to the irregularities in
% the foundary, focusing more on the overall shape of the colony.
% Maximum value of 1, the lower the value, the less circular the colony
% is.
roundness = (4*pi*area)/(convexPerimetre^2);
disp('roundness')
disp(roundness)

% Elongation is the ratio of the width and height of the bounding box
% of the object. Maximum value of 1. As the value decreases, the colony
% is more elongated.
% The BoundingBox property given by regionprops does not
% rotate according to the orientation of the colony, so the property
% 'Orientation' is used to adjust for that.
measurements = regionprops(stateLattice, 'Orientation');
% Compute angle from measurements.Orientation
angle = measurements.Orientation;
% Rotate image into upright position.
uprightImage = imrotate(stateLattice, angle);
% Find the bounding box of the rotated image.
boundingBox = regionprops(uprightImage,'BoundingBox');
for k = 1 : length(boundingBox)
     BB = boundingBox(k).BoundingBox;
end
% Calculate elongation.
widthBoundingBox = min([BB(3) BB(4)]);
lengthBoundingBox = max([BB(3) BB(4)]);
disp('elongation')
disp(lengthBoundingBox/widthBoundingBox)

% Holes measures the number of holes in the colony, giving a
% measurement of how dense the colony is.
eul = bweuler(stateLattice,8);
disp('number of holes')
disp(1-eul)

% The mean normalised radial length is the mean of all lengths between the
% centroid of each pixel along the boundary. These radial lengths are
% normalised by dividing by the maximum radial length. The mean and
% standard deviation can be used to measure boundary fluctuations.
[~,~,meanNormRadialLength,normRadialLengthStandDev] = findRadialLengths(stateLattice);
disp('mean normalised radial length')
disp(meanNormRadialLength)
disp('standard deviation of normalised radial lengths')
disp(normRadialLengthStandDev)




