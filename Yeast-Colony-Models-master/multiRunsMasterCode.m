function multiRunsMasterCode(NUMBER_OF_RUNS,DIFFUSION_STEPS,FINAL_CELL_COUNT,FINAL_TIMESTEP,TIME_OR_COUNT,NUTRS_FOR_BUDDING,MIN_ANGLE,MAX_ANGLE,MAGNETIC_FIELD,START_NUTRS,MF_STRENGTH,AXIAL_FRAC,UNIPOLAR_ON)

% Runs the simulation repeatedly for a given set of conditions.
validateattributes(NUMBER_OF_RUNS,{'numeric'},{'>',0,'size',[1,1]})

%% Preallocate arrays for all measurements.
finalCellCounts = zeros(NUMBER_OF_RUNS,1);
timeList = zeros(NUMBER_OF_RUNS,1);
perimeterList = zeros(NUMBER_OF_RUNS,1);
areaList = zeros(NUMBER_OF_RUNS,1);
convexityList = zeros(NUMBER_OF_RUNS,1);
solidityList = zeros(NUMBER_OF_RUNS,1);
densityList = zeros(NUMBER_OF_RUNS,1);
compactnessList = zeros(NUMBER_OF_RUNS,1);
roundnessList = zeros(NUMBER_OF_RUNS,1);
elongationList = zeros(NUMBER_OF_RUNS,1);
holesList = zeros(NUMBER_OF_RUNS,1);
meanLengthsList = zeros(NUMBER_OF_RUNS,1);
stdDevLengthList = zeros(NUMBER_OF_RUNS,1);

%% Run
parfor i = 1:NUMBER_OF_RUNS
    [~, cellno, time, stateLattice] = colonySimulation(DIFFUSION_STEPS,FINAL_CELL_COUNT,FINAL_TIMESTEP,TIME_OR_COUNT,START_NUTRS,NUTRS_FOR_BUDDING,AXIAL_FRAC,MAGNETIC_FIELD,MF_STRENGTH,MIN_ANGLE,MAX_ANGLE,UNIPOLAR_ON,false,0,false);

%% Save Results

    finalCellCounts(i) = cellno;
    
    timeList(i) = time;
    
    perimeter = findPerimetre(stateLattice,false);
    perimeterList(i) = perimeter;
    
    area = bwarea(stateLattice);
    areaList(i) = area;

    convexHull = bwconvhull(stateLattice);

    % Convexity gives a measurement of how irregular the boundary. Maximum
    % value of 1, the lower the number the less convex the colony is,
    % therefore the more irregular the boundary is.
    convexPerimetre = findPerimetre(convexHull,false);
    convexity = convexPerimetre/perimeter;
    convexityList(i) = convexity;
    
    % Solidity gives a measurement of the density of the colony, both in
    % terms of how many holes the colony contains as well as how irregular
    % the boundary is. Maximum value of 1, the smaller the value the less
    % solid the colony is.
    convexArea = bwarea(convexHull);
    solidity = area/convexArea;
    solidityList(i) = solidity;

    % Compactness is the ratio of the area of the colony with the area of a
    % circle with the same perimeter. Gives a measurement of much the
    % colony varies from a circle, both in the overall shape and in the
    % irregularity of the boundary.
    compactness = (4*pi*area)/(perimeter^2);
    compactnessList(i) = compactness;

    % Roundness is the ration of the area of the colony to the area of a
    % circle with the same convex perimeter. Gives a measurement of how
    % circular the colony is with less sensitivity to the irregularities in
    % the foundary, focusing more on the overall shape of the colony.
    % Maximum value of 1, the lower the value, the less circular the colony
    % is.
    roundness = (4*pi*area)/(convexPerimetre^2);
    roundnessList(i) = roundness;

    % Elongation is the ratio of the width and height of the bounding box
    % of the object. Maximum value of 1. As the value increases, the colony
    % is more elongated.
    % The BoundingBox property given by regionprops does not
    % rotate according to the orientation of the colony, so the property
    % 'Orientation' is used to adjust for that.
    measurements = regionprops(stateLattice, 'Orientation');
    % Compute angle from measurements.Orientation.
    angle = measurements.Orientation;
    % Rotate image.
    uprightImage = imrotate(stateLattice, angle);
    % Find bounding box of rotated image.
    boundingBox = regionprops(uprightImage,'BoundingBox');
    BB = 0;
    for k = 1 : length(boundingBox)
         BB = boundingBox(k).BoundingBox;
    end
    % Calculate elongation.
    widthBoundingBox = min([BB(3) BB(4)]);
    lengthBoundingBox = max([BB(3) BB(4)]);
    elongation = (lengthBoundingBox/widthBoundingBox);
    elongationList(i) = elongation;
    
    % Holes measures the number of holes in the colony, giving a
    % measurement of how dense the colony is.
    eul = bweuler(stateLattice,8);
    holes = 1-eul;
    holesList(i) = holes;

    % The mean normalised radial length is the mean of all lengths between the
    % centroid of each pixel along the boundary. These radial lengths are
    % normalised by dividing by the maximum radial length. The mean and
    % standard deviation can be used to measure boundary fluctuations.
    [~,~,meanNormRadialLength,normRadialLengthStdDev] = findRadialLengths(stateLattice);
    meanLengthsList(i) = meanNormRadialLength;
    stdDevLengthList(i) = normRadialLengthStdDev;
end

% Build table containing all measurements. Output as an .xlsx file into
% folder titled "matlab_output"
t1 = table(finalCellCounts,timeList,perimeterList,areaList,convexityList,solidityList,densityList,compactnessList,roundnessList,elongationList,holesList,meanLengthsList,stdDevLengthList);

fileName = buildFileNameStr(MAGNETIC_FIELD,DIFFUSION_STEPS,START_NUTRS,NUTRS_FOR_BUDDING,MF_STRENGTH,AXIAL_FRAC,UNIPOLAR_ON,TIME_OR_COUNT);

fileName1 = strcat('matlab_output/',fileName);

writetable(t1,fileName1);

%% Name File
function fileName = buildFileNameStr(MF,nSteps,concentration,budNutrs,strength,ploidy,unipolar_on,end_condition)
    % Creates string for the name of the file which will be output based on
    % given parameters.
    MFstring1 = num2str(MF(1));
    MFstring2 = num2str(MF(2));
    MFString = strcat('(',MFstring1,'_',MFstring2,')');
    
    if concentration > 16*budNutrs
        nutrientString = 'rich';
    elseif concentration > 5*budNutrs
        nutrientString = 'medium';
    else
        nutrientString = 'low';
    end

    if strength > 1
        strengthString = 'extraStrong';
    elseif strength == 1
        strengthString = 'strong';
    elseif strength > 0
        strengthString = 'weak';
    else
        strengthString = 'no';
    end

    if ploidy == 0.6
        ploidyString = 'haploid';
    elseif ploidy == 0
        ploidyString = 'diploid';
    else
        ploidyString1 = 'ploidynum_';
        ploidyString2 = num2str(ploidy);
        ploidyString = strcat(ploidyString1, ploidyString2);
    end
    
    if nSteps == 0
        diffusionString = 'no_diffusion';
    else
        diffusionString = 'with_diffusion';
    end
    
    if unipolar_on
        unipolarString = 'unipolar';
    else
        unipolarString = 'not_unipolar';
    end
    
    fileName = strcat(end_condition,'_',unipolarString,'_',ploidyString,'_',nutrientString,'_',strengthString,'MF_',MFString,'_',diffusionString,'_file.xlsx');
end

end

