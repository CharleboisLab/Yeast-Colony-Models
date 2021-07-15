function multiRunForBuddingAngles
%% Parameters

nSteps = 10;   % number of steps that a nutrient packet takes during its random walk, correlates to the diffusion coefficient of the media.
NUTRS_FOR_BUDDING = 1;   % number of nutrients necessary for a cell to bud.

UNIPOLAR_ON = false;   % logical value sets whether the colony will switch to filamentous growth in low nutrient conditions.

MIN_ANGLE = 30;   % the minimum angle from the magnetic field of the range of angles the MF biases budding towards.
MAX_ANGLE = 150;   % the maximum angle from the magnetic field of the range of angles the MF biases budding towards.

TIME_OR_COUNT = 'time';   % string or char to set end condition. 'count' to stop at FINAL_CELL_COUNT cells, 'time' to stop after FINAL_TIMESTEP timesteps.
FINAL_CELL_COUNT = 10000;   % number of cells that colony should reach in order fo the program to stop.
FINAL_TIMESTEP = 320000;   % number of timesteps the program will go to before stopping.

MUTATION_ON = false;   % logical value sets whether mutations will occur.
MUTATION_PROB = 0;   % set probability of mutation occuring during during budding.

DISPLAY_IMAGE = false;   % logical value sets whether the colony will be displayed graphically as it grows throughout the simulation.

ploidyArray = [0.6 0];   %% 0.6 = haploid, 0 = diploid
concentrationArray = [2 20];   %% [0,5] = low nutrients, [6,16] = some nutrients, >16 = rich nutrients
strengthArray = [0 0.5 1];   %% 0 = no MF, 0.5 = weak MF, 1 = strong MF, >1 = extra strong MF
directionArray = [1 1;1 0;1 -1;0 1;0 -1;-1 1;-1 0;-1 -1];

%% Run       
for m = 1:length(ploidyArray)
    ploidy = ploidyArray(m);
    for k = 1:length(strengthArray)
        strength = strengthArray(k);
        for j=1:length(concentrationArray)
            concentration = concentrationArray(j);
            for i=1:length(directionArray)
                x1 = directionArray(i,1);
                x2 = directionArray(i,2);
                MF = [x1 x2];
                [buddingAngles, ~, ~, ~] = colonySimulation(nSteps,FINAL_CELL_COUNT,FINAL_TIMESTEP,TIME_OR_COUNT,concentration,NUTRS_FOR_BUDDING,ploidy,MF,strength,MIN_ANGLE,MAX_ANGLE,UNIPOLAR_ON,MUTATION_ON,MUTATION_PROB,DISPLAY_IMAGE);
                %% Create Budding Angles Table
                fileName = buildFileNameStr(MF,nSteps,concentration,NUTRS_FOR_BUDDING,strength,ploidy,UNIPOLAR_ON,TIME_OR_COUNT);
                fileName1 = strcat('matlab_output/budding_angles/buddingAngles_',fileName);
                t1 = table(buddingAngles);
                writetable(t1,fileName1);
            end
        end
    end
end
                    
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