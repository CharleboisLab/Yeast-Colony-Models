function multiConditionsMasterCode

% Runs the simulation under all conditions for repeated runs

FINAL_CELL_COUNT = 10000;   % number  of cells that colony should reach in order fo the program to stop.
FINAL_TIMESTEP = 320000;   % number of time steps the program funs for.
NUTRS_FOR_BUDDING = 1;   % number of nutrients necessary for a cell to bud.
MIN_ANGLE = 30;   % the minimum angle from the magnetic field of the range of angles the MF biases budding towards.
MAX_ANGLE = 150;   % the maximum angle from the magnetic field of the range of angles the MF biases budding towards.
NUMBER_OF_RUNS = 100;

endArray = {'time','count'};   %% 'time' = ends after FINAL_TIMESTEP time steps, 'count' = ends after colony reaches FINAL_CELL_COUNT cells
unipolarArray = [false];   %% true = cells switch to filamentous growth in low nutrients, false = cells never switch to filamentous growth
ploidyArray = [0.6 0];   %% 0.6 = haploid, 0 = diploid
diffusionArray = [10];   %% 0 = no diffusion, >0 = diffusion
concentrationArray = [2 20];   %% [0,5] = low nutrients, [6,16] = some nutrients, >16 = rich nutrients
strengthArray = [0 0.5 1 2];   %% 0 = no MF, 0.5 = weak MF, 1 = strong MF, >1 = extra strong MF
directionArray = [1 1;1 0;1 -1;0 1;0 -1;-1 1;-1 0;-1 -1];

for o = 1:length(endArray)
    end_condition = endArray{o};
    disp('END CONDITION')
    disp(end_condition)
    
    for n = 1:length(unipolarArray)
        unipolar_on = unipolarArray(n);
        disp('UNIPOLAR')
        if unipolar_on
            disp('unipolar on')
        else
            disp('unipolar off')
        end
        
        for m = 1:length(ploidyArray)
            ploidy = ploidyArray(m);
            disp('PLOIDY')
            disp(ploidy)

            for l = 1:length(diffusionArray)
                nSteps = diffusionArray(l);
                disp('DIFFUSION')
                disp(nSteps)

                for k = 1:length(strengthArray)
                    strength = strengthArray(k);
                    disp('STRENGTH')
                    disp(strength)

                    for j=1:length(concentrationArray)
                        concentration = concentrationArray(j);
                        disp('CONCENTRATION')
                        disp(concentration)

                        for i=1:length(directionArray)
                            x1 = directionArray(i,1);
                            x2 = directionArray(i,2);
                            MF = [x1 x2];
                            disp('MF')
                            disp(MF)
                            multiRunsMasterCode(NUMBER_OF_RUNS,nSteps,FINAL_CELL_COUNT,FINAL_TIMESTEP,end_condition,NUTRS_FOR_BUDDING,MIN_ANGLE,MAX_ANGLE,MF,concentration,strength,ploidy,unipolar_on)
                        end
                    end
                end
            end
        end
    end
end

