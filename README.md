# Yeast-Colony-Models

A lattice-based simulation of yeast colony growth under a number of different conditions including different budding patterns, varying nutrient concentrations, neutral mutations, and the influences of static magnetic fields. 

## Dependencies 

- MATLAB code requires the [Image Processing Toolbox](https://www.mathworks.com/products/image.html)

- Some MATLAB code depends on the [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html)
  - multiRunsMasterCode.m
  - multiConditionsMaster

- Some Python code depends on the 'pandas' library:
  - combine_data.py
  - format_data.py
  - split_violin_plots.py
- Some Python code depends on the 'pathlib' library:
  - combine_data.py
  - format_data.py
  - split_violin_plots.py
- Some Python code depends on the 'os' library:
  - split_violin_plots.py
- Some Python code depends on the 'plotly.graph_objects' library:
  - split_violin_plots.py

## Description of files

#### MATLAB files:

| filename                    | description                                                  |
| --------------------------- | ------------------------------------------------------------ |
| colonySimulation.m          | Run a single simulation of yeast colony growth.              |
| masterCode.m                | Runs colonySimulation.m a single time with the given parameters. |
| multiRunsMasterCode.m       | Runs colonySimulation.m a set number of times and outputs results into an .xlsx file in matlab_output folder. |
| multiConditionsMasterCode.m | Runs multiRunsMasterCode.m multiple times for a list of different parameters. |
| CellWithNutrients.m         | Class describing the cell state and nutrient concentration at a given lattice point. |
| findPerimetre.m             | Calculates perimeter of a colony.                            |
| findRadialLengths.m         | Calculates lengths from the centre of the colony to each of the boundary points. |

#### Python files:

| filename              | description                                                  |
| --------------------- | ------------------------------------------------------------ |
| format_data.py        | Formats and collects output results into matlab_output folder into .csv files for each magnetic field direction and ploidy saved in python_output folder. |
| combine_data.py       | Collects the output files from format_data.py found in python_output folder into .csv files for each ploidy saved in python_output folder. |
| split_violin_plots.py | Dreates violin plots from .csv files output by combine_data.py. Saves plots in figures_output folder. |

## Usage

### Overview of simulation parameters

| parameter         | description                                                  |
| ----------------- | ------------------------------------------------------------ |
| nSteps            | Number of steps a single nutrient packet will take on its random walk. A higher nSteps corresponds to a more diffusive media |
| START_NUTRS       | Number of nutrient packets at every lattice point at the beginning of the simulation. Represents initial nutrient concentration. |
| NUTRS_FOR_BUDDING | Number of nutrient packets a cell must consume in order to bud. |
| AXIAL_FRAC        | Fraction of cells which will bud axially.                    |
| UNIPOLAR_ON       | true for pseudohyphal growth, false for no pseudohyphal growth condition |
| MF_STRENGTH       | Fraction of the time that the magnetic field bias is applied. Set to 0 for no magnetic field. |
| MAGNETIC_FIELD    | Vector setting the direction of the magnetic field.          |
| MIN_ANGLE         | The lower end of the range of angles into which the cell will be biased toward budding. |
| MAX_ANGLE         | The upper end of the range of angles into which the cell will be biased toward budding. |
| TIME_OR_COUNT     | Set to 'time' for simulation to end after a set number of timesteps. Set to 'count' for simulation to end when the colony reaches a set number of cells |
| FINAL_CELL_COUNT  | Number of cells for the colony to reach in order to stop the simulation, if TIME_OR_COUNT = 'count' |
| FINAL_TIMESTEP    | Number of timesteps the simulation will run for if TIME_OR_COUNT = 'time' |
| MUTATION_ON       | true for cells to mutate, false for cells not to mutate      |
| MUTATION_PROB     | Probability of a cell mutating if MUTATION_ON = true         |
| DISPLAY_IMAGE     | true for a visualisation of the colony to be shown, false otherwise |

### Running simulations

#### Running a single simulation with a visual

Set parameters in masterCode.m:

```matlab
function masterCode

%% Parameters

nSteps = 10;  % controls diffusion
START_NUTRS = 20;   % initial nutrient concentration
NUTRS_FOR_BUDDING = 1;   % number of nutrients necessary for a cell to bud.

AXIAL_FRAC = 0.6;   % 0.6 = average haploid colony.
UNIPOLAR_ON = false;   % no pseudohyphal growth.

MF_STRENGTH = 0; % no magnetic field
MAGNETIC_FIELD = [0 0];
MIN_ANGLE = 30;
MAX_ANGLE = 150;

TIME_OR_COUNT = 'time';   % simulation will run for a set number of timesteps
FINAL_CELL_COUNT = 0;
FINAL_TIMESTEP = 320000;   % number of timesteps the program will go to before stopping.

MUTATION_ON = false;   % no mutation
MUTATION_PROB = 0;

DISPLAY_IMAGE = true;   % a visualisation of the colony will be shown
```

Run masterCode.m in MATLAB Command Window:

```matlab
>> masterCode
```

#### Running multiple simulations for a single set of parameters

Either run multiRunsMasterCode.m in MATLAB Command Window directly:

​	The arguments are:  NUMBER_OF_RUNS, DIFFUSION_STEPS, FINAL_CELL_COUNT, FINAL_TIMESTEP, TIME_OR_COUNT, NUTRS_FOR_BUDDING, 				    	MIN_ANGLE, MAX_ANGLE, MAGNETIC_FIELD, START_NUTRS, MF_STRENGTH, AXIAL_FRAC, UNIPOLAR_ON

Example:
```matlab
>> multiRunsMasterCode(200,10,0,320000,'time',1,30,150,[0 0],20,0,0.6,false)
```

Or use multiConditionsMasterCode.m with only one parameter in each array:

```matlab
function multiConditionsMasterCode

FINAL_CELL_COUNT = 0;
FINAL_TIMESTEP = 320000;
NUTRS_FOR_BUDDING = 1;
MIN_ANGLE = 30;
MAX_ANGLE = 150;
NUMBER_OF_RUNS = 200;

endArray = {'time'};   %% sets TIME_OR_COUNT
unipolarArray = [false];   %% sets UNIPOLAR_ON
ploidyArray = [0.6];   %% sets AXIAL_FRAC
diffusionArray = [10];   % sets nSteps
concentrationArray = [20];   %% sets START_NUTRS
strengthArray = [0];   %% sets MF_STRENGTH
directionArray = [0 0];   %% sets MAGNETIC_FIELD
```

Run multiConditionsMasterCode.m in MATLAB Command Window: 

```matlab
>> multiConditionsMasterCode
```

#### Running multiple simulations for multiple parameter sets

Set parameters in multiConditionsMasterCode.m file:

```matlab
function multiConditionsMasterCode

FINAL_CELL_COUNT = 10000;   % number  of cells that colony should reach in order fo the program to stop.
FINAL_TIMESTEP = 320000;   % number of time steps the program funs for.
NUTRS_FOR_BUDDING = 1;   % number of nutrients necessary for a cell to bud.
MIN_ANGLE = 30;   % the minimum angle from the magnetic field of the range of angles the MF biases budding towards.
MAX_ANGLE = 150;   % the maximum angle from the magnetic field of the range of angles the MF biases budding towards.
NUMBER_OF_RUNS = 200;

% Loops through all necessary parameters.
endArray = {'time','count'};   %% 'time' = ends after FINAL_TIMESTEP time steps, 'count' = ends after colony reaches FINAL_CELL_COUNT cells
unipolarArray = [false];   %% true = cells switch to filamentous growth in low nutrients, false = cells never switch to filamentous growth
ploidyArray = [0.6 0];   %% 0.6 = haploid, 0 = diploid
diffusionArray = [10];   %% 0 = no diffusion, >0 = diffusion
concentrationArray = [20 2];   %% [0,5] = low nutrients, [6,16] = some nutrients, >16 = rich nutrients
strengthArray = [0 0.5 1];   %% 0 = no MF, 0.5 = weak MF, 1 = strong MF
directionArray = [1 1;1 0;1 -1;0 1;0 -1;-1 1;-1 0;-1 -1];
```

Run multiConditionsMasterCode.m: 

```matlab
>> multiConditionsMasterCode
```

#### File names

Results are output in .xlsx files with the naming convention:

​	*endCondition* _ *unipolar* _ *ploidy* _ *nutrientConc* _ *strength*MF _ *direction* _ *diffusion* _ file.xlsx

Where:

- *endCondition* = `TIME_OR_COUNT`
- *unipolar* = 
  - unipolar, if `UNIPOLAR_ON = true`
  - not_unipolar, if `UNIPOLAR_ON = false`
- *ploidy* = 
  - haploid, if `AXIAL_FRAC = 0.6`
  - diploid, if `AXIAL_FRAC = 0`
  - ploidynum_`AXIAL_FRAC`, otherwise
- *nutrientConc* =
  - rich, if `START_NUTRS > 16*NUTRS_TO_BUD`
  - medium, if `5*NUTRS_TO_BUD < START_NUTRS < 16*NUTRS_TO_BUD`
  - low, if `START_NUTRS < 5*NUTRS_TO_BUD`
- *strength* =
  - extraStrong, if `MF_STRENGTH > 1`
  - strong, if `MF_STRENGTH = 1`
  - weak, if `0< MF_STRENGTH < 1`
  - no, if `MF_STRENGTH = 0`
- *direction* = (`MAGNETIC_FIELD(1)`_`MAGNETIC_FIELD(2)`)
- *diffusion* =
  - no_diffusion, if `nSteps = 0`
  - with_diffusion, otherwise

### Collecting results into single .csv file and creating violin plots

1. Run format_data.py:

   Set correct strings according to the names of the .xlsx files (lines 11-16):

   ```python
   # Example    
       budding_pattern_list = ['not_unipolar']  # for unipolar
       ploidy_list = ['haploid', 'diploid']   # for ploidy
       diffusion_steps_list = ['with_diffusion']   # for diffusion
       conc_list = ['rich', 'low']   # for nutrientConc
       strengths_list = ['strong', 'weak', 'no']   # for strength
       time_count = 'time'   # for endCondition, not a list
   ```

   format_data.py outputs .csv files to python_output with the naming convention:

   ​	*endCondition* _ *unipolar* _ *ploidy* _ *nutrientConc* _ nutrient _ *diffusion* _ *strength* _MF.csv

2. Run combine_data.py:

   Set correct strings according to the names of the .csv files from format_data.py (lines 13-23):

   ```python
       ploidy_list = ['haploid', 'diploid']   # for ploidy
       concentration_list = ['rich', 'low']   # for nutrientConc
       strength_list = ['strong', 'weak', 'no']   # for strength
       diffusion_list = ['with_diffusion']   # for diffusion
       budding_pattern_list = ['not_unipolar']   # for unipolar
       time_count = 'time'   # for endCondition, not a list
   ```

   combine_data.py outputs .csv files to python_output with the naming convention:

   ​	*endCondition* _ all_results_ *ploidy* _ *unipolar*.csv

3. Run split_violin_plots.py:

   Set correct strings according to the name of the .csv files from combine_data.py (lines 63-66):

   ```python
       ploidy_list = ['haploid', 'diploid']   # for ploidy
       budding_list = ['not_unipolar']   # for unipolar
       diffusion_list = ['with_diffusion']   # for diffusion
       end_list = ['time', 'count']   # for endCondition
   ```

### Code authors
Code was written by Rebekah Hall under the supervision of Dr. Daniel Charlebois. 


   

   
