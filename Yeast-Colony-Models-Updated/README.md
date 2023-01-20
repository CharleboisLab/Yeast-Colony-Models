# Yeast-Colony-Models-Updated

Updated version of Yeast-Colony-Models. A lattice-based simulation of yeast colony growth under a number of different conditions including different budding patterns, varying nutrient concentrations, neutral mutations, and the influences of static magnetic fields. Code was written by Rebekah Hall and questions regarding the code should be directed to rebekah_hall@sfu.ca.


## Dependencies 

- MATLAB code requires the [Image Processing Toolbox](https://www.mathworks.com/products/image.html)
- Some MATLAB code depends on the [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html)
  - neutralMutationMultiRun.m

## Description of files

#### MATLAB files:

| filename                    | description                                                  |
| --------------------------- | ------------------------------------------------------------ |
| colonySimulation.m          | Run a single simulation of yeast colony growth.              |
| runSimulation.m             | Runs colonySimulation.m a single time with the given parameters.|
| neutralMutationMultiRun.m       | Runs colonySimulation.m a set number of times under neutral mutation conditions and outputs results into an .xlsx files. |
| YeastCell.m         | Class describing the cell state and nutrient concentration at a given lattice point. |
| findPerimetre.m             | Calculates perimeter of a colony.                            |
| findRadialLengths.m         | Calculates lengths from the centre of the colony to each of the boundary points. |

## Usage

### Overview of simulation parameters

| parameter          | description                                                  |
| ------------------ | ------------------------------------------------------------ |
| SEED_TYPE          | Sets the location of beginning seed cells. `'leftSide'` seeds wildtype cells along the leftmost column. `'bothSides'` seeds wildtype cells along the leftmost and rightmost columns. `'centre'` seeds wildtype cells a column of five wildtype cells in the centre of the grid |
| TIME_OR_COUNT      | Set to 'time' for simulation to end after a set number of timesteps. Set to 'count' for simulation to end when the colony reaches a set number of cells |
| FINAL_CELL_COUNT   | Number of cells for the colony to reach in order to stop the simulation, if TIME_OR_COUNT = 'count' |
| FINAL_TIMESTEP     | Number of timesteps the simulation will run for if TIME_OR_COUNT = 'time' |
| AXIAL_FRAC         | Fraction of cells which will bud axially.                    |
| UNIPOLAR_ON        | true for pseudohyphal growth, false for no pseudohyphal growth condition |
| MUTATION_PROB      | Probability of a cell mutating if MUTATION_ON = true         |
| DISPLAY_TYPE     | Sets whether colony visualisations at every 100 time steps will be shown (`'show'`), saved as .jpg files (`'save'`), or not visualised |
| NUTRIENT_STEPS            | Number of steps a single nutrient packet will take on its random walk. A higher value corresponds to a more diffusive media |
| NUTRS_FOR_BUDDING  | Number of nutrient packets a cell must consume in order to bud. |
| startNutrs        | Number of nutrient packets at every lattice point at the beginning of the simulation. Represents initial nutrient concentration |
| AGAR_WEIGHT       | Agar concentration |
| MF_STRENGTH        | Fraction of the time that the magnetic field bias is applied. Set to 0 for no magnetic field. |
| MAGNETIC_FIELD     | Vector setting the direction of the magnetic field.          |
| MIN_ANGLE          | The lower end of the range of angles into which the cell will be biased toward budding. |
| MAX_ANGLE          | The upper end of the range of angles into which the cell will be biased toward budding. |

### Code authors
Code was written by Rebekah Hall under the supervision of Dr. Daniel Charlebois. 





