# Drug-Resistance-Models

A lattice-based simulation of colony growth and the development of genetic and non-genetic drug resistance under the influence of varying nutrient and drug conditions with two different models of drug resistance. Code was written by Rebekah Hall and questions regarding the code should be directed to rebekah_hall@sfu.ca.


## Dependencies 

- Code written using MATLAB Version 9.10.0 (R2021a)
- Requires the [Image Processing Toolbox](https://www.mathworks.com/products/image.html)

## Description of files

#### MATLAB files:

| filename                    | description                                                  |
| --------------------------- | ------------------------------------------------------------ |
| colonySimulation.m          | Run a single simulation of colony growth.                    |
| runSimulation.m             | Runs colonySimulation.m a single time with the given parameters. |
| Cell.m         | Class describing the state and behaviour of a single cell. |
| findPerimetre.m             | Calculates perimeter of a colony.                            |
| findRadialLengths.m         | Calculates lengths from the centre of the colony to each of the boundary points. |

## Simulation Details

### Overview of mutation modes

There are two mutation modes which can be used: the SNG model and the mega-plate model. 

#### SNG model

In the SNG model, cells can exist in three different states: susceptible, non-genetically resistant, or genetically resistant. Cells can switch back and forth between susceptible and non-genetically resistant, or mutate to permanently become genetically resistant. 

Susceptible cells have no resistance to drugs and genetically resistant cells are completely resistant to drugs, while the resistance of non-genetically resistant cells can be varied by varying the probability a cell will divide or survive in the presence of a static or cidal drug respectively. 

`colonySimulation.m` will output the appearance time, establishment time, and fixation time (in timesteps) of the genetically resistant cells, a table containing the number of each type of cell sample at various timesteps (approximatedly 1 minute intervals), as well as the final number of each type of cell and total number of timesteps the simulation ran for.

The default colour map set by `runSimulation.m` sets susceptible cells to appear as white, non-genetically resistant cells to appear as orange, and genetically resistant cells to appear as blue.

#### Mega-Plate model

In the mega-plate models, cells are assigned a genotype number and phenotype modifier which determine that cell's resistance. A genotype number of 1 is the starting genotype and represents a fully susceptible cell. Based on the set mutation probability, a cell may mutate to increase its genotype number by one, representing an increase in resistance. A cell will not mutate down to a lower genotype number. The resistance level can then be modified by its phenotype number. A phenotype number of 0 means the resistance level a cell expresses will be the same as its genotype level. A cell can also have a phenotype number of -1 or 1 which will decrease or increase respectively the expressed resistance level relative to the genotype level. For example, a cell may have a genotype of 6 but given a phenotype number of -1, will express like a cell with a genotype of 5. 

The probability of a cell being unaffected by the drugs is given by the formula $\phi = \max[0, 1-(\frac{D}{10^R})^2]$, where $D$ is the number of drug packets located at the lattice point where the cell exists, and $R$ is the resistance, given by the sum of the cell's genotype number with its phenotype modifier. In the case of a static drug, $\phi$ sets the probability the cell will divide. In the case of a cidal drug, $\phi$ gives the probability the cell will live. With this formula, a cell with resistance $R$ will be fully susceptible as drug concentrations greater than or equal to $D = 10^R$.

Drugs do not diffuse in mega-plate mode.

`colonySimulation.m` will output the average phenotype level (genotype number + phenotype modifier) and average genotype level (genotype number only) of all the cells on the lattice at approximately 1 minute intervals, with exact timestep values given in a second row. It will also output the appearance, establishment, and fixation times (in timesteps) of each genotype level, as well as the final number of cells of each genotype level, the final total cell count, and total number of timesteps the simulation ran for. 

`colonySimulation.m` will represent a simulation of the mega-plate model with two side by side representations of the state of the lattice. On the left the phenotype level is displayed, while on the right the genotype level is displayed. In the default colour map set by `runSimulation.m`, a phenotype/genotype level of 1 is represented by white. Level 2 is red, level 3 orange, level 4 yellow, level 5 green, level 6 blue, and levels 7+ purple. When no cells are present, black indicated no drug is present. Light grey represents low drug concentrations, and darker grey bands represent increasing drug concentrations.

#### Note on time

Time is measured in timesteps. Each time a nutrient packet takes a random walk, the number of timesteps increases by one. Based on a random walk of 10, a single timestep is estimated to correspond to 0.0129 minutes. Using this conversion, the simulation time is converted to minutes and the representation of the state of the lattice is updated every ten minutes. In mega-plate mode the average phenotype and genotype levels are calculated and saved every ten minutes as well. Appearance, establishment, and fixation times are saved and output in the number of timesteps, as is the final simulation run time.

## Usage

### Overview of end condition settings

End condition settings are stored and input into colonySimulation.m (`endConditions` argument) as the structure array with the following fields.

| parameter          | description                                                  |
| ------------------ | ------------------------------------------------------------ |
| endCondition| A string which sets what type of end condition will be used. `'count'` for an end condition based on final number of cells; `'time'` for an end condition based on final number of timesteps; `'both'` for an end condition based on whichever the program reaches first, the set final number of timesteps or the set final cell count. Set in runSimulation.m by parameter `TIME_OR_COUNT` |
|finalCellNo| Cell count to reach for the program to stop. Set in runSimulation.m by parameter `FINAL_CELL_COUNT`|
|finalTime| Final number of timesteps to reach for the program to stop. Set in runSimulation.m by parameter `FINAL_TIMESTEP`|

### Overview of the mutation settings

Mutation settings are stored and input into colonySimulation.m (`mutationSettings` argument) and the colony constructor of Cell.m as the structure array with the following fields.

| parameter          | description                                                  |
| ------------------ | ------------------------------------------------------------ |
| mutationProb | A float between 0.0 and 1.0 giving the probability of a cell mutating during division. Set in runSimulation.m by parameter `MUTATION_PROB` |
| switchUpProb | A float between 0.0 and 1.0 giving the probability of a cell's phenotype increasing in resistance: in SNG mode it is the probability of a cell becoming non-genetically resistant, in megaPlate mode it is the probability of a cell's phenotype number increasing. Set in runSimulation.m by parameter `SWITCH_UP_PROB` |
|switchDownProb | A float between 0.0 and 1.0 giving the probability of a cell's phenotype decreasing in resistance: in SNG mode it is the probability of a cell losing non-genetic resistance, in megaPlate mode it is the probability of a cell's phenotype number decreasing. Set in runSimulation.m by parameter `SWITCH_DOWN_PROB` |
| mutationType | A str which sets the mutation model: `'SNG'` for the SNG model, `'megaPlate'` for the mega-plate model. Set in runSimulation.m by parameter `MUTATION_TYPE` |
| inheritPhenotype | 0 for false, 1 for true. For SNG mode, true means a non-genetically resistant cell's daughter will also be non-genetically resistant. In mega-plate mode, true means a daughter cell will inherit its mother cell's phenotype modifier. Set in runSimulation.m by parameter `INHERIT_PHENOTYPE` |
| NDivisionProb | A float between 0.0 and 1.0 giving the probability of a non-genetically resistant cell dividing in the presence of a drug in the 'SNG' mode. Set in runSimulation.m by parameter `N_DIVISION_PROB` |

### Overview of display settings

Display settings are stored and input in colonySimulation.m (`displaySettings` argument) as the structure array with the following fields.

| parameter          | description                                                  |
| ------------------ | ------------------------------------------------------------ |
| displayType | A str which determines how the colony is visualied. Three options: `'show'` `'save'` `'none'`. See detailed description of different display options below. Set in runSimulation.m by parameter `DISPLAY_TYPE`|
| folderName | A str giving the folder where results will be save. Set in runSimulation.m by `FOLDER_NAME` |
| map | A three-column matrix of values vetween 0.0 and 1.0 creating a custom colormap for visualing the colony and the different cell genotypes/phenotypes and drug concentrations. Set in runSimulation.m by parameter `map` |

If `displayType` is set to `'show'`, an updating image of the state of the colony will update as the simulation runs. The command window will print the appearance, establishment, and fixation times of mutations as they occur. runSimulation.m will print the final total cell count and final simulation time (in timesteps). In SNG mode, runSimulation.m will also display a line plot of the cell counts of the different cell types over time at the end of the simulation and print the final counts of each type of cell. In mega-plate mode, it will display the final counts of the number of cells of each genotype, the first appearance, establishment, and fixation times of each genotype (in timesteps), and line plots of the average genotype and phenotype levels of the colony over time (in timesteps).

If `displayType` is set to `'save'`, the updating images of the state of the colony will be saved as .jpg files to the specified folder under the name `imageT.jpg` where T is the number of the frame. In addition, in SNG mode a .csv file containing the cells counts of each cell type over time will be saved imder the name `cellCounts.csv`.

If `displayType` is set to `'none'`, nothing will be visualised or saved, colonySimulation.m will simply output the normal function outputs.

### Overview of diffusion settings

Diffusion settings are stored and input in colonySimulation.m (`diffusionRates` argument) as the structure array with the following fields.

| parameter          | description                                                  |
| ------------------ | ------------------------------------------------------------ |
| nutrientDiff | The number of steps a single nutrient packet will take in its random walk, correlates to the diffusion coefficient of the media. Default 10. Set in runSimulation.m by parameter `NUTRIENT_STEPS` |
| drugDiff | The number of steps a single nutrient packet will take in its random walk in SNG mode. Note in mega-plate mode drugs will not diffuse. Set in runSimulation.m by parameter `DRUG_STEPS` |

### Overview of lattice settings

Lattice size is stored and input in colonySimulation.m (`latticeSize` argument) as the structure array with the following fields.

| parameter          | description                                                  |
| ------------------ | ------------------------------------------------------------ |
| x | x dimension size, ie number of columns in the lattice. Set in runSimulation.m by parameter `xlatticeSize` |
| y | y dimension size, ie number of rows in the lattice. Set in runSimulation.m by parameter `ylatticeSize` |

There are four matrices that represent different aspects of the state of the lattice. These are stored and input in colonySimulation.m (`lattices` argument) as the structure array with the following fields.

| parameter          | description                                                  |
| ------------------ | ------------------------------------------------------------ |
| cells | This is a matrix a structure array with a single field, `'cell'`. This field contains either a a zero if no cell is present at that lattice point, or a Cell object if one is present. This Cell object stores all of the cell's information. See Cell.m (or call help(Cell)) for more details on Cell class. Note a structure array with a single field is used at each point so that the lattice can contain mixed data types.|
| nutrients | This is a matrix storing the number of nutrient packets located at each lattice point. |
| states | This is a matrix of ones and zeroes which represent whether cell is present at a given point (1) or not (0) |
| drugs | This is a matrix storing the number of drug packets located at each lattice point. |

There is one other matrix used throughout colonySimulation.m : `colouredStateLattice` This is a matrix of integers which correspond to different colours via the colormap `map` set in the display settings. The function `createColouredLattice` creates/updates this lattice, assigning numbers based on the mutation mode, drug concentrations, whether a cell is present, and the genotype/phenotype of that cell. In mega-plate mode, the function `assignDrugColours` is used to assign numbers based on drug concentration to lattice points where no cells are found. Further details on these functions can be found in their descriptions in colonySimulation.m

### Overview of initial cell layouts

The initial layout of seed cell location and state must be fed into colonySimulation.m with the cell and state lattices in `lattices` (see above). runSimulation.m is already set up to create one of six preset layouts, which can be set with the `SEED_TYPE` parameter. `SEED_TYPE` can be assigned one of the following strings.

| mode          | description                                                  |
| ------------------ | ------------------------------------------------------------ |
| `'centre'` | Centre mode seeds the initial lattice with a vertical line of five cells in the centre column of the lattice. All cells are susceptible, and have genotype = 1 and phenotype = 0. |
| `'leftSide'` | Left-side mode seeds the initial lattice with cells along the entire left-most column. All cells are susceptible, and have genotype = 1 and phenotype = 0. Ideal for mega-plate mode with the default drug layout (see below). |
| `'bothSides'` | Both-sides mode seeds the initial lattice with cells along both the left- and right-most columns. All cells are susceptible, and have genotype = 1 and phenotype = 0. Ideal for mega-plate mode with an altered drug layout where drug concentrations are low on both sides of the lattice and increase going inwards. Runs slow. |
| `'cross'` | Cross mode seeds the initial lattice with a vertical line of ten cells in the centre column crossed with a horizonal line of ten cells in the centre row. Approximately 10% of cells are non-genetically resistant, the rest are susceptible. All cells have genotype = 1 and phenotype = 0. Ideal for SNG mode. |
| `'square'` | Square mode seeds the initial lattice with a solid 20x20 square of cells in the centre of the lattice. Approximately 20% of the cells are non-genetically resistant, the rest are susceptible. All cells have genotype = 1 and phenotype = 0. Ideal for SNG mode. |
| `'checkerboard'` | Checkerboard mode seed the initial lattice with a 40x40 checkerboard of cells in  the centre of the lattice. Approximately 20% of the cells are non-genetically resistant, the rest are susceptible. All cells have genotype 1 and phenotype = 0. Ideal for SNG mode. |

### Overview of drug settings

Drug diffusion is set in diffusion settings, see above. The argument `DRUGS_TYPE` colonySimulation.m can take one of three strings, which each set the simulation to run a different drug scenario.

| mode          | description                                                  |
| ------------------ | ------------------------------------------------------------ |
| `'static'` | This will model a static drug. In this case, the presence of a drug changes the probability of the cell dividing, but the cell will not die. |
|`'cidal'` | This will model a cidal drug. In this case, the presence of a drug will kill the cell according to a given probability of survival. |
|`'none'` | This should be set when drugs are being excluded from the simulation. |

#### Note on drug concentrations

In SNG mode, the number of drug packets at a lattice point will not affect the cell's response to the drug, ie., a single drug packet will have the same effect as 100 packets. However, if drugs are set to diffuse, then a smaller number of packets a point increases the chance that all those packets will diffuse away from that point, leaving the point drug-free. 

In mega-plate mode, the number of drug packets at a lattice point will affect the cell's response to the drug (see above formula in description of mega-plate model).

#### Initial drug lattice layouts

SNG mode and mega-plate mode each have a different default layout of drug concentrations when using runSimulation.m. In SNG mode, every lattice point is given the same number of drug packets. This number is set by the parameter `startDrugs`. In mega-plate mode, the lattice is split into five bands of equal width. The first band contains no drugs, then from there the k-1th band contains $ 10^{k-1} $ drug packets at each lattice point within that band. Based on the default colormap set in runSimulation.m, these bands that contain drugs are represented by shades of grey which get progressively darker as drug concentrations go up. 

### Overview of nutrient settings

The class constructor function for Cell.m requires the argument `nutrForDivision`, which sets the number of nutrient packets a cell must consume before it is able to divide. This parameter is set in runSimulation.m by `NUTRS_FOR_DIVISION`. In runSimulation.m, the initial nutrient lattice is set to have the same number of nutrient packets at every lattice point. The number of nutrient packets is set by the parameter `startNutrs`.

### Code authors
Code was written by Rebekah Hall under the supervision of Dr. Daniel Charlebois. 





