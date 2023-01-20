classdef Cell
    % CELL Stores all necessary information and actions of a single cell in
    % colonySimulation.m
    %   A Cell object represents a single cell that can divide, consume
    %   nutrients, mutate, and change phenotypes.
    properties
        state = 1;   % 1 if the cell is alive, 0 if the cell has died and must be removed from the lattice.
        consumedNutrients = 0;   % indicates how many nutrients the cell has consumed since dividing last.
        nutrForDivision = 1;   % number of nutrients the cell must consume before dividing.
        canDivide = true;   % indicates whether the cell is able to divide at this moment. 
        mutated = false;   % bool determining whether the cell has mutated.
        NGresistant = false;   % bool indicating whether the cell is non-genetically resistant in SNG mode.
        mutationProb = 0.0;   % float between 0 and 1 giving the probability of mutation
        switchUpProb = 0.0;   % float between 0.0 and 1.0 giving the probability of a cell's phenotype increasing in resistance: in SNG mode it is the probability of a cell becoming non-genetically resistant, in megaPlate mode it is the probability of a cell's phenotype number increasing
        switchDownProb = 0.0;   % float between 0.0 and 1.0 giving the probability of a cell's phenotype increasing in resistance: in SNG mode it is the probability of a cell losing non-genetic resistance, in megaPlate mode it is the probability of a cell's phenotype number decreasing
        mutationType = 'SNG';   % str which sets the mutation model: 'SNG' = cells can be susceptible, non-genetically drug resistant, or genetically drug resistant; 'megaPlate' = cells are assigned a genotype which sets the level of resistance and a phenotype number which is added to the genotype number to either increase or decrease the resistance
        inheritPhenotype = false;   % bool setting whether a daughter cell inherits its mother's phenotype
        genotype = 1;   % genotype number which determines resistance level in megaPlate mode
        phenotype = 0;   % phenotype number, either -1, 0, 1, which is added to genotype number to vary resistance levels among cells of the same genotype
        divisionProb = 1.0;   % float betwen 0.0 and 1.0 which sets the probability the cell will divide
        NDivisionProb = 0.5;   % float between 0.0 and 1.0 which sets the probability that a non-genetically drug resistant cell will divide or survive in the presence of a static or cidal drug respectively in SNG mode
    end
    methods     
        function obj = Cell(nutrForDivision,mutationSettings)
            % CELL Creates a new Cell object. 
            %    INPUTS
            %        nutrForDivision sets the number of nutrients needed for a cell to divide
            %        mutationSettings is a structure array containing the following fields:
            %            -- 'mutationProb' is a float between 0.0 and 1.0 giving the probability
            %                of a cell mutating during division
            %            -- 'mutationType' is a str which sets the mutation model:
            %               'SNG' = cells can be susceptible, non-genetically drug
            %                resistant, or genetically drug resistant
            %               'megaPlate' = cells are assigned a genotype which sets the level
            %                of resistance and a phenotype number which is added to the
            %                genotype number to either increase or decrease the resistance
            %            -- 'switchUpProb' is a float between 0.0 and 1.0 giving the probability
            %                of a cell's phenotype increasing in resistance: in SNG mode it is the
            %                probability of a cell becoming non-genetically resistant, 
            %                in megaPlate mode it is the probability of a cell's phenotype 
            %                number increasing
            %            -- 'switchDownProb' is a float between 0.0 and 1.0 giving the probability
            %                of a cell's phenotype decreasing in resistance: in SNG mode it is the
            %                probability of a cell losing non-genetic resistance, 
            %                in megaPlate mode it is the probability of a cell's phenotype
            %                number decreasing
            %            -- 'inheritPhenotype' is an int. 0 if the daughter cell does not
            %                inherit the phenotype of the mother, 1 if the daughter cell does
            %                inherit the phenotype of the mother.
            %            -- 'NDivisionProb' is a float between 0.0 and 1.0 giving the
            %                probability of a non-genetically resistant cell dividing in the
            %                presence of a drug in the 'SNG' mode.
           
             obj.nutrForDivision = nutrForDivision;
             obj.mutationProb = mutationSettings.('mutationProb');
             obj.switchUpProb = mutationSettings.('switchUpProb');
             obj.switchDownProb = mutationSettings.('switchDownProb');
             obj.mutationType = mutationSettings.('mutationType');
             obj.inheritPhenotype = mutationSettings.('inheritPhenotype');
             obj.NDivisionProb = mutationSettings.('NDivisionProb');
        end
        
        function obj = consumeNutrient(obj)
            % CONSUMENUTRIENT Updates consumedNutrients based on the cell consuming one nutrient
            
            obj.consumedNutrients = obj.consumedNutrients + 1;
        end
        
        
        function obj = updateDivisibility(obj)
            % UPDATEDIVISIBILITY Updates whether the cell can bud or not
            
            if obj.consumedNutrients < obj.nutrForDivision
                obj.canDivide = false;
            else
                obj.canDivide = true;
            end
        end
        
        function obj = respondToDrug(obj,drugConc,drugType)
            % RESPONDTODRUG Update cell state according to drug type, concentration, and mutation mode
            %   INPUTS
            %       drugConc is the number of drug packets at the lattice
            %           site where the cell is located
            %       drugType is the str indicating what type of drug is present:
            %           'cidal' = cidal drug
            %           'static' = static drug
            
            % calculate resistance based on genotype and phenotype for
            % megaPlate mode
            resistance = obj.genotype + obj.phenotype;
            phi = max([0,1-(drugConc/10^resistance)^2]);
            
            if strcmp(drugType,'cidal')
                if strcmp(obj.mutationType,'SNG')
                    if drugConc > 0
                        if ~obj.mutated
                            if obj.NGresistant
                                p = rand();
                                if p > obj.NDivisionProb
                                    obj.state = 0;
                                end
                            else
                                obj.state = 0;
                            end
                        end
                    end
                elseif strcmp(obj.mutationType,'megaPlate')
                    p = rand();
                    if p > phi
                        obj.state = 0;
                    end
                end
            elseif strcmp(drugType,'static')
                if strcmp(obj.mutationType,'SNG')
                    if drugConc > 0
                        if obj.mutated
                            obj.divisionProb = 1;
                        elseif obj.NGresistant
                            obj.divisionProb = obj.NDivisionProb;
                        else
                            obj.divisionProb = 0;
                        end
                    else
                        obj.divisionProb = 1;
                    end
                elseif strcmp(obj.mutationType,'megaPlate')
                    obj.divisionProb = phi;
                end
            end
        end
        
        function [divided,lattices,cellNumbers] = divide(obj,lattices,coordinates,cellNumbers)
            % DIVIDE Attempt to create daughter cell.
            %   INPUTS
            %       lattices is a structure array containing the following fields:
            %            -- 'cells' is a matrix of structure arrays with the following field:
            %                      + 'cell' is 0 if no cell is present, a Cell object if
            %                         a cell is present
            %            -- 'nutrients' is a matrix containing the number of nutrient packets
            %                at each point on the lattice
            %            -- 'states' is a matrix of ones and zeros, where a 1 states that a
            %                cell is a present and 0 states that no cell is present
            %            -- 'drugs' is a matrix containing the number of drug packets at each
            %                point on the lattice
            %       coordinates is 1x2 array containing the row then column where the mother cell is located
            %       cellNumbers is a structrue array containing the total number of cells
            %           and for 'SNG' mode, the number of cells of each SNG type.
            %           Contains the following fields:
            %               -- 'cellsno' is the total number of living cells
            %               -- 'susceptibles' is the number of susceptible cells
            %               -- 'nongenet' is the number of non-genetically drug
            %                   resistant cells
            %               -- 'genet' is the number of genetically drug resistant cells
            %   OUTPUTS
            %       divided is a bool indicating whether the cell divided
            %       lattices is the updated struct of that from INPUTS
            %       cellNumbers is the updated struct of that from OUTPUTS
            %
            %   See also EMPTYCELLARRAY, UPDATEMOTHERDAUGHTER, UPDATEDIVISIBILITY, RESPONDTODRUG.
            
            divided = false;
            if obj.canDivide
                p = rand();
                if p < obj.divisionProb
                    i = coordinates(1);
                    j = coordinates(2);
                    emptyCells = Cell.emptyCellArray(lattices.('states'),coordinates);
                    [numRows,~] = size(emptyCells);
                    randInd = randi(numRows,1);
                    chosenSitei = emptyCells(randInd,1);
                    chosenSitej = emptyCells(randInd,2);
                    
                    [lattices,cellNumbers] = obj.updateMotherDaughter(lattices,cellNumbers,i,j,chosenSitei,chosenSitej);
                    divided = true;
                    cellNumbers.('cellsno') = cellNumbers.('cellsno') + 1;
                end
            end
        end
        
        function obj = updatePhenotype(obj)
            % UPDATEPHENOTYPE Updates the phenotype of a cells randomly.
            
            % SNG mode
            if strcmp(obj.mutationType,'SNG')
                p = rand();
                if ~obj.NGresistant
                    if p < obj.switchUpProb
                        obj.NGresistant = true;
                    end
                elseif obj.NGresistant
                    if p < obj.switchDownProb
                        obj.NGresistant = false;
                    end
                end
            % megaPlate mode
            elseif strcmp(obj.mutationType,'megaPlate')
                if obj.phenotype == 0
                    p1 = rand();
                    if p1 < obj.switchUpProb
                        obj.phenotype = 1;
                    elseif p1 < obj.switchDownProb + obj.switchUpProb
                        obj.phenotype = -1;
                    end
                elseif obj.phenotype == -1
                    p2 = rand();
                    if p2 < obj.switchUpProb
                        obj.phenotype = 0;
                    end
                else
                    p3 = rand();
                    if p3 < obj.switchDownProb
                        obj.phenotype = 0;
                    end
                end
            end
        end
        
        function obj = mutate(obj,motherCell)
            % MUTATE Decides whether a mutation will occur in the daughter cell.
            %   INPUTS
            %       motherCell is the Cell object representing the daughter cell's mother
            %
            %   See also UPDATEMOTHERDAUGHTER.
            
            if strcmp(obj.mutationType,'SNG')
                if motherCell.mutated
                    obj.mutated = true;
                elseif obj.NGresistant
                    p = rand();
                    if p < obj.mutationProb
                        obj.mutated = true;
                    end
                else
                    p2 = rand();
                    if p2 < obj.mutationProb/2
                        obj.mutated = true;
                    end
                end
            elseif strcmp(obj.mutationType,'megaPlate')
                p1 = rand();
                if p1<obj.mutationProb
                    obj.genotype = motherCell.genotype + 1;
                else
                    obj.genotype = motherCell.genotype;
                end
            end
        end   
        
        function [lattices,cellNumbers] = updateMotherDaughter(obj,lattices,cellNumbers,i,j,chosenSitei,chosenSitej)
            % UPDATEMOTHERDAUGHTER Updates a new daughter cell and its mother after division.
            %   INPUTS
            %       lattices is a structure array containing the following fields:
            %            -- 'cells' is a matrix of structure arrays with the following field:
            %                      + 'cell' is 0 if no cell is present, a Cell object if
            %                         a cell is present
            %            -- 'nutrients' is a matrix containing the number of nutrient packets
            %                at each point on the lattice
            %            -- 'states' is a matrix of ones and zeros, where a 1 states that a
            %                cell is a present and 0 states that no cell is present
            %            -- 'drugs' is a matrix containing the number of drug packets at each
            %                point on the lattice
            %       cellNumbers is a structrue array containing the total number of cells
            %           and for 'SNG' mode, the number of cells of each SNG type.
            %           Contains the following fields:
            %               -- 'cellsno' is the total number of living cells
            %               -- 'susceptibles' is the number of susceptible cells
            %               -- 'nongenet' is the number of non-genetically drug
            %                   resistant cells
            %               -- 'genet' is the number of genetically drug resistant cells
            %       i is the row number of the mother cell
            %       j is the column number of the mother cell
            %       chosenSitei is the row number of the new daughter cell
            %       chosenSitej is the column number of the new daughter cell
            %   OUTPUTS
            %       lattices is the updated struct of that from INPUTS
            %       cellNumbers is the updated struct of that from OUTPUTS
            %
            %   See also MUTATE, UPDATEPHENOTYPE, DIVIDE.
            
            % Set mother cell's consumedNutrients to 0.
            lattices.('cells')(i,j).('cell').consumedNutrients = 0;
            
            % Copy mother cell object to create new daughter cell.
            lattices.('cells')(chosenSitei,chosenSitej).('cell') = obj;
            lattices.('states')(chosenSitei,chosenSitej) = 1;
            
            if obj.inheritPhenotype == 0
                lattices.('cells')(chosenSitei,chosenSitej).('cell').phenotype = 0;
                lattices.('cells')(chosenSitei,chosenSitej).('cell').NGresistant = false;
            end
            
            lattices.('cells')(chosenSitei,chosenSitej).('cell') = lattices.('cells')(chosenSitei,chosenSitej).('cell').mutate(obj);
            
            lattices.('cells')(chosenSitei,chosenSitej).('cell') = lattices.('cells')(chosenSitei,chosenSitej).('cell').updatePhenotype();
            
            % Update SNG cell counts
            if lattices.('cells')(chosenSitei,chosenSitej).('cell').mutated
                cellNumbers.('genet') = cellNumbers.('genet') + 1;
            elseif lattices.('cells')(chosenSitei,chosenSitej).('cell').NGresistant
                cellNumbers.('nongenet') = cellNumbers.('nongenet') + 1;
            else
                cellNumbers.('susceptibles') = cellNumbers.('susceptibles') + 1;
            end
        end
    end
    
    methods(Static)
        
        function emptyCells = emptyCellArray(stateLattice,coordinates)
            % EMPTYCELLARRAY Gives an array containing the indices of every empty lattice point surrounding a cell at (i,j). 
            %   INPUTS
            %       stateLattice is the lattice of 0s and 1s indicating where living cells are located
            %       coordinates is 1x2 array containing the row then column where the mother cell is located
            %   OUTPUTS
            %       emptyCells is an array containing the indices of every empty 
            %       lattice point surrounding a cell at (i,j) with dimensions emptySpots x 2,
            %       where emptySpots is the number of empty lattice points.
            %   
            %   See also DIVIDE.
            
            latticeSizes = size(stateLattice);
            xlatticeSize = latticeSizes(2);
            ylatticeSize = latticeSizes(1);
            
            ind = 1; 
            
            i = coordinates(1);
            j = coordinates(2);
            
            % Find the number of empty sites surrounding the cell.
            emptySpots = 0;
            for I = -1:1
                for J = -1:1
                    if I ~= 0 || J ~= 0
                        if i+I <= ylatticeSize && i+I > 0 && j+J <= xlatticeSize && j+J > 0
                            if stateLattice(i+I,j+J) == 0
                                emptySpots = emptySpots + 1;
                                ind = ind + 1;
                            end
                        end
                    end
                end
            end  
            
            emptyCells = zeros(emptySpots,2);
            ind = 1;

            % Create array of the coordinates of empty lattice points.
            for I = -1:1
                for J = -1:1
                    if I ~= 0 || J ~= 0
                        if i+I <= ylatticeSize && i+I > 0 && j+J <= xlatticeSize && j+J > 0
                            if stateLattice(i+I,j+J) == 0
                                emptyCells(ind,1) = i+I;
                                emptyCells(ind,2) = j+J;
                                ind = ind + 1;
                            end
                        end
                    end
                end
            end    
        end
    end 
end
