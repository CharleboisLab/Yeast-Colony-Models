classdef YeastCell
    properties
        state = 1
        budScar = [0,0]   % gives a vector in the direction of the bud scar of a cell at that lattice site.
        consumedNutrients = 0   % indicates how many nutrients the cell has consumed since budding last.
        nutrForBudding = 1;
        unipolarFrac = 1   % indicates the fraction of the time the cell will bud regularly.
        axialFrac = 0.6
        deltaUnipolarFrac = 0.125   % indicates by how much unipolarFrac will change as the number of nutrient packets changes.
        buddable = true;   % indicates whether the cell is able to bud at this moment. 
        budded = false;   % bool indicates whether the cell has previously budded.
        budProb = 1.0;   % probability that the cell will bud.
        mutant = 1;   % strain number. 1 indicates FLO11 WT. 2 indicates flo11 knockout.
        mutationProb = 0.1;   % probability of mutation.
    end
    methods     
        function obj = YeastCell(budScar,axialFrac,nutrForBudding,mutant,mutationProb)
           % class constructor
             obj.budScar = budScar;
             obj.axialFrac = axialFrac;
             obj.nutrForBudding = nutrForBudding;
             obj.deltaUnipolarFrac = 1/(nutrForBudding*8);
             obj.mutant = mutant;
             obj.mutationProb = mutationProb;
        end
        
        function cellArray = interactions(obj,cellArray,coordinates,magnet)
            % magnetic field interactions
            strength = magnet.strength;
            p = rand();
            if p < strength
                cellArray = YeastCell.EMFBias(cellArray,coordinates,magnet);
            end
        end
        
        function obj = mutate(obj,lattices,motherCoordi,motherCoordj)
            mutationP = lattices.('cells')(motherCoordi,motherCoordj).('cell').mutationProb;
            p = rand();
            if p < mutationP
                obj.mutant = 2;
            end
        end
        
        function obj = consumeNutrient(obj)
            % updateNutrients updates the nutrientConcentration and
            % consumedNutrients based on the cell consuming one nutrient
            obj.consumedNutrients = obj.consumedNutrients + 1;
        end
        
        function obj = updateBudProb(obj,agarWeight)
            % updateBudProb updates the probability of budding based on the
            % flo11 genotype. Can edit to include dependance on agar weight
            
            % WT case
            if obj.mutant == 1
                obj.budProb = 0.75;
            % mutant case
            elseif obj.mutant == 2
                obj.budProb = 1.0;
            end
       
        end
        
        function obj = updateUnipolarFrac(obj,nutrients)
            % updateUnipolarFrac updates unipolarFrac based on
            % nutrientConcentration
            obj.unipolarFrac = obj.deltaUnipolarFrac * nutrients;
        end
        
        function obj = updateBuddability(obj)
            % updateBuddability updates whether the cell can bud or not
            if isequal(obj.budScar,[0 0])
                obj.buddable = false;
            elseif obj.consumedNutrients < obj.nutrForBudding
                obj.buddable = false;
            else
                obj.buddable = true;
            end
        end
        
        function [unipolar,bud,lattices,cellsno,genoCellsno] = bud(obj,lattices,coordinates,cellsno,genoCellsno,magnet)
            bud = false;
            unipolar = false;
            if obj.buddable
                p0 = rand();
                if p0 < obj.budProb
                    p1 = rand();
                    if p1 < obj.unipolarFrac
                        if obj.budded == false
                            p2 = rand();
                            if p2 < obj.axialFrac
                                [bud,genoNo,lattices] = obj.axialBud(lattices,coordinates,magnet);
                            else
                                [bud,genoNo,lattices] = obj.polarBud(lattices,coordinates,magnet);
                            end
                        else
                            p3 = rand();
                            if p3 < (1+obj.axialFrac)/2
                                [bud,genoNo,lattices] = obj.axialBud(lattices,coordinates,magnet);
                            else
                                [bud,genoNo,lattices] = obj.polarBud(lattices,coordinates,magnet);
                            end
                        end
                    else
                        [bud,genoNo,lattices] = obj.unipolarBud(lattices,coordinates,magnet);
                        unipolar = true;
                    end
                end
                if bud
                    cellsno = cellsno + 1;
                    totalGenotypes = length(genoCellsno);
                    for G = 1:totalGenotypes
                        if genoNo == G
                            genoCellsno(G) = genoCellsno(G) + 1;
                        end
                    end
                end
            end
        end
        
        function [bud,genoNo,lattices] = axialBud(obj,lattices,coordinates,magnet)
            bud = false;
            i = coordinates(1);
            j = coordinates(2);
            genoNo = lattices.('cells')(i,j).('cell').mutant;
            lattices = obj.assignAngles(lattices,coordinates);
            nearbyEmpties = YeastCell.findNearby(lattices,[i,j],true);
            if ~isempty(nearbyEmpties)
                
                nearbyEmpties = obj.interactions(nearbyEmpties,coordinates,magnet);

                % Randomly select one of these instances as the location for the
                % new daughter cell.
                [numRows,~] = size(nearbyEmpties);
                randInd = randi(numRows,1);
                chosenSitei = nearbyEmpties(randInd,1);
                chosenSitej = nearbyEmpties(randInd,2);
                
                lattices = obj.updateMotherDaughter(lattices,i,j,chosenSitei,chosenSitej);
                bud = true;
                genoNo = lattices.('cells')(chosenSitei,chosenSitej).('cell').mutant;
                
            % If both sites 45 degrees away from the pevious bud site are full,
            % the cell will bud into one of the sites, pushing the existing cell out
            % of the way. If there are no empty site surround the existing cell,
            % nothing will happen.
            else
                % Determine the indices for every instance of the minimum angle.
                nearbySites = YeastCell.findNearby(lattices,[i,j],false);
                if ~isempty(nearbyEmpties)
                    nearbySites = obj.interactions(nearbySites,coordinates,magnet);

                    % Randomly select one of these instances as the location for the
                    % new daughter cell.
                    [numRows,~] = size(nearbySites);
                    randInd = randi(numRows,1);
                    chosenSitei = nearbySites(randInd,1);
                    chosenSitej = nearbySites(randInd,2);
                    % Randomly select an empty site surround the preexisting cell to
                    % move into.
                    stateLattice = lattices.('states');
                    otherEmptyCoords = YeastCell.emptyCellArray(stateLattice,[chosenSitei,chosenSitej]);

                    if ~isempty(otherEmptyCoords)
                        otherEmptyCoords = obj.interactions(otherEmptyCoords,[chosenSitei,chosenSitej],magnet);
                        [numRows,~] = size(otherEmptyCoords);
                        randInd = randi(numRows,1);
                        pushedOutSitei = otherEmptyCoords(randInd,1);
                        pushedOutSitej = otherEmptyCoords(randInd,2);

                        lattices.('cells')(pushedOutSitei,pushedOutSitej).('cell') = lattices.('cells')(chosenSitei,chosenSitej).('cell');
                        lattices.('states')(pushedOutSitei,pushedOutSitej) = lattices.('states')(chosenSitei,chosenSitej);

                        lattices = obj.updateMotherDaughter(lattices,i,j,chosenSitei,chosenSitej);
                        bud = true;
                        genoNo = lattices.('cells')(chosenSitei,chosenSitej).('cell').mutant;
                    end
                end
            end
        end
        
        function [bud,genoNo,lattices] = polarBud(obj,lattices,coordinates,magnet)
            bud = false;
            i = coordinates(1);
            j = coordinates(2);
            genoNo = lattices.('cells')(i,j).('cell').mutant;
            lattices = obj.assignAngles(lattices,coordinates);
            polarEmpties = YeastCell.findPolar(lattices,[i,j],true);
            if ~isempty(polarEmpties)
                
                polarEmpties = obj.interactions(polarEmpties,coordinates,magnet);

                % Randomly select one of these instances as the location for the
                % new daughter cell.
                [numRows,~] = size(polarEmpties);
                randInd = randi(numRows,1);
                chosenSitei = polarEmpties(randInd,1);
                chosenSitej = polarEmpties(randInd,2);
                
                lattices = obj.updateMotherDaughter(lattices,i,j,chosenSitei,chosenSitej);
                bud = true;
                genoNo = lattices.('cells')(chosenSitei,chosenSitej).('cell').mutant;
            % If both sites 45 degrees away from the pevious bud site are full,
            % the cell will bud into one of the sites, pushing the existing cell out
            % of the way. If there are no empty site surround the existing cell,
            % nothing will happen.
            else
                % Determine the indices for every instance of the minimum angle.
                polarSites = YeastCell.findPolar(lattices,[i,j],false);
                if ~isempty(polarSites)

                    polarSites = obj.interactions(polarSites,coordinates,magnet);

                    % Randomly select one of these instances as the location for the
                    % new daughter cell.
                    [numRows,~] = size(polarSites);
                    randInd = randi(numRows,1);
                    chosenSitei = polarSites(randInd,1);
                    chosenSitej = polarSites(randInd,2);
                    % Randomly select an empty site surround the preexisting cell to
                    % move into.
                    stateLattice = lattices.('states');
                    otherEmptyCoords = YeastCell.emptyCellArray(stateLattice,[chosenSitei,chosenSitej]);

                    if ~isempty(otherEmptyCoords)
                        otherEmptyCoords = obj.interactions(otherEmptyCoords,[chosenSitei,chosenSitej],magnet);
                        [numRows,~] = size(otherEmptyCoords);
                        randInd = randi(numRows,1);
                        pushedOutSitei = otherEmptyCoords(randInd,1);
                        pushedOutSitej = otherEmptyCoords(randInd,2);

                        lattices.('cells')(pushedOutSitei,pushedOutSitej).('cell') = lattices.('cells')(chosenSitei,chosenSitej).('cell');
                        lattices.('states')(pushedOutSitei,pushedOutSitej) = lattices.('states')(chosenSitei,chosenSitej);

                        lattices = obj.updateMotherDaughter(lattices,i,j,chosenSitei,chosenSitej);
                        bud = true;
                        genoNo = lattices.('cells')(chosenSitei,chosenSitej).('cell').mutant;
                    end
                end
            end
        end
        
        function [bud, genoNo, lattices] = unipolarBud(obj,lattices,coordinates,magnet)
            bud = false;
            i = coordinates(1);
            j = coordinates(2);
            genoNo = lattices.('cells')(i,j).('cell').mutant;
            lattices = obj.assignAngles(lattices,coordinates);
            unipolarEmpties = YeastCell.findUnipolar(lattices,[i,j]);
            if ~isempty(unipolarEmpties)
                unipolarEmpties = obj.interactions(unipolarEmpties,coordinates,magnet);
                
                % Randomly select one of these instances as the location for the
                % new daughter cell.
                [numRows,~] = size(unipolarEmpties);
                randInd = randi(numRows,1);
                chosenSitei = unipolarEmpties(randInd,1);
                chosenSitej = unipolarEmpties(randInd,2);
                
                lattices = obj.updateMotherDaughterUnipolar(lattices,i,j,chosenSitei,chosenSitej);
                bud = true;
                genoNo = lattices.('cells')(chosenSitei,chosenSitej).('cell').mutant;
            else
                
            end
        end
        
        function lattices = updateMotherDaughter(obj,lattices,i,j,chosenSitei,chosenSitej)
            lattices.('cells')(i,j).('cell').consumedNutrients = 0;
            
            daughterBudScar = [i - chosenSitei,j - chosenSitej];
            lattices.('cells')(chosenSitei,chosenSitej).('cell') = obj;
            lattices.('cells')(chosenSitei,chosenSitej).('cell') = lattices.('cells')(chosenSitei,chosenSitej).('cell').mutate(lattices,i,j);
            lattices.('cells')(chosenSitei,chosenSitej).('cell').budScar = daughterBudScar;
            lattices.('cells')(chosenSitei,chosenSitej).('cell').budded = false;
            lattices.('states')(chosenSitei,chosenSitej) = 1;
            lattices.('mutants')(chosenSitei,chosenSitej) = lattices.('cells')(chosenSitei,chosenSitej).('cell').mutant;

            lattices.('cells')(i,j).('cell').budScar = [chosenSitei - i,chosenSitej - j];
            lattices.('cells')(i,j).('cell').budded = true;
        end
        
        function lattices = updateMotherDaughterUnipolar(obj,lattices,i,j,chosenSitei,chosenSitej)
            lattices.('cells')(i,j).('cell').consumedNutrients = 0;
            
            secondSitei = 2*chosenSitei - i;
            secondSitej = 2*chosenSitej - j;
            daughterBudScar = [i - chosenSitei,j - chosenSitej];
            cellEnd1 = obj;
            cellEnd1 = cellEnd1.mutate(lattices,i,j);
            cellEnd1.budScar = [0 0];
            cellEnd1.buddable = false;
            cellEnd1.budded = false;
            
            cellEnd2 = cellEnd1;
            cellEnd2.budScar = daughterBudScar;
            lattices.('cells')(chosenSitei,chosenSitej).('cell') = cellEnd1;
            lattices.('cells')(secondSitei,secondSitej).('cell') = cellEnd2;
            lattices.('states')(chosenSitei,chosenSitej) = 1;
            lattices.('states')(secondSitei,secondSitej) = 1;
            lattices.('mutants')(chosenSitei,chosenSitej) = cellEnd1.mutant;
            lattices.('mutants')(secondSitei,secondSitej) = cellEnd2.mutant;

            lattices.('cells')(i,j).('cell').budded = true;
        end
        
        function lattices = assignAngles(obj,lattices,coordinates)

            % createAngleArray creates an array containing the angle between each
            % lattice site in cellArray and the lattice site (i,j).

            % create array containing angle between the vector of the
            % mother's bud scar and the vector from the mother cell to each
            % lattice site in cellArray
            latticeSizes = size(lattices.('states'));
            xlatticeSize = latticeSizes(2);
            ylatticeSize = latticeSizes(1);
            i = coordinates(1);
            j = coordinates(2);
            vector = obj.budScar;
            for I = -1:1
                for J = -1:1
                    if I ~= 0 || J ~= 0
                        if i+I <= ylatticeSize && i+I > 0 && j+J <= xlatticeSize && j+J > 0
                            dotProd = vector*[I;J];
                            norm1 = norm(vector);
                            norm2 = norm([I J]);
                            norms = norm1*norm2;
                            normDot = dotProd/norms;
                            angle = acosd(normDot);
                            lattices.('cells')(i+I,j+J).('angle') = angle;
                        end
                    end
                end
            end
        end
    end
    
    methods(Static)
        function possibleCells = EMFBias(cellArray,coordinates,magnet)
            magneticField = magnet.direction;
            minAngle = magnet.minAngle;
            maxAngle = magnet.maxAngle;
            [numRows,~] = size(cellArray);
    
            if ~isequal(magneticField,[0 0])

                newCellArray = YeastCell.assignEMFAngles(cellArray,coordinates,magneticField);

                EMFCellno = 0;
                for x = 1:numRows
                    angle = newCellArray(x,3);
                    if (angle > minAngle)&&(angle < maxAngle)
                        EMFCellno = EMFCellno + 1;
                    end
                end
                
                if EMFCellno > 0
                    possibleCells = zeros(EMFCellno,2);
                    ind = 1;
                    for x = 1:numRows
                        i = newCellArray(x,1);
                        j = newCellArray(x,2);
                        angle = newCellArray(x,3);
                        if (angle > minAngle)&&(angle < maxAngle)
                            possibleCells(ind,1) = i;
                            possibleCells(ind,2) = j;
                            ind = ind + 1;
                        end
                    end
                else
                    possibleCells = cellArray;
                end
            else
                possibleCells = cellArray;
            end
        end
        
        function possibleCoords = findNearby(lattices,coordinates,needEmpty)
            latticeSizes = size(lattices.('states'));
            xlatticeSize = latticeSizes(2);
            ylatticeSize = latticeSizes(1);
            possibleSitesNo = 0;
            i = coordinates(1);
            j = coordinates(2);
            for I = -1:1
                for J = -1:1
                    if I ~= 0 || J ~= 0
                        if i+I <= ylatticeSize && i+I > 0 && j+J <= xlatticeSize && j+J > 0
                            neighbourAngle = lattices.('cells')(i+I,j+J).('angle');
                            neighbourState = lattices.('states')(i+I,j+J);
                            if needEmpty
                                if neighbourAngle > 44 && neighbourAngle < 46 && neighbourState == 0
                                    possibleSitesNo = possibleSitesNo + 1;
                                end
                            else
                                if neighbourAngle > 44 && neighbourAngle < 46
                                    possibleSitesNo = possibleSitesNo + 1;
                                end
                            end
                        end
                    end
                end
            end
            
            if possibleSitesNo > 0
                possibleCoords = zeros(possibleSitesNo,2);
                ind = 1;
                for I = -1:1
                    for J = -1:1
                        if I ~= 0 || J ~= 0
                            if i+I <= ylatticeSize && i+I > 0 && j+J <= xlatticeSize && j+J > 0
                                neighbourAngle = lattices.('cells')(i+I,j+J).('angle');
                                neighbourState = lattices.('states')(i+I,j+J);
                                if needEmpty
                                    if neighbourAngle > 44 && neighbourAngle < 46 && neighbourState == 0
                                        possibleCoords(ind,1) = i+I;
                                        possibleCoords(ind,2) = j+J;
                                        ind = ind + 1;
                                    end
                                else
                                    if neighbourAngle > 44 && neighbourAngle < 46
                                        possibleCoords(ind,1) = i+I;
                                        possibleCoords(ind,2) = j+J;
                                        ind = ind + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            else
                possibleCoords = [];
            end                        
        end
        
        function possibleCoords = findPolar(lattices,coordinates,needEmpty)
            latticeSizes = size(lattices.('states'));
            xlatticeSize = latticeSizes(2);
            ylatticeSize = latticeSizes(1);
            possibleSitesNo = 0;
            i = coordinates(1);
            j = coordinates(2);
            for I = -1:1
                for J = -1:1
                    if I ~= 0 || J ~= 0
                        if i+I <= ylatticeSize && i+I > 0 && j+J <= xlatticeSize && j+J > 0
                            neighbourAngle = lattices.('cells')(i+I,j+J).('angle');
                            neighbourState = lattices.('states')(i+I,j+J);
                            if needEmpty
                                if neighbourAngle > 134 && neighbourAngle < 181 && neighbourState == 0
                                    possibleSitesNo = possibleSitesNo + 1;
                                end
                            else
                                if neighbourAngle > 134 && neighbourAngle < 181
                                    possibleSitesNo = possibleSitesNo + 1;
                                end
                            end
                        end
                    end
                end
            end
            
            if possibleSitesNo > 0
                possibleCoords = zeros(possibleSitesNo,2);
                ind = 1;
                for I = -1:1
                    for J = -1:1
                        if I ~= 0 || J ~= 0
                            if i+I <= ylatticeSize && i+I > 0 && j+J <= xlatticeSize && j+J > 0
                                neighbourAngle = lattices.('cells')(i+I,j+J).('angle');
                                neighbourState = lattices.('states')(i+I,j+J);
                                if needEmpty
                                    if neighbourAngle > 134 && neighbourAngle < 181 && neighbourState == 0
                                        possibleCoords(ind,1) = i+I;
                                        possibleCoords(ind,2) = j+J;
                                        ind = ind + 1;
                                    end
                                else
                                    if neighbourAngle > 134 && neighbourAngle < 181
                                        possibleCoords(ind,1) = i+I;
                                        possibleCoords(ind,2) = j+J;
                                        ind = ind + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            else
                possibleCoords = [];
            end 
        end
        
        function possibleCoords = findUnipolar(lattices,coordinates)
            latticeSizes = size(lattices.('states'));
            xlatticeSize = latticeSizes(2);
            ylatticeSize = latticeSizes(1);
            possibleSitesNo = 0;
            i = coordinates(1);
            j = coordinates(2);
            for I = -1:1
                for J = -1:1
                    if I ~= 0 || J ~= 0
                        if i+I+I <= ylatticeSize && i+I+I > 0 && j+J+J <= xlatticeSize && j+J+J > 0
                            neighbourAngle = lattices.('cells')(i+I,j+J).('angle');
                            neighbourState = lattices.('states')(i+I,j+J);
                            secondNeighbourState = lattices.('states')(i+I+I,j+J+J);
                            if neighbourAngle > 134 && neighbourAngle < 181 && neighbourState == 0 && secondNeighbourState == 0
                                possibleSitesNo = possibleSitesNo + 1;
                            end
                        end
                    end
                end
            end
            
            if possibleSitesNo > 0
                possibleCoords = zeros(possibleSitesNo,2);
                ind = 1;
                for I = -1:1
                    for J = -1:1
                        if I ~= 0 || J ~= 0
                            if i+I+I <= ylatticeSize && i+I+I > 0 && j+J+J <= xlatticeSize && j+J+J > 0
                                neighbourAngle = lattices.('cells')(i+I,j+J).('angle');
                                neighbourState = lattices.('states')(i+I,j+J);
                                secondNeighbourState = lattices.('states')(i+I+I,j+J+J);
                                if neighbourAngle > 134 && neighbourAngle < 181 && neighbourState == 0 && secondNeighbourState == 0
                                    possibleCoords(ind,1) = i+I;
                                    possibleCoords(ind,2) = j+J;
                                    ind = ind + 1;
                                end
                            end
                        end
                    end
                end
            else
                possibleCoords = [];
            end 
        end
        
        function newCellArray = assignEMFAngles(cellArray,coordinates,vector)

            % createAngleArray creates an array containing the angle between each
            % lattice site in cellArray and the lattice site (i,j).

            [numRows,~] = size(cellArray);
                
            newCellArray = zeros(numRows,3);

            % create array containing angle between the vector of the
            % mother's bud scar and the vector from the mother cell to each
            % lattice site in cellArray
            i_1 = coordinates(1);
            j_1 = coordinates(2);
            
            for x = 1:numRows
                i_2 = cellArray(x,1);
                j_2 = cellArray(x,2);
                vector2 = [i_2-i_1;j_2-j_1];
                dotProd = vector*vector2;
                norm1 = norm(vector);
                norm2 = norm(vector2);
                norms = norm1*norm2;
                normDot = dotProd/norms;
                angle = acosd(normDot);
                newCellArray(x,1) = i_2;
                newCellArray(x,2) = j_2;
                newCellArray(x,3) = angle;
            end
        end
        
        function emptyCells = emptyCellArray(stateLattice,coordinates)

            % emptyCellArray gives an array containing the indices of every empty
            % cell surrounding a cell at (i,j) with dimensions emptySpots x 2.
            latticeSizes = size(stateLattice);
            xlatticeSize = latticeSizes(2);
            ylatticeSize = latticeSizes(1);
            
            ind = 1; 
            
            i = coordinates(1);
            j = coordinates(2);
            
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
        
        function aliveCellsno = countLiveNeighbours(stateLattice,coordinates)
            latticeSizes = size(stateLattice);
            xlatticeSize = latticeSizes(2);
            ylatticeSize = latticeSizes(1);
            aliveCellsno = 0;
            i = coordinates(1);
            j = coordinates(2);
            for I = -1:1
                for J = -1:1
                    if I ~= 0 || J ~= 0
                        if i+I <= ylatticeSize && i+I > 0 && j+J <= xlatticeSize && j+J > 0
                            aliveCellsno = aliveCellsno + stateLattice(i+I,j+J);
                        end
                    end
                end
            end
        end
    end 
end
