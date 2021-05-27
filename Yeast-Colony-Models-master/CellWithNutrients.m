%% Class for Cell (with nutrients)
classdef CellWithNutrients
    
    properties
        state = 0   % indicates whether a cell is located at that site.
        budScar = [0,0]   % gives a vector in the direction of the bud scar of a cell at that lattice site.
        nutrientConcentration = 5   % gives the number of nutrient packets at that lattice site.
        consumedNutrients = 0   % indicates how many nutrients the cell has consumed since budding last.
        unipolarFrac = 0.8   % indicates the fraction of the time the cell will bud regularly.
        deltaUnipolarFrac = 0.1   % indicates by how much unipolarFrac will change as the number of nutrient packets changes.
        buddable = true;   % indicates whether the cell is able to bud at this moment. 
        mutated = false;   % bool determining whether the cell has mutated.
        budded = false;   % bool indicates whether the cell has previously budded.
    end
    methods
        
        function obj = CellWithNutrients(state,budScar,nutrientConcentration,deltaUnipolarFrac)
           % class constructor
             obj.state = state;
             obj.budScar = budScar;
             obj.nutrientConcentration = nutrientConcentration;
             obj.deltaUnipolarFrac = deltaUnipolarFrac;
        end
        
        function obj = updateNutrients(obj)
            % updateNutrients updates the nutrientConcentration and
            % consumedNutrients based on the cell consuming one nutrient
            
            if obj.state == 1 && obj.nutrientConcentration > 0
                obj.nutrientConcentration = obj.nutrientConcentration - 1;
                obj.consumedNutrients = obj.consumedNutrients + 1;
            end
        end
        
        function obj = updateUnipolarFrac(obj)
            % updateUnipolarFrac updates unipolarFrac based on
            % nutrientConcentration
            obj.unipolarFrac = obj.deltaUnipolarFrac * obj.nutrientConcentration;
        end
        
        function obj = updateBuddability(obj,nutrForBudding)
            % updateBuddability updates whether the cell can bud or not
            if isequal(obj.budScar,[0 0])
                obj.buddable = false;
            elseif obj.consumedNutrients < nutrForBudding
                obj.buddable = false;
            else
                obj.buddable = true;
            end
        end
    end
end