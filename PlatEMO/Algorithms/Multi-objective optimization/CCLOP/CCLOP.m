classdef CCLOP < ALGORITHM
% <multi/many> <real/binary/permutation> <constrained/none>
% Cosinus pair based multiobjective Genetic Algorithm

%------------------------------- Reference --------------------------------
%------------------------------- Copyright --------------------------------
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            
            %% Generate random population
            Population = Problem.Initialization();
            
            IdealPoint = min(Population.objs, [], 1);
            OldNadir = max(Population.objs, [], 1);

            [~,FrontNo,Fitness,~,~] = EnvironmentalSelection(Population,Problem.N,IdealPoint);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
              MatingPool = TournamentSelection(2,Problem.N,FrontNo,Fitness);
              Offspring  = OperatorGA(Population(MatingPool));
              IdealPoint = min(cat(1, IdealPoint, Offspring.objs), [], 1);
              [Population,FrontNo,Fitness,~,~] = EnvironmentalSelection([Population,Offspring],Problem.N,IdealPoint);
            end
        end
    end
end