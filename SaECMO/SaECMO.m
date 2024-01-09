classdef SaECMO < ALGORITHM  
    % <multi> <real/binary/permutation> <constrained>
    %SaECMO
%------------------------------- Reference --------------------------------
% S. Song, K. Zhang, L. Zhang, and N. Wu, “A dual-population algorithm based 
% on self-adaptive epsilon method for constrained multi-objective optimization,”
% Information Sciences, vol. 655, pp. 119906–119906, Jan. 2024, 
% doi: https://doi.org/10.1016/j.ins.2023.119906.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population1 = Problem.Initialization();
            Fitness1    = CalFitness(Population1.objs,Population1.cons);
            Population2 = Problem.Initialization();
            Fitness2 = CalFitness(Population2.objs);
            cons = [Population1.cons;Population2.cons];
            cons(cons<0) = 0;
            CV0 = max(sum(cons,2));
            if CV0 == 0
                CV0 = 1;
            end
            X=0;
            cp=(-log(CV0)-6)/log(1-0.5);
            N = Problem.N;  
         
            %% Optimization
            while Algorithm.NotTerminated(Population1)
                VAR = 1 - (1-X)^cp;
                MatingPool1 = TournamentSelection(2,2*N,Fitness1);
                Offspring1  = OperatorDE(Problem,Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                MatingPool2 = TournamentSelection(2,2*N,Fitness2);
                Offspring2 = OperatorDE(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
                %% Environmental selection
                [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring2,Offspring1], N, VAR, false);
                [Population1,Fitness1] = EnvironmentalSelection([Population1, Offspring1,Offspring2],N, VAR, true);                                 
                X=X+1/(Problem.maxFE/Problem.N);
            end          
        end
    end
end