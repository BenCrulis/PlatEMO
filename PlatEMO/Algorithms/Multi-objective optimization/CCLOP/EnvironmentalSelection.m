function [Population,FrontNo,Fitness,D,Nadir] = EnvironmentalSelection(Population,N,IdealPoint)
% The environmental selection of C-CLOP


% This function is written by Ben Crulis

    %% let's round the objective values
    objs = round(Population.objs, 6);
    %objs = Population.objs;
    
    %% Non-dominated sorting
    [FrontNo, Last] = NDSort(objs,Population.cons,N);
    %[FrontNo, Last] = NDSort(objs, N);
    Front = FrontNo == 1;

    [nInd, ~]    = size(objs);
    
    Infeas = any(Population.cons > 0, 2);
    ConsF = (Population.cons > 0) .* Population.cons;
    Cons = sum(ConsF, 2);
        
    %% Compute the normalized scores
    F1 = objs(Front,:);
    
    Nadir = max(F1, [], 1);
    %Nadir = mean(F1, 1); % using max is unstable for some problems such as IDTLZ1
    %Nadir = median(F1, 1);
    
    Normalization = Nadir - IdealPoint;
    Normalization(Normalization == 0) = 1;

    F = (objs - repmat(IdealPoint, nInd,1))./repmat(Normalization, nInd,1);
        
    %% compute the fitness as the distance to the ideal point

    Fitness = squeeze(sum(F.^2,2));
    %Fitness = sum(F, 2);
    
    %% compute the distance between solutions
    
    %F = WS(F + 1e-7); % transform the vectors using MOEA/D-AWA WS transformation
    
    Fu = F - 1e-6;
    
    D = pdist2(Fu,Fu,'cosine');
    
    D(isnan(D)) = 0;
    
    D(logical(eye(nInd))) = inf;
    
    %% Truncate solutions
    Keep = false(1,nInd);
    Keep(FrontNo <= Last) = true;
    
    %Keep(Infeas) = false;
    
    %if sum(Keep) < N
        %Keep(~Infeas) = true;
    %end
    
    Corner = FindCornerSolutions(F);
    %Corner = [];
    subsetD(Corner) = 1;
    %Fitness(Corner) = -inf;
    
    % create variables for the subproblem
    subsetD = D(Keep, Keep);
    subsetFitness = Fitness(Keep);
    KeepSubset = logical(ones(1, size(subsetD, 1)));
    
    KeepSubset = Closest_selection(KeepSubset, N, subsetD, subsetFitness);
    
    Keep(Keep == true) = KeepSubset;
        
    %% Population for next generation
    Population = Population(Keep);
    FrontNo    = FrontNo(Keep);
    %Fitness = sum(mink(D, size(objs,2), 2), 2); % use angle fitness
    Fitness    = Fitness(Keep);
    D = D(Keep,Keep);
end