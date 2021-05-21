function [Population,FrontNo,Fitness,D,Nadir] = EnvironmentalSelection(Population,N,IdealPoint)
    % The environmental selection of C-CLOP

    % This function is written by Ben Crulis

    %% let's round the objective values
    objs = round(Population.objs, 6);
    %objs = Population.objs;
    
    %% Non-dominated sorting
    [FrontNo, Last] = NDSort(objs,Population.cons,N);

    Front = FrontNo == 1;

    [nInd, ~] = size(objs);
    
    Infeas = any(Population.cons > 0, 2);
    ConsF = (Population.cons > 0) .* Population.cons;
    Cons = sum(ConsF, 2);
        
    %% Compute the normalized scores
    F1 = objs(Front,:);
    
    Nadir = max(F1, [], 1);
    
    Normalization = Nadir - IdealPoint;
    Normalization(Normalization == 0) = 1;

    F = (objs - repmat(IdealPoint, nInd,1))./repmat(Normalization, nInd,1);
        
    %% compute the fitness as the distance to the ideal point

    Fitness = squeeze(sum(F.^2,2));
    
    %% compute the distance between solutions
        
    Fu = F - 1e-6;
    
    D = pdist2(Fu,Fu,'cosine');
    
    D(isnan(D)) = 0;
    
    D(logical(eye(nInd))) = inf;
    
    %% Truncate solutions
    Keep = false(1,nInd);
    Keep(FrontNo <= Last) = true;
    
    % create variables for the subpopulation
    subsetD = D(Keep, Keep);
    subsetFitness = Fitness(Keep);
    KeepSubset = logical(ones(1, size(subsetD, 1)));
    
    KeepSubset = Closest_selection(KeepSubset, N, subsetD, subsetFitness);
    
    Keep(Keep == true) = KeepSubset;
        
    %% Population for next generation
    Population = Population(Keep);
    FrontNo    = FrontNo(Keep);
    Fitness    = Fitness(Keep);
    D = D(Keep,Keep);
end