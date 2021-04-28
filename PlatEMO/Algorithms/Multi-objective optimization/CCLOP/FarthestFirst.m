function index = FarthestFirst(K,N,D,rank,fitness)
    FrontInd = rank == 1;

    FrontD = D(FrontInd, FrontInd);
    
    Dist = sum(mink(FrontD, K, 2), 2);
    
    %Dist(isinf(Dist)) = 1;
    %Dist(isnan(Dist)) = 1;
    
    SumDist = sum(Dist);
    
    NormDist = Dist ./ repmat(SumDist, size(Dist));

    NbC = ceil(NormDist * N);
    
    Parents = repelem(find(FrontInd), NbC);
    
    index = Parents(randperm(size(Parents,2), N));
end