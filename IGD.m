function IGDValue = IGD(truePF, PF)

% The PopObj denotes the population obtained from MOEAs
% PopObj" always denotes the matrix of the objective values of the population
% trueFront is the ideal Pareto front of the benchmark problem

% <metric> <min>
% Inverted generational distance
%PF is a Pareto front ontained by a MOEA, and truePF is the ideal Pareto front of a MOP

IGD=0;
%A=repmat(truePF,size(PF,1),1);
    for i=1:size(PF,2)
        diff=truePF(:,i)-PF(:,i);
        dist=sqrt(sum(diff.^2,2));         
        IGD=IGD+min(dist);
    end
IGDValue=IGD/size(PF,2);
end