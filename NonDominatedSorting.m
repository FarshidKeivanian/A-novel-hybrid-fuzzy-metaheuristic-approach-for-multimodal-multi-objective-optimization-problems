function [country, F]=NonDominatedSorting(country)

    nPop=numel(country);

    for i=1:nPop
        country(i).DominationSet = [];
        country(i).DominatedCount = 0;
    end
    
    F{1}=[];
    
    for i=1:nPop
        for j=i+1:nPop
            p=country(i);
            q=country(j);
            
            if Dominates(p,q)
                p.DominationSet=[p.DominationSet j];
                q.DominatedCount=q.DominatedCount+1;
            end
            
            if Dominates(q.Cost,p.Cost)
                q.DominationSet=[q.DominationSet i];
                p.DominatedCount=p.DominatedCount+1;
            end
            
            country(i)=p;
            country(j)=q;
        end
        
        if country(i).DominatedCount==0
            F{1}=[F{1} i];
            country(i).Rank=1;
        end
    end
    
    k=1;
    
    while true
        
        Q=[];
        
        for i=F{k}
            p=country(i);
            
            for j=p.DominationSet
                q=country(j);
                
                q.DominatedCount = q.DominatedCount-1;
                
                if q.DominatedCount==0
                    Q=[Q j];
                    q.Rank=k+1;
                end
                
                country(j)=q;
            end
        end
        
        if isempty(Q)
            break;
        end
        
        F{k+1}=Q;
        
        k=k+1;
        
    end
    

end