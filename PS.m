function country = PS(country,F)
    nF=numel(F);
    Dmax = zeros(1, nF);
    Dmin = zeros(1, nF);
    for k=1:nF
        keep = zeros(1,1);
        for i = 1 : numel(F{k})
            for j = 1 : numel(F{k})
                if j ~= i                   
                    %A = sqrt((country(F{k}(1,i)).Cost(1) - country(F{k}(1,j)).Cost(1))^2 + (country(F{k}(1,i)).Cost(2) - country(F{k}(1,j)).Cost(2))^2 + (country(F{k}(1,i)).Cost(3) - country(F{k}(1,j)).Cost(3))^2);
                    A = sqrt((country(F{k}(1,i)).Cost(1) - country(F{k}(1,j)).Cost(1))^2 + (country(F{k}(1,i)).Cost(2) - country(F{k}(1,j)).Cost(2))^2);
                    keep = [keep; A]; %#ok<*AGROW>
                end
            end
        end
    if numel(keep)>1 && keep(1) == 0
        keep(1) = [];
    else        
    end
    Dmax(k) = max(keep);
    Dmin(k) = min(keep);
    end
        
    for k=1:nF
        for i = 1 : numel(F{k})
            NormalizedDistance = 0;
            Sum_Distance = 0;
            Mu = 0;
            for j = 1 : numel(F{k})
                if j ~= i                   
                    %NormalizedDistance = (sqrt((country(F{k}(1,i)).Cost(1) - country(F{k}(1,j)).Cost(1))^2 + (country(F{k}(1,i)).Cost(2) - country(F{k}(1,j)).Cost(2))^2 + (country(F{k}(1,i)).Cost(3) - country(F{k}(1,j)).Cost(3))^2)) - (Dmax(k)-Dmin(k));
                    NormalizedDistance = sqrt((country(F{k}(1,i)).Cost(1) - country(F{k}(1,j)).Cost(1))^2 + (country(F{k}(1,i)).Cost(2) - country(F{k}(1,j)).Cost(2))^2);
                    Sum_Distance = NormalizedDistance + Sum_Distance;
                end
            end
            Mu = Sum_Distance/(numel(F{k})-1);
            PenalizationTerm = 0;
            for j = 1 : numel(F{k})
                if j ~= i                   
                    B = (Dmax(k)-Dmin(k))/(sqrt((country(F{k}(1,i)).Cost(1) - country(F{k}(1,j)).Cost(1))^2 + (country(F{k}(1,i)).Cost(2) - country(F{k}(1,j)).Cost(2))^2));
                    PenalizationTerm = PenalizationTerm + B;
                end
            end
            PM = Mu + PenalizationTerm;
            country(F{k}(1,i)).PS = PM + (country(F{k}(1,i)).Rank-1)*numel(country(F{k}(1,i)).Cost);
            TT = isnan(country(F{k}(1,i)).PS);
            if TT == 1
                country(F{k}(1,i)).PS = max(Dmax) + (country(F{k}(1,i)).Rank-1)*numel(country(F{k}(1,i)).Cost);
            end
        end
    end
end