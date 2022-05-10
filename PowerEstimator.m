function country = PowerEstimator(country)

    for i=1:numel(country)
        country(i).Power = 1/(country(i).PS);
    end
    
end