% This function computes the Shannon and Simpson indices

function [shannon,simpson] = diversityMeasures(cell_props)

    N_types = length(cell_props); % The number of species is the number of distinct genotypes of cells considered
    shannon = 0;
    simpson = 0;

    for i=1:N_types
        if(cell_props(i)~=0)
            shannon = shannon + (cell_props(i) * log(cell_props(i)));
        end
        simpson = simpson + (cell_props(i) * cell_props(i));
    end
    shannon = shannon * (-1);
    simpson = 1 - simpson; % As stated in the text, we consider 1-Simpson index, to ensure comparabiliity with other indices

end


