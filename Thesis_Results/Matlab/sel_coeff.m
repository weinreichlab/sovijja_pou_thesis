function selection = sel_coeff(wi, wj)
    
    selection = (wj - wi)/wi;
end

% Calculates the selection coefficient from allele i to allele j. 
% Input = The Wrightian fitness of allele i and allele j. 
% Output = selection coefficient (i -> j)
