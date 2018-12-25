function en = computeEnergy(f2var, factors, lab)

en = 0;

if strcmp(class(factors{1}), 'double')
    for k = 1:numel(factors)        
        t = factors{k};    
        if ~isempty(t)
            tmp = num2cell(lab(f2var{k}));
            val = t(sub2ind(size(t), tmp{:}));   
            en = en + -log(val);
        end
    end  

else % assume bayesNet structure

    for k = 1:numel(factors)    

        t = convert_to_table(factors{k}, f2var{k}, cell(numel(lab), 1));    
        tmp = num2cell(lab(f2var{k}));
        val = t(sub2ind(size(t), tmp{:}));

        en = en + -log(val);
    end
end    
    