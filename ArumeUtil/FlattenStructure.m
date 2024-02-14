function s = FlattenStructure(s)
    if (~isstruct(s))
        error('parameter s should be a struct');
    end
    fields = fieldnames(s);
    
    for i=1:length(fields)
        s2 = s.(fields{i});
        
        if ( isstruct(s2))
            
            s2 = FlattenStructure(s2);
            
            fields2 = fieldnames(s2);
            
            for i2 = 1:length(fields2)
                s.([fields{i} '__' fields2{i2}]) = s2.(fields2{i2});
            end
            
            s = rmfield(s, fields{i});
        end
    end
end