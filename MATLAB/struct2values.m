function valuesArray = struct2values(structEx,names)
if nargin > 1
    valuesArray = zeros(size(names,1),1);
    for i=1:size(names,1)
        name = names(i,1);
        valuesArray(i,1) = getfield(structEx,name{1});
    end
else
    namesStruct = fieldnames(structEx);
    valuesArray = zeros(size(namesStruct,1),1);
    for i=1:size(namesStruct,1)
        name = namesStruct(i,1);
        valuesArray(i,1) = getfield(structEx,name{1});
    end
end
