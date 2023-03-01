function out = values2struct(values,names)

for i=1:size(names,1)
   nn = names(i);
   out.(nn{1}) = values(i,1);
end

end

