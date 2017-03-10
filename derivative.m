function [hz, ht] = derivative(fctHandle, x0, step)

nz = length(x0);
hz = zeros(nz);

h0 = feval(fctHandle,x0);

for i = 1:nz
    
end

end