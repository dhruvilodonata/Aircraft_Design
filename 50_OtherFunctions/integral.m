function result = integral(fx,x)
if length(fx) ~= length(x)
    result = NaN;
    return
end
result = zeros(length(x),1);
for i = 1:length(x)-1
    index = length(x) - i + 1;
    fx_i = (fx(index) + fx(index-1)) / 2;
    dx_i = x(index) - x(index - 1);
    result(index-1) = result(index) + fx_i * dx_i;
end
end