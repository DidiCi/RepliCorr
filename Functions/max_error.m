function b=max_error(a)
    b=(max(a)-min(a))/2;
    if isempty(a) b=NaN;
end