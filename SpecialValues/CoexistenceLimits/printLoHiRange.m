function [] = printLoHiRange(range,name)
    fprintf(['Smallest ',name,': %+23.16E\n'],min(range));
    fprintf(['Largest  ',name,': %+23.16E\n'],max(range));
end