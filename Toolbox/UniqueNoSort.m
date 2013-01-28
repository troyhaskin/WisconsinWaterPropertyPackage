function [UniqueArray,UniqueMask] = UniqueNoSort(Array)
    
    [UniqueArray,UndoSort,UniqueMask] = unique(Array);
    UniqueArray = UniqueArray(UndoArray);
    
end