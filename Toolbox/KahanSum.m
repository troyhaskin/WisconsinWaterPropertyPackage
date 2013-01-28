function [GoodSum,Error] = KahanSum(Sum,Part,Error)
    Shift = Part - Error;
    GoodSum = Sum + Shift;
    Error        = (GoodSum - Sum) - Shift;
%     
%     GoodSum = Sum + Part;
end