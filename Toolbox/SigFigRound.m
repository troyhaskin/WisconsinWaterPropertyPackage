function Rounded = SigFigRound(Array,NumberOfSigFigs)
    
    if (nargin < 2)
        NumberOfSigFigs = 0;
    elseif (NumberOfSigFigs < 0) || not(IsInteger(NumberOfSigFigs))
        error('Thermodynamics:SigFigRound:BadSigFigInput',...
              'The NumberOfSigFigs argument must be a positive integer.');
    end
    
%     Scale       = 10^NumberOfSigFigs            ;
%     ScaledArray = Scale * Array                 ;
%     Rounded     = fix(ScaledArray) / Scale    ;
    Rounded = str2num(num2str(Array,['%0.',num2str(NumberOfSigFigs-1),'E'])); %#ok<ST2NM>
    
end

function TrueFalse = IsInteger(Number)
    TrueFalse = eq(round(Number),Number);
end