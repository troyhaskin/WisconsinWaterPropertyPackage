function [] = VerifySinglePhaseValues()
    
% Data load and calculations
    FileNameSinglePhase = 'IAPWSCheckValues_OnePhase.csv';
    Data = importdata(FileNameSinglePhase, ',', 1);
    
    % Number of data sets
    Nsets = size(Data.data,1);
    
    % State variables
    T   = Data.data(:,1);
    rho = Data.data(:,2);
    
    % IAPWS Values
    IAPWS_P  = Data.data(:,3)   ;
    IAPWS_cv = Data.data(:,4)   ;
    IAPWS_w  = Data.data(:,5)   ;
    IAPWS_s  = Data.data(:,6)   ;
    IAPWS    = Data.data        ;
    
    % Package Values
    WWPP_P  = SigFigRound(PressureOne(rho,T),9);
    WWPP_cv = SigFigRound(HeatCapacityIsochoricOne(rho,T),9);
    WWPP_w  = SigFigRound(SoundSpeedOne(rho,T),9);
    WWPP_s  = SigFigRound(EntropyOne(rho,T),9);
    WWPP    = [T,rho,WWPP_P,WWPP_cv,WWPP_w,WWPP_s];
    
    % Absolute Differences
    AbsDiff_P  = IAPWS_P - WWPP_P  ;
    AbsDiff_cv = IAPWS_cv- WWPP_cv ;
    AbsDiff_w  = IAPWS_w - WWPP_w  ;
    AbsDiff_s  = IAPWS_s - WWPP_s  ;
    AbsDiff    = [T,rho,AbsDiff_P,AbsDiff_cv,AbsDiff_w,AbsDiff_s];
    
    % Relative Differences (w.r.t to IAPWS)
    RelDiff_P  = AbsDiffP  ./ IAPWS_P    ;
    RelDiff_cv = AbsDiffcv ./ IAPWS_cv   ;
    RelDiff_w  = AbsDiffw  ./ IAPWS_w    ;
    RelDiff_s  = AbsDiffs  ./ IAPWS_s    ;
    RelDiff    = [T,rho,RelDiff_P,RelDiff_cv,RelDiff_w,RelDiff_s];
    
% Table 1: IAPWS values

    % Buffer size
    BufferSize = 120;
    
    % title
    Title     = 'IAPWS Verification Values';
    
    % data
       DataLine = ['   %+4G','     ','%+13.8E','    ',...
                   repmat('      %+13.8E',1,4),'\n'];
    
    % Manipulate data for printing
    Data = cell(Nsets,1);
    for k = 1:Nsets
        Data{k} = [T(k),rho(k),IAPWS_P(k),IAPWS_cv(k),IAPWS_w(k),IAPWS_s(k)];
    end
    
    % Write table
    WriteTable(Title,DataLine,Data,BufferSize);
    
    
% Table 2: WWPP values ===================
    % title
    Title     = 'WWPP Calculated Values';
    
    % data
       DataLine = ['   %+4G','     ','%+13.8E','    ',...
                   repmat('      %+13.8E',1,4),'\n'];
    
    % Manipulate data for printing
    Data = cell(Nsets,1);
    for k = 1:Nsets
        Data{k} = [T(k),rho(k),WWPP_P(k),WWPP_cv(k),WWPP_w(k),WWPP_s(k)];
    end
    
    % Write table
    WriteTable(Title,DataLine,Data,BufferSize);
  
    
% Table 3: Absolute Differences ===================
    % title
    Title     = 'Absolute Differences in Verification Values';


    % data
       DataLine = ['   %+4G','     ','%+13.8E','    ',...
                   repmat('      %+13.8E',1,4),'\n'];
    
    % Manipulate data for printing
    Data = cell(Nsets,1);
    for k = 1:Nsets
        Data{k} = [T(k),rho(k),AbsDiffP(k),AbsDiffcv(k),AbsDiffw(k),AbsDiffs(k)];
    end
    
    % Write table
    WriteTable(Title,DataLine,Data,BufferSize);
    
    
% Table 4: Relative Differences ===================
    % title
    Title     = 'Relative Differences in Verification Values';
    
    % data
       DataLine = ['   %+4G','     ','%+13.8E','    ',...
                   repmat('      %+13.8E',1,4),'\n'];
    
    % Manipulate data for printing
    Data = cell(Nsets,1);
    for k = 1:Nsets
        Data{k} = [T(k),rho(k),RelDiffP(k),RelDiffcv(k),RelDiffw(k),RelDiffs(k)];
    end
    
    % Write table
    WriteTable(Title,DataLine,Data,BufferSize);
    
end

% Subfunctions
function [] = WriteTable(TitleLine,DataLine,Data,BufferSize)
    TopRule(BufferSize);
    WritePropertyHeaders(TitleLine);
    MiddleRule(BufferSize);
    for k = 1:length(Data)
        fprintf(DataLine,Data{k}(:));
    end
    BottomRule(BufferSize);
end


function [] = WritePropertyHeaders(Title)
    fprintf([GetCenteredString('State',28),repmat(' ',1,8)]);
    fprintf([GetCenteredString(Title,80),'\n']);
    fprintf([repmat('-',1,28),repmat(' ',1,8),repmat('-',1,80),'\n']);
    fprintf(['   T [K]      rho[kg/m^3]  ',repmat(' ',1,10),...
             '     P [Pa]           '               ,...
             ' Cv [J/kg-K]          '               ,...
             '  w [m/s]            '                ,...
             ' s[J/kg-K]\n'                         ] );
end

% ========================== Table rules ========================== 
    function [] = TopRule(BufferSize)
        fprintf([repmat('=',1,BufferSize),'\n']);
    end
    
    function [] = MiddleRule(BufferSize)
        fprintf([repmat('-',1,BufferSize),'\n']);
    end
    
    function [] = BottomRule(BufferSize)
        fprintf([repmat('=',1,BufferSize),'\n']);
    end

    function CenteredString = GetCenteredString(String,BufferSize)
        StringLength = length(String);
        LengthPadLeft  = ceil ((BufferSize-StringLength)/2);
        LengthPadRight = floor((BufferSize-StringLength)/2);
        
        CenteredString = [repmat(' ',1,LengthPadLeft)   ,...
                                   String               ,...
                          repmat(' ',1,LengthPadRight)  ];
    end
