function [] = VerifySinglePhaseValues()
    
%% Data load and calculations
    FileNameSinglePhase = 'IAPWSCheckValues_OnePhase.csv';
    Data = importdata(FileNameSinglePhase, ',', 1);
    
    % Number of data sets
    Nsets = size(Data.data,1);
    
    % State variables
    T   = Data.data(:,1);
    rho = Data.data(:,2);
    
    % IAPWS Values
    IAPWS_P  = Data.data(:,3);
    IAPWS_cv = Data.data(:,4); 
    IAPWS_w  = Data.data(:,5);
    IAPWS_s  = Data.data(:,6);
    
    % Package Values
%     WWPP_P  = PressureOne(rho,T);
%     WWPP_cv = HeatCapacityIsochoricOne(rho,T);
%     WWPP_w  = SoundSpeedOne(rho,T);
%     WWPP_s  = EntropyOne(rho,T);
    WWPP_P  = SigFigRound(PressureOne(rho,T),9);
    WWPP_cv = SigFigRound(HeatCapacityIsochoricOne(rho,T),9);
    WWPP_w  = SigFigRound(SoundSpeedOne(rho,T),9);
    WWPP_s  = SigFigRound(EntropyOne(rho,T),9);
    
    % Absolute Differences
    AbsDiffP  = IAPWS_P - WWPP_P  ;
    AbsDiffcv = IAPWS_cv- WWPP_cv ;
    AbsDiffw  = IAPWS_w - WWPP_w  ;
    AbsDiffs  = IAPWS_s - WWPP_s  ;
    
    % Relative Differences (w.r.t to IAPWS)
    RelDiffP  = AbsDiffP  ./ IAPWS_P    ;
    RelDiffcv = AbsDiffcv ./ IAPWS_cv   ;
    RelDiffw  = AbsDiffw  ./ IAPWS_w    ;
    RelDiffs  = AbsDiffs  ./ IAPWS_s    ;
    
%% Table 1: IAPWS values

    % Buffer size
    BufferSize = 120;
    
    % title
    Title     = {'IAPWS Verification Values'};
    OuterPad  = round((BufferSize - length(Title{1}))/2);
    TitleLine = MakeTableRow(Title,BufferSize,0,OuterPad);
    
    % header
    HeaderTitles = {'T [K]','rho[kg/m^3]','P [Pa]','Cv [J/kg-K]','w [m/s]','s[J/kg-K]'};
    InnerPad   = [5,15,14,12,13];
    HeaderLine = MakeTableRow(HeaderTitles,BufferSize,InnerPad,3);
    
    % data
    DataFormats = [{'%+G'},repmat({'%+13.5E'},1,5)];
    InnerPad    = [5,10,9,9,8];
    DataLine    = MakeTableRow(DataFormats,BufferSize,InnerPad,3);
    
    % Manipulate data for printing
    Data = cell(Nsets,1);
    for k = 1:Nsets
        Data{k} = [T(k),rho(k),IAPWS_P(k),IAPWS_cv(k),IAPWS_w(k),IAPWS_s(k)];
    end
    
    % Write table
    WriteTable(TitleLine,HeaderLine,DataLine,Data,BufferSize);
    fprintf('\n\n');
    
    
%% Table 2: WWPP values ===================
    % title
    Title     = {'WWPP Verification Values'};
    OuterPad  = round((BufferSize - length(Title{1}))/2);
    TitleLine = MakeTableRow(Title,BufferSize,0,OuterPad);
    
    % header
    HeaderTitles = {'T [K]','rho[kg/m^3]','P [Pa]','Cv [J/kg-K]','w [m/s]','s[J/kg-K]'};
    InnerPad   = [5,15,14,12,13];
    HeaderLine = MakeTableRow(HeaderTitles,BufferSize,InnerPad,3);
    
    % data
    DataFormats = [{'%+G'},repmat({'%+13.5E'},1,5)];
    InnerPad    = [5,10,9,9,8];
    DataLine    = MakeTableRow(DataFormats,BufferSize,InnerPad,3);
    
    % Manipulate data for printing
    Data = cell(Nsets,1);
    for k = 1:Nsets
        Data{k} = [T(k),rho(k),WWPP_P(k),WWPP_cv(k),WWPP_w(k),WWPP_s(k)];
    end
    
    % Write table
    WriteTable(TitleLine,HeaderLine,DataLine,Data,BufferSize);
  
    
%% Table 3: Absolute Differences ===================
    % title
    Title     = {'Absolute Differences in Verification Values'};
    OuterPad  = round((BufferSize - length(Title{1}))/2);
    TitleLine = MakeTableRow(Title,BufferSize,0,OuterPad);
    
    % header
    HeaderTitles = {'T [K]','rho[kg/m^3]','P [Pa]','Cv [J/kg-K]','w [m/s]','s[J/kg-K]'};
    InnerPad   = [5,15,14,12,13];
    HeaderLine = MakeTableRow(HeaderTitles,BufferSize,InnerPad,3);
    
    % data
    DataFormats = [{'%+G'},repmat({'%+13.5E'},1,5)];
    InnerPad    = [5,10,9,9,8];
    DataLine    = MakeTableRow(DataFormats,BufferSize,InnerPad,3);
    
    % Manipulate data for printing
    Data = cell(Nsets,1);
    for k = 1:Nsets
        Data{k} = [T(k),rho(k),AbsDiffP(k),AbsDiffcv(k),AbsDiffw(k),AbsDiffs(k)];
    end
    
    % Write table
    WriteTable(TitleLine,HeaderLine,DataLine,Data,BufferSize);
    
    
%% Table 4: Relative Differences ===================
    % title
    Title     = {'Relative Differences in Verification Values'};
    OuterPad  = round((BufferSize - length(Title{1}))/2);
    TitleLine = MakeTableRow(Title,BufferSize,0,OuterPad);
    
    % header
    HeaderTitles = {'T [K]','rho[kg/m^3]','P [Pa]','Cv [J/kg-K]','w [m/s]','s[J/kg-K]'};
    InnerPad   = [5,15,14,12,13];
    HeaderLine = MakeTableRow(HeaderTitles,BufferSize,InnerPad,3);
    
    % data
    DataFormats = [{'%+G'},repmat({'%+13.5E'},1,5)];
    InnerPad    = [5,10,9,9,8];
    DataLine    = MakeTableRow(DataFormats,BufferSize,InnerPad,3);
    
    % Manipulate data for printing
    Data = cell(Nsets,1);
    for k = 1:Nsets
        Data{k} = [T(k),rho(k),RelDiffP(k),RelDiffcv(k),RelDiffw(k),RelDiffs(k)];
    end
    
    % Write table
    WriteTable(TitleLine,HeaderLine,DataLine,Data,BufferSize);
    
end

%% Subfunctions
function [] = WriteTable(TitleLine,HeaderLine,DataLine,Data,BufferSize)
    fprintf(TitleLine);
    TopRule(BufferSize);
    fprintf(HeaderLine);
    MiddleRule(BufferSize);
    for k = 1:length(Data)
        fprintf(DataLine,Data{k}(:));
    end
    BottomRule(BufferSize);
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

% ======================== Make Table Row =======================
    function TableRow = MakeTableRow(Data, BufferSize,InnerPad,OuterPad)
        TableRow = repmat(' ',1,BufferSize);
        
        
        Ntitles = length(Data);
        Pad     = [OuterPad,InnerPad.*ones(1,Ntitles-1),OuterPad];
        Top     = 0;
        for k = 1 : length(Data)
            Datum  = Data{k};
            Bottom = Top + Pad(k) + 1;
            Top    = Bottom + length(Datum) - 1;
            TableRow(Bottom:Top) = Datum;
        end
        TableRow = [TableRow,'\n'];
    end
