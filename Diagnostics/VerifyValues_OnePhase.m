function [] = VerifyValues_OnePhase(FileName)
    
    
    %                                Input Checking
    % ================================================================================
    if (nargin > 0) && not(isempty(FileName)) && ischar(FileName)
        OutputFileName  = [CallerDirectory(),FileName];
        FileID          = fopen(OutputFileName,'w','native','UTF-8');
        
        if (FileID ~= -1)
            Print = @(Format,varargin) fprintf(FileID,Format,varargin{:});
        else
            fprintf('\n');
            warning('WWPP:VerifyTwoPhaseValues:UnableToOpenFile',...
                ['**** The requested file could not be created, ',...
                'so the verification is being output to the command window. ****']);
            fprintf('\n');
            Print = @(Format,varargin) fprintf(Format,varargin{:});
        end
    
    else
        Print = @(Format,varargin) fprintf(Format,varargin{:});
        FileID = -1;
    end

    
% Data load and calculations
    FileNameSinglePhase = 'IAPWSValues_OnePhase.csv';
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
    RelDiff_P  = AbsDiff_P  ./ IAPWS_P    ;
    RelDiff_cv = AbsDiff_cv ./ IAPWS_cv   ;
    RelDiff_w  = AbsDiff_w  ./ IAPWS_w    ;
    RelDiff_s  = AbsDiff_s  ./ IAPWS_s    ;
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
        Data{k} = [T(k),rho(k),AbsDiff_P(k),AbsDiff_cv(k),AbsDiff_w(k),AbsDiff_s(k)];
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
        Data{k} = [T(k),rho(k),RelDiff_P(k),RelDiff_cv(k),RelDiff_w(k),RelDiff_s(k)];
    end
    
    % Write table
    WriteTable(Title,DataLine,Data,BufferSize);
    

    if (FileID ~= -1)
        fclose(FileID);
    end
    
    
    
    
    
    
% Subfunctions
function [] = WriteTable(TitleLine,DataLine,Data,BufferSize)
    TopRule(BufferSize);
    WritePropertyHeaders(TitleLine);
    MiddleRule(BufferSize);
    for m = 1:length(Data)
        Print(DataLine,Data{m}(:));
    end
    BottomRule(BufferSize);
end


function [] = WritePropertyHeaders(Title)
    Print([GetCenteredString('State',28),repmat(' ',1,8)]);
    Print([GetCenteredString(Title,80),'\n']);
    Print([repmat('-',1,28),repmat(' ',1,8),repmat('-',1,80),'\n']);
    Print(['   T [K]      rho[kg/m^3]  ',repmat(' ',1,10),...
             '     P [Pa]           '               ,...
             ' Cv [J/kg-K]          '               ,...
             '  w [m/s]            '                ,...
             ' s[J/kg-K]\n'                         ] );
end

% ========================== Table rules ========================== 
    function [] = TopRule(BufferSize)
        Print([repmat('=',1,BufferSize),'\n']);
    end
    
    function [] = MiddleRule(BufferSize)
        Print([repmat('-',1,BufferSize),'\n']);
    end
    
    function [] = BottomRule(BufferSize)
        Print([repmat('=',1,BufferSize),'\n']);
    end

    function CenteredString = GetCenteredString(String,BufferSize)
        StringLength = length(String);
        LengthPadLeft  = ceil ((BufferSize-StringLength)/2);
        LengthPadRight = floor((BufferSize-StringLength)/2);
        
        CenteredString = [repmat(' ',1,LengthPadLeft)   ,...
                                   String               ,...
                          repmat(' ',1,LengthPadRight)  ];
    end
end