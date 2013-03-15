function [] = VerifyOnePhaseValues(FileName)
    
    
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
    
    
    
    %                              Load, Allocate, and Assign
    % ================================================================================
    
    % CSV Import
    FileName = [CurrentFileDirectory(),'IAPWSCheckValues_OnePhase.csv'];
    Data     = importdata(FileName, ',', 1);
    
    % Number of data sets
    Nsets  = size(Data.data,1)      ;
    Nprops = size(Data.data,2) - 1  ;
    
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




    %                                Printing Setup
    % ================================================================================
    
    % Table (horizontal) padding
    PaddingLeft    = 3;
    PaddingInner   = 6;
    PaddingRight   = 3;
    OverHangHeader = 1;

    % DataWidth
    DataWidthInteger = 4;
    DataWidthReal    = 15;  %   Corresponds to scientific notation with 8 decimal places,
                            %   2 exponential places, and persistent sign indication

    WidthState = 2*OverHangHeader + DataWidthInteger + PaddingInner + DataWidthReal   ;
    WidthData  = 2*OverHangHeader + 4*DataWidthReal  + 3*PaddingInner                 ;

    % BufferSize
    BufferSize = PaddingLeft + WidthState + PaddingInner + WidthData + PaddingRight;


    % Pads
    PadLeft       = RepeatedString(' ',PaddingLeft );
    PadIn         = RepeatedString(' ',PaddingInner);
    PadRight      = RepeatedString(' ',PaddingRight);
    OverHangSpace = RepeatedString(' ',OverHangHeader);
    OverHangRule  = RepeatedString('-',OverHangHeader);
    
    
    % Helpful handles for clean/manageable code
    PrintTitle       = @(Title) Print('%s\n',CenteredString(Title,BufferSize));
    TopRule          = @(N) RepeatedString('=',N);  
    MidRule          = @(N) RepeatedString('-',N); 
    BottomRule       = @(N) TopRule(N)           ; 
    PrintTopRule     = @(N) Print('%s',TopRule   (N));
    PrintMidRule     = @(N) Print('%s',MidRule   (N));
    PrintBottomRule  = @(N) Print('%s',BottomRule(N));
    NewLine          = @()  Print('\n',[]);


    % Header handles
    HeaderLine1 = @() Print('%s\n%s\n%s\n',...
                    [   PadLeft                                             ,...
                        CenteredString('State',WidthState)                  ,...
                        PadIn                                               ,...
                        RepeatedString(' ',WidthData)                       ,...
                        PadRight                                            ],...
                    [   PadLeft                                             ,...
                        RepeatedString('-',WidthState)                      ,...
                        PadIn                                               ,...
                        OverHangSpace                                       ,...
                        CenteredString('P [Pa]'    ,DataWidthReal),PadIn    ,...
                        CenteredString('cv [Pa]'   ,DataWidthReal),PadIn    ,...
                        CenteredString('w [m/s]'   ,DataWidthReal),PadIn    ,...
                        CenteredString('s [J/kg-K]',DataWidthReal)          ]);
    
    HeaderLine2 = @() Print('%s\n',[PadLeft,OverHangSpace                               ,...
                                    CenteredString('Property',SatPropertyColumnWidth)   ,...
                                    OverHangSpace,PadIn,OverHangRule                    ,...
                                    MidRule(SatTempColumnWidth)                         ,...
                                    OverHangRule,PadRight]                              );
    
    HeaderLine3 = @() Print('%s\n',[PadLeft,OverHangSpace                       ,...
                                    RepeatedString(' ',SatPropertyColumnWidth)  ,...
                                    OverHangSpace,PadIn,OverHangSpace           ,...
                                    SatTempHeaderString                         ,...
                                    OverHangSpace,PadRight]                     );




    %                          Table 1: IAPWS values
    % ================================================================================
    
    PrintTitle('Table 1: IAPWS Reference Values for Saturated Water');
    PrintTopRule(BufferSize);NewLine();     %   Toprule
    HeaderLine1();                          %       Header
    HeaderLine2();                          %       Header
    HeaderLine3();                          %       Header
    PrintMidRule(BufferSize);NewLine();     %   MidRule
    PrintDataBlock(IAPWS');                 %       Data
    PrintBottomRule(BufferSize);NewLine();  %   BottomRule
    NewLine();                              %   Table Separation




    %                             Table 2: WWPP values
    % ================================================================================
    
    PrintTitle('Table 2: WWPP Calculated Values for Saturated Water');
    PrintTopRule(BufferSize);NewLine();
    HeaderLine1();
    HeaderLine2();
    HeaderLine3();
    PrintMidRule(BufferSize);NewLine();
    PrintDataBlock(WWPP');
    PrintBottomRule(BufferSize);NewLine();
    NewLine();




	%                   Table 3: Absolute Difference values
    % ================================================================================
    
    PrintTitle('Table 3: IAPWS-WWPP Absolute Differences for Saturated Water');
    PrintTopRule(BufferSize);NewLine();
    HeaderLine1();
    HeaderLine2();
    HeaderLine3();
    PrintMidRule(BufferSize);NewLine();
    PrintDataBlock(AbsDiff');
    PrintBottomRule(BufferSize);NewLine();
    NewLine();




    %                      Table 4: Absolute Difference values
    % ================================================================================
    
    PrintTitle('Table 4: IAPWS-WWPP Relative Differences for Saturated Water');
    PrintTopRule(BufferSize);NewLine();
    HeaderLine1();
    HeaderLine2();
    HeaderLine3();
    PrintMidRule(BufferSize);NewLine();
    PrintDataBlock(RelDiff')
    PrintBottomRule(BufferSize);NewLine();

    NewLine();


    % Close the file if one is open.
    if (FileID ~= -1)
        fclose(FileID);
    end
    
    
    
    
    
    function [] = PrintDataBlock(DataToPrint)
        for n = 1:Nprops
            Print('%s',PadLeft);
            Print('%s',[OverHangSpace,SatPropertyColumn(n,:),OverHangSpace]);
            Print('%s',PadIn);
            Print('%s',OverHangSpace);
            for p = 1:Nsets-1
                Print('%+15.8E',DataToPrint(n,p));
                Print('%s'     ,PadIn     );
            end
            Print('%+15.8E',DataToPrint(n,Nsets));
            Print('%s',OverHangSpace);
            Print('%s\n',PadRight);
        end
    end
    
end

function TheRepeatedString = RepeatedString(String,Length)
    TheRepeatedString = repmat(String,1,Length);
end

function TheCenteredString = CenteredString(String,BufferSize)
    StringLength = length(String);
    LengthPadLeft  = ceil ((BufferSize-StringLength)/2);
    LengthPadRight = floor((BufferSize-StringLength)/2);
    
    TheCenteredString = [repmat(' ',1,LengthPadLeft)   ,...
        String               ,...
        repmat(' ',1,LengthPadRight)  ];
end
