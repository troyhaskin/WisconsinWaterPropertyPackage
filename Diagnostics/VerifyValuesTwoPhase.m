function [] = VerifyValuesTwoPhase(FileName)
    
    
    %                                Input Checking
    % ================================================================================
    if (nargin > 0) && not(isempty(FileName)) && ischar(FileName)
        OutputFileName  = [CallerDirectory(),FileName];
        FileID          = fopen(OutputFileName,'w','native','UTF-8');
        
        if (FileID ~= -1)
            Print = @(Format,Data) fprintf(FileID,Format,Data);
        else
            fprintf('\n');
            warning('WWPP:VerifyTwoPhaseValues:UnableToOpenFile',...
                ['**** The requested file could not be created, ',...
                'so the verification is being output to the command window. ****']);
            fprintf('\n');
            Print = @(Format,Data) fprintf(Format,Data);
        end
    
    else
        Print = @(Format,Data) fprintf(Format,Data);
        FileID = -1;
    end
    
    
    
    %                              Load, Allocate, and Assign
    % ================================================================================
    
    % CSV Import
    FileName = [CurrentFileDirectory(),'IAPWSCheckValues_TwoPhase.csv'];
    Data     = importdata(FileName, ',', 1);
    
    % Number of data sets
    Nsets  = size(Data.data,1)      ;
    Nprops = size(Data.data,2) - 1  ;
    
    % State variables
    Tsat   = Data.data(:,1);
    
    % IAPWS Values
    IAPWS_P    = Data.data(:,2)     ;
    IAPWS_rhoL = Data.data(:,3)     ;
    IAPWS_rhoG = Data.data(:,4)     ;
    IAPWS_hL   = Data.data(:,5)     ;
    IAPWS_hG   = Data.data(:,6)     ;
    IAPWS_sL   = Data.data(:,7)     ;
    IAPWS_sG   = Data.data(:,8)     ;
    IAPWS      = Data.data(:,2:8)   ;
    
    
    
    %                                Calculations
    % ================================================================================
    
    % Package calculations
    [WWPP_P,WWPP_rhoL,WWPP_rhoG] = SaturationStateGivenTsat(Tsat);
    
    % Round the calculations to the SigFigs in the IAPWS Standard
    Nfigs     = 9;
    WWPP_P    = SigFigRound(WWPP_P,Nfigs);
    WWPP_rhoL = SigFigRound(WWPP_rhoL,Nfigs);
    WWPP_rhoG = SigFigRound(WWPP_rhoG,Nfigs);
    WWPP_hL   = SigFigRound(EnthalpyOne(WWPP_rhoL,Tsat),Nfigs);
    WWPP_hG   = SigFigRound(EnthalpyOne(WWPP_rhoG,Tsat),Nfigs);
    WWPP_sL   = SigFigRound(EntropyOne (WWPP_rhoL,Tsat),Nfigs);
    WWPP_sG   = SigFigRound(EntropyOne (WWPP_rhoG,Tsat),Nfigs);
    WWPP      = [WWPP_P,WWPP_rhoL,WWPP_rhoG,WWPP_hL,WWPP_hG,WWPP_sL,WWPP_sG];
    
    % Absolute Differences
    AbsDiff_P    = abs(IAPWS_P    - WWPP_P   )  ;
    AbsDiff_rhoL = abs(IAPWS_rhoL - WWPP_rhoL)  ;
    AbsDiff_rhoG = abs(IAPWS_rhoG - WWPP_rhoG)  ;
    AbsDiff_hL   = abs(IAPWS_hL   - WWPP_hL  )  ;
    AbsDiff_hG   = abs(IAPWS_hG   - WWPP_hG  )  ;
    AbsDiff_sL   = abs(IAPWS_sL   - WWPP_sL  )  ;
    AbsDiff_sG   = abs(IAPWS_sG   - WWPP_sG  )  ;
    AbsDiff      = [AbsDiff_P,AbsDiff_rhoL,AbsDiff_rhoG,AbsDiff_hL,AbsDiff_hG,AbsDiff_sL,AbsDiff_sG];
    
    
    % Relative Differences (w.r.t to IAPWS)
    RelDiff_P    = AbsDiff_P    ./ IAPWS_P    ;
    RelDiff_rhoL = AbsDiff_rhoL ./ IAPWS_rhoL ;
    RelDiff_rhoG = AbsDiff_rhoG ./ IAPWS_rhoG ;
    RelDiff_hL   = AbsDiff_hL   ./ IAPWS_hL   ;
    RelDiff_hG   = AbsDiff_hG   ./ IAPWS_hG   ;
    RelDiff_sL   = AbsDiff_sL   ./ IAPWS_sL   ;
    RelDiff_sG   = AbsDiff_sG   ./ IAPWS_sG   ;
    RelDiff      = [RelDiff_P,RelDiff_rhoL,RelDiff_rhoG,RelDiff_hL,RelDiff_hG,RelDiff_sL,RelDiff_sG];




    %                                Printing Setup
    % ================================================================================
    
    % Table (horizontal) padding
    PaddingLeft    = 3;
    PaddingInner   = 6;
    PaddingRight   = 3;
    OverHangHeader = 1;
    
    % DataWidth
    DataWidth = 15;     %   Corresponds to scientific notation with 8 decimal places,
    %   2 exponential places, and persistent sign indication
    
    % First column for the data block
    SatPropertyColumn = char('P_sat [Pa]    ' ,...
        'rho_L [kg/m^3]' ,...
        'rho_G [kg/m^3]' ,...
        'h_L   [J/kg]  ' ,...
        'h_G   [J/kg]  ' ,...
        's_L   [J/kg-K]' ,...
        's_G   [J/kg-K]' );
    
    SatPropertyColumnWidth = size(SatPropertyColumn,2);
    
    
    % Temperature header width
    SatTempColumnWidth = Nsets * DataWidth + (Nsets - 1) * PaddingInner;
    
    % BufferSize
    BufferSize = SatPropertyColumnWidth + 2 * OverHangHeader    + ...
                 PaddingLeft + PaddingInner + PaddingRight      + ...
                 SatTempColumnWidth     + 2 * OverHangHeader    ;
    
    
    % Pads
    PadLeft       = RepeatedString(' ',PaddingLeft );
    PadIn         = RepeatedString(' ',PaddingInner);
    PadRight      = RepeatedString(' ',PaddingRight);
    OverHangSpace = RepeatedString(' ',OverHangHeader);
    OverHangRule  = RepeatedString('-',OverHangHeader);
    
    % Temperature numeric portion of header
    SatTempHeaderString = RepeatedString(' ',SatTempColumnWidth);
    Start = 1;
    for k = 1:Nsets
        if (k < Nsets)
            CenteredTemp = [CenteredString([num2str(Tsat(k),'%+3G'),' [K]'],DataWidth),PadIn];
        else
            CenteredTemp = CenteredString([num2str(Tsat(k),'%+3G'),' [K]'],DataWidth);
        end
        End   = Start + length(CenteredTemp) - 1;
        SatTempHeaderString(Start:End) = CenteredTemp;
        Start = End + 1;
    end
    
    
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
    HeaderLine1 = @() Print('%s\n',[PadLeft,OverHangSpace                       ,...
                                    RepeatedString(' ',SatPropertyColumnWidth)  ,...
                                    OverHangSpace,PadIn,OverHangSpace           ,...
                                    CenteredString('Saturation Temperature',...
                                                SatTempColumnWidth)             ,...
                                    OverHangSpace,PadRight]                     );
    
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
