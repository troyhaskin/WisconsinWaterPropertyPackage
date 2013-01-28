%% One-Phase DataFile Generatation
Data = [300    0.9965560E+3   0.992418352E-1     0.413018112E+1     0.150151914E+4     0.393062643E+0;
        300    0.1005308E+4   0.200022515E+2     0.406798347E+1     0.153492501E+4     0.387405401E+0;
        300    0.1188202E+4   0.700004704E+3     0.346135580E+1     0.244357992E+4     0.132609616E+0;
        500    0.4350000E+0   0.999679423E-1     0.150817541E+1     0.548314253E+3     0.794488271E+1;
        500    0.4532000E+1   0.999938125E+0     0.166991025E+1     0.535739001E+3     0.682502725E+1;
        500    0.8380250E+3   0.100003858E+2     0.322106219E+1     0.127128441E+4     0.256690919E+1;
        500    0.1084564E+4   0.700000405E+3     0.307437693E+1     0.241200877E+4     0.203237509E+1;
        647    0.3580000E+3   0.220384756E+2     0.618315728E+1     0.252145078E+3     0.432092307E+1;
        900    0.2410000E+0   0.100062559E+0     0.175890657E+1     0.724027147E+3     0.916653194E+1;
        900    0.5261500E+2   0.200000690E+2     0.193510526E+1     0.698445674E+3     0.659070225E+1;
        900    0.8707690E+3   0.700000006E+3     0.266422350E+1     0.201933608E+4     0.417223802E+1];

T   = Data(:,1);
rho = Data(:,2);
P   = Data(:,3);
cv  = Data(:,4);
w   = Data(:,5);
s   = Data(:,6);

Format   = '%+16.8E';
FileName = 'IAPWSCheckValues_OnePhase.csv';
FileID   = fopen(FileName,'w+');
fprintf(FileID,['Temperature [K]  , Density [kg/m^3] , Pressure [Pa]    , '     ,...
                'Isochoric Heat Capacity [J/kg-K] , Sound Speed [m/s] , '  ,...
                ' Entropy [J/kg-K]\n']);
            
Space = @(N) repmat(' ',1,N);
            
for k = 1:length(T)
    fprintf(FileID,[Format,          ' , ' ], T  (k)      );
    fprintf(FileID,[Format,          ' , ' ], rho(k)      );
    fprintf(FileID,[Format,          ' , ' ], P  (k) * 1E6); % Convert MPa => Pa
    fprintf(FileID,[Format,Space(16),' , ' ], cv (k) * 1E3); % Convert kJ => J
    fprintf(FileID,[Format,Space(1), ' , ' ], w  (k)      );
    fprintf(FileID,[Format,            '\n'], s  (k) * 1E3); % Convert kJ => J
end
fclose(FileID);


         
%% Two-Phase DataFile Generatation

Data = [0.275000000E+3 0.450000000E+3 0.625000000E+3;
        0.698451167E-3 0.932203564E+0 0.169082693E+2;
        0.999887406E+3 0.890341250E+3 0.567090385E+3;
        0.550664919E-2 0.481200360E+1 0.118290280E+3;
        0.775972202E+1 0.749161585E+3 0.168626976E+4;
        0.250428995E+4 0.277441078E+4 0.255071625E+4;
        0.283094670E-1 0.210865845E+1 0.380194683E+1;
        0.910660121E+1 0.660921221E+1 0.518506121E+1];

Tsat = Data(1,:); %[K]
Psat = Data(2,:); %[MPa]
rhol = Data(3,:); %[kg/m^3]
rhog = Data(4,:); %[kg/m^3]
hl   = Data(5,:); %[kJ/kg]
hg   = Data(6,:); %[kJ/kg]
sl   = Data(7,:); %[kJ/kg-K]
sg   = Data(8,:); %[kJ/kg-K]

Format   = '%+16.8E';
FileName = 'IAPWSCheckValues_TwoPhase.csv';
FileID   = fopen(FileName,'w+');
fprintf(FileID,['Temperature [K]  , Pressure [Pa]    , ',...
                'Liquid Density [kg/m^3] , Gas Density [kg/m^3] , '     ,...
                'Liquid Enthalpy [J/kg] , Gas Enthalpy [J/kg] , '  ,...
                'Liquid Entropy [J/kg-K] , Gas Entropy [J/kg-K] \n']);
            
for k = 1:length(Tsat)
    fprintf(FileID,[Format,         ' , ' ], Tsat(k)      );
    fprintf(FileID,[Format,         ' , ' ], Psat(k) * 1E6); % Convert MPa => Pa
    fprintf(FileID,[Format,Space(7),' , ' ], rhol(k)      );
    fprintf(FileID,[Format,Space(4),' , ' ], rhog(k)      );
    fprintf(FileID,[Format,Space(6),' , ' ], hl  (k) * 1E3); % Convert kJ => J
    fprintf(FileID,[Format,Space(3),' , ' ], hg  (k) * 1E3); % Convert kJ => J
    fprintf(FileID,[Format,Space(7),' , ' ], sl  (k) * 1E3); % Convert kJ => J
    fprintf(FileID,[Format,         '\n'  ], sg  (k) * 1E3); % Convert kJ => J
end
fclose(FileID);
