% clc;
clear('all');



Tc  = CriticalTemperature();
Tt  = TriplePointTemperature();

Nsamp  = 25E3 * [1;1;1];
Norder = 13   * [1;1;1];

ZoneEdges = [1,1.001,1.4,Tc./Tt];

% AbscissaG = {@(x) log10(x) , @(x) log10(x) , @(x) log10(x)};
% OrdinateG = {@(x)    x     , @(x)    x     , @(x)    x    };

AbscissaG = {@(x) x , @(x) x , @(x) x };
OrdinateG = {@(x) x , @(x) x , @(x) x };

for k = 1 : 1%(length(ZoneEdges)-1)
    tau = ChebyshevSpace(ZoneEdges(k),ZoneEdges(k+1),Nsamp(k))';
    [~,delG,delL] = SaturationStateGivenTauRRND(tau);
    
    % tau-delG fit
    [pG,SG,muG] = polyfit(AbscissaG{k}(delG),OrdinateG{k}(tau),Norder(k));
    delGfit     = logspace(AbscissaG{k}(min(delG)),AbscissaG{k}(max(delG)),30)';
    tauFit      = polyval(pG,AbscissaG{k}(delGfit),SG,muG);
    
%     % tau-delG fit
%     [pG,SG,muG] = polyfit(AbscissaG{k}(delG),OrdinateG{k}(tau),Norder(k));
%     delGfit     = logspace(AbscissaG{k}(min(delG)),AbscissaG{k}(max(delG)),30)';
%     tauFit      = polyval(pG,AbscissaG{k}(delGfit),SG,muG);
    
    % tau-delG plot
    figure(k);
    semilogx(delG,tau,delGfit,tauFit,'o');
    
    % tau-delG fit error norm
    Linf = norm(interp1(delG,tau,delGfit,'pchip') - tauFit,Inf);
    
    % Data table
    fprintf('\n');
    fprintf([repmat(' ',1,31),'tau-del Fit Data',repmat(' ',1,32),'\n']);
    fprintf([repmat('=',1,80),'\n']);
    fprintf(repmat(' ',1,25));
    fprintf('tau Values:\n');
    fprintf(repmat(' ',1,25));
    fprintf('     Min: %+9.5E\n',ZoneEdges(k));
    fprintf(repmat(' ',1,25));
    fprintf('     Max: %+9.5E\n',ZoneEdges(k+1));
    %
    fprintf('\n');
    fprintf(['              delL',repmat(' ',1,29),'              delG\n']);
    fprintf([repmat('-',1,80),'\n']);
    %
    fprintf(['Fit Parameter Functions:',repmat(' ',1,24),'Fit Parameter Functions:\n']);
    %
    StringAbscissaG = func2str(AbscissaG{k});
    fprintf('     Abscissa: %s',StringAbscissaG);
    fprintf(repmat(' ',1,48-15-length(StringAbscissaG)));
    fprintf('     Abscissa: %s\n',StringAbscissaG);
    %
    StringOrdinateG = func2str(OrdinateG{k});
    fprintf('     Ordinate: %s',StringOrdinateG);
    fprintf(repmat(' ',1,48-15-length(StringOrdinateG)));
    fprintf('     Ordinate: %s\n',StringOrdinateG);
    %
    fprintf(['L_inf Norm:',repmat(' ',1,37),'L_inf Norm:\n']);
    fprintf('     %+9.5E',Linf);
    fprintf(repmat(' ',1,46-5-9));
    fprintf('     %+9.5E\n',Linf);
    %
    fprintf(['Fit Constants:',repmat(' ',1,34),'Fit Constants:\n']);
    %
    fprintf(['     mu:'          ,repmat(' ',1,39),'     mu:\n'                    ]);
    fprintf(['          %+23.16E',repmat(' ',1,14),'          %+23.16E\n'],[muG,muG]');
    fprintf(['     p:'           ,repmat(' ',1,40),'     p:\n'                    ]);
    fprintf(['          %+23.16E',repmat(' ',1,14),'          %+23.16E\n'],[pG,pG]');
    fprintf([repmat('=',1,80),'\n']);
    fprintf('\n');
    fprintf('\n');
end


