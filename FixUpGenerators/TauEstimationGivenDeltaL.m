% clc;
clear('all');



Tc  = CriticalTemperature();
Tt  = TriplePointTemperature();

Nsamp  = 30E3 * [1;1;1;10];
Norder = [10;10;10;10];

ZoneEdges = [1,1.001,1.4,2.3,Tc./Tt];

Abscissa = {@(x) x , @(x) x , @(x) x , @(x) x };
Ordinate  = {@(x) x , @(x) x , @(x) x , @(x) x };

figure(1);
clf;
set(gca(),'NextPlot','add');
H = gca();

for k = 4%1 : (length(ZoneEdges)-1)
    tau = ChebyshevSpace(ZoneEdges(k),ZoneEdges(k+1),Nsamp(k))';
    [~,delL,~] = SaturationStateGivenTauRRND(tau);
    
%     tau-delG fit
%     [pG,SG,muG] = polyfit(AbscissaG{k}(delG),OrdinateG{k}(tau),NorderG(k));
%     delGfit     = logspace(log10(min(delG)),log10(max(delG)),30)';
%     tauGFit     = polyval(pG,AbscissaG{k}(delGfit),SG,muG);
    
    % tau-delL fit
    [p,S,mu] = polyfit(Abscissa{k}(delL),Ordinate{k}(tau),Norder(k));
    delLfit     = linspace(min(delL),max(delL),30)';
    tauFit      = polyval(p,Abscissa{k}(delLfit),S,mu);
    
     % tau-delL plot
     semilogx(delL,tau,delLfit,tauFit,'o');
    
    % tau-delL fit error norm
    Linf = norm(interp1(delL,tau,delLfit,'pchip') - tauFit,Inf);
    
    % Data table
    fprintf('\n');
    fprintf([repmat(' ',1,31),'tau-del Fit Data',repmat(' ',1,32),'\n']);
    fprintf([repmat('=',1,80),'\n']);
    fprintf('Edge Values from fit:\n');
    fprintf('     Min: {del,tau} = {%+9.5E,%+9.5E}\n',min(delLfit),min(tauFit));
    fprintf('     Max: {del,tau} = {%+9.5E,%+9.5E}\n',max(delLfit),max(tauFit));
    %
    fprintf('Fit Parameter Functions:\n');
    %
    StringAbscissa = func2str(Abscissa{k});
    fprintf('     Abscissa: %s\n',StringAbscissa);
    %
    StringOrdinate = func2str(Ordinate{k});
    fprintf('     Ordinate: %s\n',StringOrdinate);
    %
    fprintf('L_inf Norm:\n');
    fprintf('     %+9.5E\n',Linf);
    %
    fprintf('Fit Constants:\n');
    %
    fprintf('     mu:\n');
    fprintf('          %+23.16E\n',mu);
    fprintf('     p:\n');
    fprintf('          %+23.16E\n',p);
    fprintf([repmat('=',1,80),'\n']);
    fprintf('\n');
    fprintf('\n');
end


