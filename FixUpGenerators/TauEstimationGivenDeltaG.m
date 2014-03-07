% clc;
clear('all');



Tc  = CriticalTemperature();
Tt  = TriplePointTemperature();

Nsamp  = 30E3 * [1;2;3];
Norder = [12;10;8];

ZoneEdges = [1,1.001,1.4,Tc./Tt];

Abscissa = {@(x) log10(x) , @(x) log10(x) , @(x) log10(x)};
Ordinate = {@(x)    x     , @(x)    x     , @(x)    x    };

figure(1);
clf;
set(gca(),'NextPlot','add');
H = gca();

for k = 1 : (length(ZoneEdges)-1)
    tau = ChebyshevSpace(ZoneEdges(k),ZoneEdges(k+1),Nsamp(k))';
    [~,delL,delG] = SaturationStateGivenTauRRND(tau);
    
    % tau-delG fit
    [p,S,mu] = polyfit(Abscissa{k}(delG),Ordinate{k}(tau),Norder(k));
    delGfit     = logspace(log10(min(delG)),log10(max(delG)),30)';
    tauFit     = polyval(p,Abscissa{k}(delGfit),S,mu);
    
     % tau-delG plot
     semilogx(delG,tau,delGfit,tauFit,'o');
    
    % tau-delG fit error norm
    Linf = norm(interp1(delG,tau,delGfit,'pchip') - tauFit,Inf);
    
    % Data table
    % Data table
    fprintf('\n');
    fprintf([repmat(' ',1,31),'tau-del Fit Data',repmat(' ',1,32),'\n']);
    fprintf([repmat('=',1,80),'\n']);
    fprintf('Edge Values from fit:\n');
    fprintf('     Min: {del,tau} = {%+9.5E,%+9.5E}\n',min(delGfit),min(tauFit));
    fprintf('     Max: {del,tau} = {%+9.5E,%+9.5E}\n',max(delGfit),max(tauFit));
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


