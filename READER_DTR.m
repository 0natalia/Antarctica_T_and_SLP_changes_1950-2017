%__________________________________________________________________________
%                  DTR (Tmax - Tmin) at Antarctic stations                %
%          Calculates DTR and its trends at each selected station         %
% Data access: https://legacy.bas.ac.uk/met/READER/data.html

% Natalia Silva - natalia3.silva@usp.br
% (2020)
%__________________________________________________________________________

close all; clear all

%%
d = dir('EST*.txt');

for i = 7
    disp(d(i).name);           % read file title
    est = load(d(i).name);     % load file
    est_names = d(i).name(4:end-20); % select file neme characteres
    tmin = []; tmax = [];
    anos = est(:,1);	  % years
    
    for j = 1:length(anos)
        tmax(end+1) = max(est(j,2:13));	   % year max temp
        
        tmin(end+1) = min(est(j,2:13));    % year min temp
    end
    %% Interpolate and remove outliers
    
    mmax = nanmean(tmax); stdmax = nanstd(tmax);
    mmin = nanmean(tmin); stdmin = nanstd(tmin);
    
    condicaomax = (tmax > mmax+2*stdmax | tmax < mmax-2*stdmax);
    tmax(condicaomax) = NaN;
    condicaomin = (tmin > mmin+2*stdmin | tmin < mmin-2*stdmin);
    tmin(condicaomin) = NaN;
    
    % Interpolate nan data
    % x(~int) nan values, tmax(~int) not-nan values,
    % x(int) nan indices
    %Tmax
    x = 1:(length(tmax)); int = isnan(tmax);
    tmax(int) = interp1(x(~int),tmax(~int),x(int));
    %Tmin
    y = 1:(length(tmin)); inte = isnan(tmin);
    tmin(inte) = interp1(y(~inte),tmin(~inte),y(inte));
   
    %% Calculo de DTR (diurnal temperature range)
     
    dtr = [];
    for i = 1:length(tmax)
        dtr = [dtr tmax(i)-tmin(i)];
    end
    
    %% Normalize - standard units
%     
%     meanmax = nanmean(tmax); dpmax = nanstd(tmax); maxnormal = [];
%     meanmin = nanmean(tmin); dpmin = nanstd(tmin); minnormal = [];
%     
%     for m = 1:length(anos)
%         normax = ((tmax(m)-meanmax)/(dpmax)); maxnormal(end+1) = normax;
%         normin = ((tmin(m)-meanmin)/(dpmin)); minnormal(end+1) = normin;
%     end
  
    %% TREND -- liner regression
    % y = ax + b (intercept + x1)
    
    linearmax = fitlm(anos,tmax); yomax = linearmax.Coefficients.Estimate(1);
    trend_max = linearmax.Coefficients.Estimate(2);
    erromax = linearmax.Coefficients.SE(2);
    linemax = yomax+trend_max*anos; % p/ plotar linha de tendencia
    trend_max = (trend_max*10); yomax = yomax*10; 
    erromax = erromax*10; %^oC/dec
    p_valuemax = linearmax.Coefficients.pValue(2);
    
    linearmin = fitlm(anos,tmin); yomin = linearmin.Coefficients.Estimate(1);
    trend_min = linearmin.Coefficients.Estimate(2);
    erromin = linearmin.Coefficients.SE(2);
    linemin = yomin+trend_min*anos;
    p_valuemin = linearmin.Coefficients.pValue(2);
    trend_min = (trend_min*10); yomin = yomin*20; erromin = erromin*10;
    
    lineardtr = fitlm(anos,dtr); yodtr = lineardtr.Coefficients.Estimate(1);
    trend_dtr = lineardtr.Coefficients.Estimate(2);
    errodtr = lineardtr.Coefficients.SE(2);
    linedtr = yodtr+trend_dtr*anos;
    p_valuedtr = lineardtr.Coefficients.pValue(2);
    trend_dtr = (trend_dtr*10); yodtr = yodtr*10; errodtr = errodtr*10; %^oC/dec
    
    if p_valuemax < 0.01
        SIG = ('p < 0.01');
    elseif p_valuemax < 0.05
        SIG = ('p < 0.05');
    elseif p_valuemax < 0.1
        SIG = ('p < 0.1');
    else
        SIG = ('p > 0.1');
    end
    
    if p_valuemin < 0.01
        SI = ('p < 0.01');
    elseif p_valuemin < 0.05
        SI = ('p < 0.05');
    elseif p_valuemin < 0.1
        SI = ('p < 0.1');
    else
        SI = ('p > 0.1');
    end
    
    if p_valuedtr < 0.01
        SDTR = ('p < 0.01');
    elseif p_valuedtr < 0.05
        SDTR = ('p < 0.05');
    elseif p_valuedtr < 0.1
        SDTR = ('p < 0.1');
    else
        SDTR = ('p > 0.1');
    end
    %% PLOT
    
    strmax = mat2str(trend_max); strmax = strmax(1:5);
    strErmax = mat2str(erromax); strErmax = strErmax(1:5);
    textmax = ['Tend Max: ',strmax,' \pm ',strErmax,'(^oC/decada)'];
    
    strmin = mat2str(trend_min); strmin = strmin(1:5);
    strErmin = mat2str(erromin); strErmin = strErmin(1:5);
    textmin = ['Tend Min: ',strmin,' \pm ',strErmin,'(^oC/decada)'];
    
    strdtr = mat2str(trend_dtr); strdtr = strdtr(1:5);
    strErdtr = mat2str(errodtr); strErdtr = strErdtr(1:5);
    textdtr = ['Tend DTR: ',strdtr,' \pm ',strErdtr,'(^oC/decada)'];
    
    figure('color','w','position',[108 305 814 677])
    subplot(3,1,1);plot(anos,tmax,'r','linewidth',1.3);
    xlim([min(anos) max(anos)]);ylim([(min(tmax)-2) (max(tmax)+2)])
    xlabel('Anos'); ylabel('Tmax (^oC)');legend('Tmax')
    title([est_names]);text(1990,3,'(a)')
    text(min(anos),min(tmax),textmax); text(min(anos),min(tmax)-1,SIG)
    hold on; plot(anos,linemax,'--k')
    %
    subplot(3,1,2),plot(anos,tmin,'b','linewidth',1.3);
    xlim([min(anos) max(anos)]);ylim([(min(tmin)-2) max(tmin)+2])
    xlabel('Anos'); ylabel('Tmin (^oC)');legend('Tmin')
    text(min(anos),min(tmin)+1,textmin); text(min(anos),min(tmin)-1,SI)
    text(1990,-15.8,'(b)')
    hold on; plot(anos,linemin,'--k')
    %
    subplot(3,1,3);plot(anos,dtr,'k');
    xlim([min(anos) max(anos)]); ylim([min(dtr)-1 max(dtr)+1])
    xlabel('Anos'); ylabel('DTR');legend('DTR');text(1990,22,'(c)')
    text(min(anos),min(dtr)+1.8,textdtr); text(min(anos),min(dtr),SDTR)
    hold on; plot(anos,linedtr,'--k')
end

