%__________________________________________________________________________
%                 Temperature trends at Antarctic stations                %
% Data access: https://legacy.bas.ac.uk/met/READER/data.html

% Natalia Silva - natalia3.silva@usp.br
% (2020)


% subtightplot.m & MMap toolbox: https://github.com/snatalias/Toolboxes.git
%__________________________________________________________________________

close all; clear all
gap=[.2,.05]; marg_h=[.1,.05]; marg_w=[.8,.1]; % subtightplot parameters

d = dir('EST*.txt');
for i = 1:length(d)
    
	disp(d(i).name);
	est = load(d(i).name);
    est_names = d(i).name(4:end-20); 
    anual = [];	ver = []; out = []; inv = []; pri = [];
	anos = est(:,1);
	
    for j = 1:length(est(:,1));     
        % T annual mean for each station
        anual = [anual nanmean(est(j,2:13))];       
	
        % Sazonal mean
    	ver = [ver nanmean(est(j,[13,2,3]))]; out = [out nanmean(est(j,4:6))];      
        inv = [inv nanmean(est(j,7:9))]; pri = [pri nanmean(est(j,10:12))];
    end
    %%
    % Find and interpolate outliers (2*sigma)
    ma = nanmean(anual); dpa=nanstd(anual);
    mv = nanmean(ver); dpv = nanstd(ver);
    mo = nanmean(out); dpo = nanstd(out);
    mi = nanmean(inv); dpi = nanstd(inv);
    mp = nanmean(pri); dpp=nanstd(pri);
    %
    condia = find(anual > ma+2.5*dpa | anual < ma-2.5*dpa); anual(condia) = NaN;
    condiv = find(ver > mv+2.5*dpv | ver < mv-2.5*dpv); ver(condiv) = NaN;
    condio = find(out > mo+2.5*dpo | out < mo-2.5*dpo); out(condio) = NaN;
    condii = find(inv > mi+2.5*dpi | inv < mi-2.5*dpi); inv(condii) = NaN;
    condip = find(pri > mp+2.5*dpp | pri < mp-2.5*dpp); pri(condip) = NaN;
    
	x = 1:(length(anual)); inta = isnan(anual);
	anual(inta) = interp1(x(~inta),anual(~inta),x(inta));
    %
	ve = 1:(length(ver)); int = isnan(ver);
	ver(int) = interp1(ve(~int),ver(~int),ve(int));       
	%
	ou = 1:(length(out)); into = isnan(out);
	out(into) = interp1(ou(~into),out(~into),ou(into));
    %
	inve = 1:(length(inv)); intv = isnan(inv);
	inv(intv) = interp1(inve(~intv),inv(~intv),inve(intv));
    %
    pr = 1:(length(pri)); intp = isnan(pri);
	pri(intp) = interp1(pr(~intp),pri(~intp),pr(intp)); 
    
    %% Calculate trend, erro and valor_p
    A = fitlm(anos,anual), trendA = A.Coefficients.Estimate(2);
    erroA = A.Coefficients.SE(2); p_valueA = A.Coefficients.pValue(2);
	yoA = A.Coefficients.Estimate(1); lineanual = yoA+trendA*anos;
    trendA = (trendA*10); yoA = yoA*10; erroA = erroA*10;
    %
    V = fitlm(anos,ver); trendV = V.Coefficients.Estimate(2); 
    erroV = V.Coefficients.SE(2); p_valueV = V.Coefficients.pValue(2);
    yoV = V.Coefficients.Estimate(1); linev = yoV+trendV*anos;
    trendV = (trendV*10); yoV = yoV*10; erroV = erroV*10;
	%
    O = fitlm(anos,out); trendO = O.Coefficients.Estimate(2); 
    erroO = O.Coefficients.SE(2); p_valueO = O.Coefficients.pValue(2);
    yoO = O.Coefficients.Estimate(1); lineo = yoO+trendO*anos;
    trendO = (trendO*10); yoO = yoO*10; erroO = erroO*10;
	%
    I = fitlm(anos,inv); trendI = I.Coefficients.Estimate(2); 
    erroI = I.Coefficients.SE(2);p_valueI = I.Coefficients.pValue(2);
    yoI = I.Coefficients.Estimate(1); linei = yoI+trendI*anos;
    trendI = (trendI*10); yoI = yoI*10; erroI = erroI*10;
    %
	P = fitlm(anos,pri); trendP = P.Coefficients.Estimate(2); 
    erroP = P.Coefficients.SE(2); p_valueP = P.Coefficients.pValue(2);
    yoP = P.Coefficients.Estimate(1); linep = yoP+trendP*anos;
    trendP = (trendP*10); yoP = yoP*10; erroP = erroP*10;
    %%
    if p_valueA < 0.01
        SIGA = ('p < 0.01')
    elseif p_valueA < 0.05
        SIGA = ('p < 0.05')
    elseif p_valueA < 0.1
        SIGA = ('p < 0.1')
    else
        SIGA = ('p > 0.1')
    end
    
    if p_valueV < 0.01
        SIGV = ('p < 0.01')
    elseif p_valueV < 0.05
        SIGV = ('p < 0.05')
    elseif p_valueV < 0.1
        SIGV = ('p < 0.1')
    else
        SIGV = ('p > 0.1')
    end
    
    if p_valueO < 0.01
        SIGO = ('p < 0.01')
    elseif p_valueO < 0.05
        SIGO = ('p < 0.05')
    elseif p_valueO < 0.1
        SIGO = ('p < 0.1')
    else
        SIGO = ('p > 0.1')
    end
    
    if p_valueI < 0.01
        SIGI = ('p < 0.01')
    elseif p_valueI < 0.05
        SIGI =('p < 0.05')
    elseif p_valueI < 0.1
        SIGI =('p < 0.1')
    else
        SIGI =('p > 0.1')
    end

    if p_valueP < 0.01
        SIGP = ('p < 0.01')
    elseif p_valueP < 0.05
        SIGP = ('p < 0.05')
    elseif p_valueP < 0.1
        SIGP = ('p < 0.1')
    else
        SIGP = ('p > 0.1')
    end
   
%     %% PLOT
%     
    strA = mat2str(trendA); strA = strA(1:6);
    strErA = mat2str(erroA); strErA = strErA(1:6);
    textA = ['Tend: ',strA,' \pm ',strErA,'(^oC/decada)']
    
	strV = mat2str(trendV); strV = strV(1:6);
    strErV = mat2str(erroV); strErV = strErV(1:6);
    textV = ['Tend: ', strV, '\pm ', strErV ,'(^oC/decada)']
    
    strO = mat2str(trendO); strO = strO(1:6);
    strErO = mat2str(erroO); strErO = strErO(1:6);
    textO = ['Tend: ', strO, ' \pm ', strErO, '(^oC/decada)']
    
    strI = mat2str(trendI); strI = strI(1:6);
    strErI = mat2str(erroI); strErI = strErI(1:6);
    textI = ['Tend: ', strI, ' \pm ', strErI, '(^oC/decada)']
    
    strP = mat2str(trendP); strP = strP(1:6);
    strErP = mat2str(erroP); strErP = strErP(1:6);
    textP = ['Tend: ', strP, ' \pm ', strErP, '(^oC/decada)']
    
    figure('color','w','position',[108 305 1100 677]) 
	subtightplot(3,1,1);plot(anos,anual);
    xlim([min(anos) max(anos)]); ylim([min(anual)-1 max(anual+1)])
	title(est_names);%text(1952,-1.8,'(a)')
	xlabel('Anos'); ylabel('Temperatura (^oC)');legend('Anual')
    hold on; plot(anos,lineanual,'--k')
    text(min(anos),min(anual)+0.2,textA); text(min(anos),min(anual)-0.5,SIGA)
    hold on; line([2000 2000],[min(anual)-1 max(anual)+1],'Color', [0.6 0.6 0.6],'LineStyle',':')

	subplot(3,2,3); plot(anos,ver','r');
    xlim([min(anos) max(anos)]); ylim([min(ver)-1 max(ver)+1])
	xlabel('Anos'); ylabel('Temperatura (^oC)');legend('Ver')
    hold on; plot(anos,linev,'--k');%text(1952,2,'(b)')
    text(min(anos),min(ver)-0.2,textV,'fontsize',8); text(min(anos),min(ver)-0.7,SIGV,'fontsize',8)
    hold on; line([2000 2000],[min(ver)-2 max(ver)+2],'Color', [0.6 0.6 0.6],'LineStyle',':')
    
	subplot(3,2,4); plot(anos,out','m')
    xlim([min(anos) max(anos)]); ylim([min(out)-1 max(out)+1])
	xlabel('Anos'); ylabel('Temperatura (^oC)');legend('Out')
    hold on; plot(anos,lineo,'--k');%text(1952,-0.5,'(c)')
    text(min(anos),min(out)+1,textO,'fontsize',8); text(min(anos),min(out),SIGO,'fontsize', 8)
    hold on; line([2000 2000],[min(out)-2 max(out)+2],'Color', [0.6 0.6 0.6],'LineStyle',':')
    
	subplot(3,2,5); plot(anos,inv','c')
    xlim([min(anos) max(anos)]); ylim([min(inv)-1 max(inv)+1])
	xlabel('Anos'); ylabel('Temperatura (^oC)');legend('Inv')
    hold on; plot(anos,linei,'--k');%text(1952,-4.5,'(d)')
    text(min(anos),min(inv)+1,textI,'fontsize',8); text(min(anos),min(inv)-0.1,SIGI,'fontsize',8)
    hold on; line([2000 2000],[min(inv)-2 max(inv)+2],'Color', [0.6 0.6 0.6],'LineStyle',':')
    
    subplot(3,2,6); plot(anos,pri','g')
    xlim([min(anos) max(anos)]); ylim([min(pri)-1 max(pri)+1])
	xlabel('Anos'); ylabel('Temperatura (^oC)');legend('Pri')
    hold on; plot(anos,linep,'--k');%text(1952,-2,'(e)')
    text(min(anos),min(pri)+0.8,textP,'fontsize',8); text(min(anos),min(pri)-0.2,SIGP,'fontsize',8)
    hold on; line([2000 2000],[min(pri)-2 max(pri)+2],'Color', [0.6 0.6 0.6],'LineStyle',':')
end
