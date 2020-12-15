%__________________________________________________________________________
%                    Regional trends (T, SLP and WS)                
%                       Annual and seasonal trends        
% Data access: https://legacy.bas.ac.uk/met/READER/data.html

% Antarctic regions:
% A - Bellingshausen (2), Esperanza(6), Faraday (7), Great Wall (9),      
%     Marambio (12), Marsh (14), Orcadas (20), Rothera(21), San Martin(22)
% B - Halley (10), Neumayer (18), Novolazarevskaya (19), Syowa (23)       
% C - Casey (3), Davis (4), Dummond Durvile (5), Mawson (15),             
%     Mirny (17), Zhongshan (25)                                          
% D - McMurdo (16)                                                        
% E - Amundsen Scott (1), Vostok (24)           

% Natalia Silva - natalia3.silva@usp.br
% (2020)
%__________________________________________________________________________
%%
close all; clear all
gap=[.2,.05]; marg_h=[.1,.05]; marg_w=[.8,.1];
d = dir('EST*.txt');

A = [2, 6, 7, 9, 12, 14, 20, 21, 22]; B = [10, 18, 19, 23];
C = [3, 4, 5, 15, 17, 25] ;D = 16; E = [1, 24];
e = {A B C D E}; y = {(1903:2017) (1957:2017) (1954:2017) (1956:2017) (1957:2017)};
reg_names = {'Peninsula','Oc Indico/Mar de Wedell','Oc Pacifico','Mar de Ross',...
    'Continental'};

for w = 3
    yrange = length(y{w});
    tempmax = repmat(-99999,1,yrange); tempmin = repmat(99999,1,yrange);
    anual = []; ver = []; out = []; inv = []; pri = [];
    p = 1; k = 1; % position count
 
    for i = e{w}
        disp(d(i).name);
        est = load(d(i).name);
        est_names = d(i).name(4:end-20);
        anos = est(:,1);
        % initial yr
        inicio = yrange - length(anos);   
        % each loop adds space for new data/station
        anual = [anual NaN(yrange,12)]; ver = [ver NaN(yrange,3)]; 
        out = [out NaN(yrange,3)]; inv = [inv NaN(yrange,3)]; pri = [pri NaN(yrange,3)];
      
        for j = 1:length(anos)
            temax = max(est(j,2:13)); temin = min(est(j,2:13)); 
            if temax > tempmax(j+inicio)
                tempmax(j+inicio) = temax;
            end
            if temin < tempmin(j+inicio)
                tempmin(j+inicio) = temin;
            end
            
            anual(j+inicio,p:p+11) = est(j,2:end);
            % adjust data vec
            ver(j+inicio,k) = est(j,13); ver(j+inicio,k+1:k+2) = est(j,2:3);
            out(j+inicio,k:k+2) = est(j,4:6); inv(j+inicio,k:k+2) = est(j,7:9); 
            pri(j+inicio,k:k+2) = est(j,10:12);
        end
        p = p+12; k = k+3;  
    end
    
    % mean temperature vec
    manu = []; mver = []; mout = []; minv = []; mpri = [];
    for j = 1:yrange
        ma = nanmean(anual(j,:)); mv = nanmean(ver(j,:)); mo = nanmean(out(j,:));
        mi = nanmean(inv(j,:)); mp = nanmean(pri(j,:));
        manu = [manu ma]; mver = [mver mv]; mout = [mout mo]; 
        minv = [minv mi]; mpri = [mpri mp];
    end
 
    % Remove outliers
    cmax = find((tempmax > mean(tempmax)+2.5*std(tempmax))|tempmax < mean(tempmax)-2.5*std(tempmax));
    cmin = find((tempmin > mean(tempmin)+2.5*std(tempmin))|tempmin < mean(tempmin)-2.5*std(tempmin));
    tempmax(cmax) = NaN; tempmin(cmin) = NaN;
     
    condia = find(manu > ma+2.5*std(manu) | manu < ma-2.5*std(manu)); manu(condia) = NaN;
    condiv = find(mver > mv+2.5*std(mver) | mver < mv-2.5*std(mver)); mver(condiv) = NaN;
    condio = find(mout > mo+2.5*std(mout) | mout < mo-2.5*std(mout)); mout(condio) = NaN;
    condii = find(minv > mi+2.5*std(minv) | minv < mi-2.5*std(minv)); minv(condii) = NaN;
    condip = find(mpri > mp+3.5*std(mpri) | mpri < mp-3.5*std(mpri)); mpri(condip) = NaN;
%     
    % Interpolate
    mx = 1:(length(tempmax)); intmx = isnan(tempmax);
	tempmax(intmx) = interp1(mx(~intmx),tempmax(~intmx),mx(intmx));
    
    mn = 1:(length(tempmin)); intmn = isnan(tempmin);
	tempmin(intmn) = interp1(mn(~intmn),tempmin(~intmn),mn(intmn));
   
	x = 1:(length(manu)); inta = isnan(manu);
	manu(inta) = interp1(x(~inta),manu(~inta),x(inta));
    
	ve = 1:(length(mver)); int = isnan(mver);
	mver(int) = interp1(ve(~int),mver(~int),ve(int));       
	
	ou = 1:(length(mout)); into = isnan(mout);
	mout(into) = interp1(ou(~into),mout(~into),ou(into));
    
	inve = 1:(length(minv)); intv = isnan(minv);
	minv(intv) = interp1(inve(~intv),minv(~intv),inve(intv));
    
    pr = 1:(length(mpri)); intp = isnan(mpri);
	mpri(intp) = interp1(pr(~intp),mpri(~intp),pr(intp)); 
% %%
%     % Normalize standard units
%     meanmax = nanmean(tempmax); dpmax = nanstd(tempmax);
%     meanmin = nanmean(tempmin); dpmin = nanstd(tempmin);
%     meanmean = mean(manu); dpmean = std(manu);
%     maxnormal = []; minnormal = []; meannormal = [];
%     for z = 1:length(tempmax)
%         normax = ((tempmax(z)-meanmax)/(dpmax));
%         maxnormal = [maxnormal normax];
%             
%         normin = ((tempmin(z)-meanmin)/(dpmin));
%         minnormal = [minnormal normin];
%     
%         normean = ((manu(z)-meanmean)/dpmean);
%         meannormal = [meannormal normean];
%     end

%%
% DTR
dtr = [];
for i = 1:length(tempmax)
    dtr = [dtr tempmax(i)-tempmin(i)];
end
%% Trens
    % Max e Min
    linearmax = fitlm(y{w},tempmax);  p_valuemax = linearmax.Coefficients.pValue(2);
    trend_max = linearmax.Coefficients.Estimate(2); erromax = linearmax.Coefficients.SE(2);
    yomax = linearmax.Coefficients.Estimate(1); linemax = yomax+trend_max*y{w}; 
    trend_max = (trend_max*10); erromax = erromax*10; 

    linearmin = fitlm(y{w},tempmin); p_valuemin = linearmin.Coefficients.pValue(2);
    trend_min = linearmin.Coefficients.Estimate(2); erromin = linearmin.Coefficients.SE(2);
    yomin = linearmin.Coefficients.Estimate(1); linemin = yomin+trend_min*y{w};
    trend_min = (trend_min*10); erromin = erromin*10;
    
    lineardtr = fitlm(y{w},dtr); p_valuedtr = lineardtr.Coefficients.pValue(2);
    trend_dtr = lineardtr.Coefficients.Estimate(2); errodtr = lineardtr.Coefficients.SE(2);
    yodtr = lineardtr.Coefficients.Estimate(1); linedtr = yodtr+trend_dtr*y{w};
    trend_dtr = (trend_dtr*10); errodtr = errodtr*10;
    
    % Means
    linan = fitlm(y{w},manu);  tan = linan.Coefficients.Estimate(2); 
    erroan = linan.Coefficients.SE(2);  p_a = linan.Coefficients.pValue(2); 
    yoa = linan.Coefficients.Estimate(1); linean = (yoa+(tan*y{w}));
    tan = (tan*10); erroan = erroan*10;
    %
    linv = fitlm(y{w},mver);  tver = linv.Coefficients.Estimate(2); 
    errover = linv.Coefficients.SE(2); p_v = linv.Coefficients.pValue(2);
    yov = linv.Coefficients.Estimate(1); linev = (yov+(tver*y{w}));
    tver = (tver*10); errover = errover*10;
    %
    lino = fitlm(y{w},mout); tout = lino.Coefficients.Estimate(2); 
    erroout = lino.Coefficients.SE(2); p_o = lino.Coefficients.pValue(2);
    yoo = lino.Coefficients.Estimate(1); lineo = (yoo+(tout*y{w})); 
    tout = (tout*10); erroout = erroout*10;
    %
    lini = fitlm(y{w},minv); tinv = lini.Coefficients.Estimate(2); 
    erroinv = lini.Coefficients.SE(2); p_i = lini.Coefficients.pValue(2); 
    yoi = lini.Coefficients.Estimate(1); linei = (yoi+(tinv*y{w}));
    tinv = (tinv*10); erroinv = erroinv*10;
    %
    linp = fitlm(y{w},mpri); tpri = linp.Coefficients.Estimate(2); 
    erropri = linp.Coefficients.SE(2); p_p = linp.Coefficients.pValue(2);
    yop = linp.Coefficients.Estimate(1); linep = (yop+(tpri*y{w}));
    tpri = (tpri*10); erropri = erropri*10;

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

    if p_a < 0.01
        SIGA = ('p < 0.01');
    elseif p_a < 0.05
        SIGA = ('p < 0.05');
    elseif p_a < 0.1
        SIGA = ('p < 0.1');
    else
        SIGA = ('p > 0.1');
    end
    
    if p_v < 0.01
        SIGV = ('p < 0.01');
    elseif p_v < 0.05
        SIGV = ('p < 0.05');
    elseif p_v < 0.1
        SIGV = ('p < 0.1');
    else
        SIGV = ('p > 0.1');
    end
    
    if p_o < 0.01
        SIGO = ('p < 0.01');
    elseif p_o < 0.05
        SIGO = ('p < 0.05');
    elseif p_o < 0.1
        SIGO = ('p < 0.1');
    else
        SIGO = ('p > 0.1');
    end
    
    if p_i < 0.01
        SIGI = ('p < 0.01');
    elseif p_i < 0.05
        SIGI =('p < 0.05');
    elseif p_i < 0.1
        SIGI =('p < 0.1');
    else
        SIGI =('p > 0.1');
    end

    if p_p < 0.01
        SIGP = ('p < 0.01');
    elseif p_p < 0.05
        SIGP = ('p < 0.05');
    elseif p_p < 0.1
        SIGP = ('p < 0.1');
    else
        SIGP = ('p > 0.1');
     end
% % PLOT 1
%     strmax = mat2str(trend_max); strmax = strmax(1:5); strErmax = mat2str(erromax); 
%     strErmax=strErmax(1:5); textmax=['Tend Max: ',strmax,' \pm ',strErmax,' (^oC/decada)'];
%     
%     strmin = mat2str(trend_min); strmin = strmin(1:5); strErmin = mat2str(erromin); 
%     strErmin=strErmin(1:5);textmin=['Tend Min: ',strmin,' \pm ',strErmin,' (^oC/decada)'];
%     
%     strdtr = mat2str(trend_dtr); strdtr = strdtr(1:5); strErdtr = mat2str(errodtr); 
%     strErdtr=strErdtr(1:5);textdtr=['DTR: ',strdtr,' \pm ',strErdtr,' (^oC/decada)'];
%     
%     % Tmax
%     figure('color','w','position',[108 305 814 677]);
%     subplot(3,1,1);plot(y{w},tempmax,'r','linewidth',1.3);
%     xlim([min(y{w}) 2017]);ylim([min(tempmax)-1 max(tempmax)+1]);
%     xlabel('Anos'); ylabel('Tmax (^oC)');legend('Tmax')
%     title('Peninsula');text(1904,3.5,'(a)')
%     text(1910,mean(tempmax)+2,textmax); text(1910,mean(tempmax)+1.2,SIG)
%     hold on; plot(y{w},linemax,'--k')
%     
%     % Tmin
%     subplot(3,1,2);plot(y{w},tempmin,'b','linewidth',1.3);
%     xlim([min(y{w}) 2017]); ylim([min(tempmin)-1 max(tempmin)+1])
%     xlabel('Anos'); ylabel('Tmin (^oC)');legend('Tmin')
%     text(1910,min(tempmin)+2,textmin); text(1910,min(tempmin),SI)
%     hold on; plot(y{w},linemin,'--k');text(1904,-9,'(b)')
%     
%     % DTR
%     subplot(3,1,3); plot(y{w},dtr,'k');xlim([min(y{w}) 2017]);
%     xlabel('Anos'); ylabel('DTR');legend('DTR')
%     text(1910,min(dtr)-1.5,textdtr); text(1910,min(dtr)-6,SDTR)
%     hold on; plot(y{w},linedtr,'--k');text(1904,28,'(c)')
%     
%% PLOT 2
    stra = mat2str(tan); stra = stra(1:5); 
    strEra = mat2str(erroan); strEra = strEra(1:5);
    texta = ['Tend: ',stra,' \pm ',strEra,' (^oC/decada)'];

    strV = mat2str(tver); strV = strV(1:6);
    strErV = mat2str(errover); strErV = strErV(1:6);
    textV = ['Tend: ', strV, '\pm ', strErV ,'(^oC/decada)'];
    
    strO = mat2str(tout); strO = strO(1:6);
    strErO = mat2str(erroout); strErO = strErO(1:6);
    textO = ['Tend: ', strO, ' \pm ', strErO, '(^oC/decada)'];
    
    strI = mat2str(tinv); strI = strI(1:6);
    strErI = mat2str(erroinv); strErI = strErI(1:6);
    textI = ['Tend: ', strI, ' \pm ', strErI, '(^oC/decada)'];
    
    strP = mat2str(tpri); strP = strP(1:6);
    strErP = mat2str(erropri); strErP = strErP(1:6);
    textP = ['Tend: ', strP, ' \pm ', strErP, '(^oC/decada)'];
    
    % Anual
    figure('color','w','position',[108 305 1100 677]) 
	subtightplot(3,1,1);plot(y{w},manu);
    xlim([min(y{w}) max(y{w})]); ylim([min(manu)-1 max(manu)+1]);legend('Anual')
	title([reg_names{w}])
	xlabel('Anos'); ylabel('T (^oC)')
    hold on; plot(y{w},linean,'--k');text(1958,-9,'(a)')
    text(min(y{w}),min(manu),texta); text(min(y{w}),min(manu)-0.5,SIGA)
    hold on; line([2000 2000],[min(manu)-1 max(manu)+1],'Color', [0.6 0.6 0.6],'LineStyle',':')

    % Sazonal
	subplot(3,2,3); plot(y{w},mver','r');
    xlim([min(y{w}) max(y{w})]); ylim([min(mver)-1 max(mver)+1])
	xlabel('Anos'); ylabel('T (^oC)');legend('Ver')
    hold on; plot(y{w},linev,'--k');text(1958,0,'(b)')
    text(min(y{w}),min(mver),textV,'fontsize',8); text(min(y{w}),min(mver)-0.5,SIGV,'fontsize',8)
    hold on; line([2000 2000],[min(mver)-2 max(mver)+2],'Color', [0.6 0.6 0.6],'LineStyle',':')
    
	subplot(3,2,4); plot(y{w},mout','m')
    xlim([min(y{w}) max(y{w})]); ylim([min(mout)-1 max(mout)+1])
	xlabel('Anos'); ylabel('T (^oC)');legend('Out')
    hold on; plot(y{w},lineo,'--k');text(1958,-10,'(c)')
    text(min(y{w}),min(mout),textO,'fontsize',8); text(min(y{w}),min(mout)-0.5,SIGO,'fontsize', 8)
    hold on; line([2000 2000],[min(mout)-2 max(mout)+2],'Color', [0.6 0.6 0.6],'LineStyle',':')
    
	subplot(3,2,5); plot(y{w},minv','c')
    xlim([min(y{w}) max(y{w})]); ylim([min(minv)-1 max(minv)+1])
	xlabel('Anos'); ylabel('T (^oC)');legend('Inv')
    hold on; plot(y{w},linei,'--k');text(1958,-15,'(d)')
    text(min(y{w}),min(minv)+1,textI,'fontsize',8); text(min(y{w}),min(minv),SIGI,'fontsize',8)
    hold on; line([2000 2000],[min(minv)-2 max(minv)+2],'Color', [0.6 0.6 0.6],'LineStyle',':')
    
    subplot(3,2,6); plot(y{w},mpri','g')
    xlim([min(y{w}) max(y{w})]); ylim([min(mpri)-1 max(mpri)+1])
	xlabel('Anos'); ylabel('T (^oC)');legend('Pri')
    hold on; plot(y{w},linep,'--k');text(1958,-9,'(e)')
    text(min(y{w}),min(mpri)+1,textP,'fontsize',8); text(min(y{w}),min(mpri),SIGP,'fontsize',8)
    hold on; line([2000 2000],[min(mpri)-2 max(mpri)+2],'Color', [0.6 0.6 0.6],'LineStyle',':')
 
 %% PLOT 3
%     figure('color','w','position',[100 300 850 380])
%     plot(y{w},maxnormal,'m');xlim([min(y{w}) 2017]);
%     title(['Variabilidade de Temperatura: ', reg_names{w}])
%     xlabel('Anos'); ylabel('T (^oC)')
%     hold on; plot(y{w},minnormal)
%     hold on; plot(y{w},meannormal,'k','linewidth',1.5)
%     legend('Tmax','Tmin','Tmedia')
end
