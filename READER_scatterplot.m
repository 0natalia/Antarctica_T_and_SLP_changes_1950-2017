%__________________________________________________________________________
%                 Scatter plot READER data results with m_map
% Data access: https://legacy.bas.ac.uk/met/READER/data.html

% Natalia Silva - natalia3.silva@usp.br
% (2020)
%__________________________________________________________________________

close all; clear all

%% 
addpath('home/natalia/READER_Data')
addpath('home/natalia/rotinas_MATLAB/m_map')
clear all; close all

% Each column refers to previous calculated results (DTR, trends, max, min,
% mean, etc.....)
[~,lat,~] = xlsread('plotar_QGis.xlsx',2,'B2:B25'); lat = str2double(lat);
[~,lon,~] = xlsread('plotar_QGis.xlsx',2,'C2:C25'); lon = str2double(lon);
[~,temp,~] = xlsread('plotar_QGis.xlsx',2,'I2:I25'); temp = str2double(temp);

% vetor c/ positivos
find(temp< 0);

% negativos 
tempneg = temp(ans); lonneg = lon(ans);
latneg = lat(ans);

temp(ans) = 0; temp(isnan(temp)) = 0; tempneg = -tempneg;

temp = temp*1000; tempneg = tempneg*1000;

m_proj('stereographic','lat',-90,'long',0,'radius',30);
m_scatter(lon,lat,temp,'r','filled');hold on
m_scatter(lonneg,latneg,tempneg,'b','filled')
m_grid('ytick',[-70 -80 -90],'xtick',12,'radius',60,'tickdir','out', 'xaxisLocation', 'top','yaxisLocation',...
    'middle','box','on', 'fontsize',7,'linewidth',0.8);
m_coast('color','k','linewidth',0.8);
title('READER T Anual');

