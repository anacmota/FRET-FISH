clear all
close all
clc

%% Input information
figure('units','normalized','outerposition',[0 0 0.5 0.5])

i = 0;
gene = {'Magix', 'Kdm5c'};
LL = 1000;

for g = 1:2
    FRET_G1 = [];
    FRET_G2 = [];
    
for p = 1:2
    
if strcmp(gene{g}, 'Magix')
    filenum = {'iAM600'; 'iAM600_G2'; 'iAM594'; 'iAM594_G2'};
elseif strcmp(gene{g}, 'Kdm5c')
    filenum = {'iAM602'; 'iAM602_G2'; 'iAM682'; 'iAM682_G2'};
end

%% Extract data from file
for con_treat = 1:2
    
    file = p*2-2 + con_treat;
    
    try
        ID = readtable(['data/fret_' filenum{file} '_autoPick.csv']);
    catch
        disp('file does not exist')
        continue
    end
    
r = 0;
D = nan(LL, 2);
A = nan(LL, 2);
FRET = nan(LL, 2);
for Fields = 1:ID.fov(end)
    NucleiAuto = unique(ID.original_nuclei(ID.fov == Fields));

    for n = 1:length(NucleiAuto)
        Indx = (ID.fov == Fields & ID.original_nuclei == NucleiAuto(n));
        
        r = r + 1;
        D(r, :) = unique(ID.Donor(Indx));
        A(r, :) = unique(ID.Acceptor(Indx));
        FRET(r, :) = unique(ID.Fret(Indx));
    end
end

EfficiencyFRET = FRET./(FRET+D)*100;


i = i + 1;
subplot(2, 4, i)
plotDist(EfficiencyFRET)
if rem(i, 2) == 0
    cycle = '(non-G1)';
else
    cycle = '(G1)';
end
title([gene{g} ' ' cycle], 'Color', 'k')
    
end
end
end






function plotDist(EfficiencyFRET)
EfficiencyFRET = EfficiencyFRET(~isnan(EfficiencyFRET));

pts = linspace(0,100,100);
[f,xi] = ksdensity(EfficiencyFRET, pts);

th = df_ThDialog(xi, f);

LowCondensation = sum(EfficiencyFRET < th)/length(EfficiencyFRET) *100;
HighCondensation = sum(EfficiencyFRET > th)/length(EfficiencyFRET) *100;

area(xi(xi<th), f(xi<th), 'FaceColor', [0.3010 0.7450 0.9330])
hold on
area(xi(xi>th), f(xi>th), 'FaceColor', [0.8500 0.3250 0.0980])
text(15, 0.015, [num2str(round(LowCondensation)) '%'],'FontSize', 8, 'Color', 'k')
text(40, 0.015, [num2str(round(HighCondensation)) '%'],'FontSize', 8, 'Color', 'k')


xlim([0 100])
ylim([0 0.065])

ylabel('Probability Density', 'Color', 'k');
xlabel('FRET efficiency (%)', 'Color', 'k');
set(gca,'XColor', 'k')
set(gca,'YColor', 'k')
set(gca,'linew',1)
set(gca,'FontSize', 8)
text(80, 0.035, ['n = ' num2str(sum(sum(~isnan(EfficiencyFRET))))],'HorizontalAlignment', 'center','FontSize', 8, 'Color', 'k')
end







function th = df_ThDialog(D, f)
% function th = df_ThDialog(D)
% Pick an upper threshold for intensities
%


if ~exist('th', 'var')
    th = median(D);
end

gui.isok = 0; % set to 1 if 'Ok' was pressed at exit
gui.f = figure('Name', 'Set upper threshold for FRET efficiency', 'NumberTitle', 'off');
gui.a = axes('Position', [0.1,0.25,.8,.7]);
% gui.h = histogram('Parent', gui.a, D, 50);
gui.h = plot(D, f);
hold on
ax = axis();
gui.thLine = plot([th, th], [ax(3), ax(4)], 'LineWidth', 2);

set(gui.f, 'WindowButtonDownFcn', @interStart);

gui.thValue = uicontrol('Style', 'text', ...
    'String', '', ...
    'Units', 'Normalized', ...
    'Position', [0.1,0,.8,.2], ...
    'Callback', @ok, ...
    'Parent', gui.f, ...
    'HorizontalAlignment','left', ...
    'FontName', get(0,'FixedWidthFontName'));

gui.ok = uicontrol('Style', 'pushbutton', ...
    'String', 'Ok', ...
    'Units', 'Normalized', ...
    'Position', [0.85,0.05,.1,.1], ...
    'Callback', @ok, ...
    'Parent', gui.f);

setTh(th);

uiwait(gui.f);

if gui.isok == 0
    th = [];
end

if isvalid(gui.f)
    close(gui.f);
end

function ok(varargin)
    gui.isok = 1;
    uiresume();
end

    function interStart(varargin)
        gco
        if gco == gui.h | gco == gui.a
            x = get(gui.a, 'CurrentPoint'); x = x(1);        
            setTh(x);          
        end
        if gco == gui.thLine
            set(gui.f, 'WindowButtonMotionFcn', @lineDrag);  
            set(gui.f, 'WindowButtonUpFcn', @stopDrag);
        end
    end

    function stopDrag(varargin)
            set(gui.f, 'WindowButtonMotionFcn', []);  
    end

    function lineDrag(varargin)
           x = get(gui.a, 'CurrentPoint'); x = x(1);
           setTh(x);
    end

    function setTh(x)
        gui.thLine.XData = ones(1,2)*x;
        th = x;
        gui.thValue.String = sprintf('Nuclei: %d\nTh: %.2e\nAbove: %d\nBelow: %d', numel(D), th, sum(D>th), sum(D<th));     
    end


end