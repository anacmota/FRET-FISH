clear all
close all
clc

replicate = 1;

color={'Cy5'; 'tmr'; 'tmrfret'};
design = {'Oligo1 Space50', 'Oligo1 Space150', 'Oligo2 Space50', 'Oligo2 Space300', 'Oligo4 Space50', 'Oligo4 Space300'};

figure('units','normalized','outerposition',[0 0 0.25 0.7])
for i = 1:6
    % There are 6 different designs tested
    if replicate == 1
        ID = readtable(['datasets/MEF/wCentr.out.dilate0.iAM295_00', num2str(i),'.csv']);
    elseif replicate == 2
        ID = readtable(['datasets/NIH3T3/iAM324_00', num2str(i),'.csv']);
    end
    
Donor = [];
FRET = [];

% Reading the number of the last field from the long extension
if iscell(ID.File(end))
    Temp = [];
    for c = 1:numel(ID.File)
        Temp = [Temp; str2double(ID.File{c}(end-5:end-3))];
    end
    % Simplifies the calling
    ID.File = Temp;
end
LastIndex = ID.File(end);

for Fields = 1:LastIndex
    % Extract information by each field of view
    
    PositionFile =  find(ID.File==Fields); % Select current field of view
    nuclei = unique(ID.Nuclei(PositionFile)); % Identify the Nuclei to be analysed
    
    for n = 1:length(nuclei)
        % Extract information by nuclei
        PositionNuclei = PositionFile(eq(nuclei(n), ID.Nuclei(PositionFile))); % Select current nuclei
        
        [F, D] = get_Intensity(ID.Channel(PositionNuclei), ID.Value(PositionNuclei), ...
            [ID.x(PositionNuclei)*0.13, ID.y(PositionNuclei)*0.13, ID.z(PositionNuclei)*0.3], color);

        FRET = [FRET, F];
        Donor = [Donor, D];
    end
end

%% Calculate FRET efficiency
FRET_eff = FRET./(Donor + FRET)*100;

subplot(3, 2, i)
yyaxis left
histogram(FRET_eff, 20)
ylabel('# FRET pairs', 'Color', 'k')
set(gca,'YColor', 'k')
text(60, 10, ['n = ' num2str(length(FRET_eff))],'FontSize', 8, 'Color', 'k')

hold on

yyaxis right
plotDist(FRET_eff, '-k');
xlim([0 85])
set(gca,{'linew'},{1})
xlabel('FRET efficiency (%)', 'Color', 'k')
ylabel('Probability Density', 'Color', 'k')
set(gca,'XColor', 'k')
set(gca,'YColor', 'k')
title(design{i}, 'Color', 'k')
set(gca,'FontSize', 8)
end




function plotDist(EfficiencyFRET, color)
% probability density of FRET efficiency values on top of the histogram
[f,xi] = ksdensity(EfficiencyFRET);
plot(xi,f, color, 'LineWidth', 1);
end




function [F, D] = get_Intensity(Channel, Values, position, color)
% Remove values higher than 10^4
% Associate the closest dots from different channels with less than 0.5µm
% Extract the Donor and FRET intensities
if any(Values > 10^4)
    Channel(Values > 10^4) =  [];
    position(Values > 10^4, :) =  [];
    Values(Values > 10^4) =  [];
    
    if numel(unique(Channel)) ~= 3
        F = [];
        D = [];
        return
    end
end

F = nan(1, 9);
D = nan(1, 9);

thr = 0.5;

A_pos = position(strcmp(Channel, color{1}), :);
D_idx = find(strcmp(Channel, color{2}));
F_idx = find(strcmp(Channel, color{3}));
for dd = 1:size(A_pos, 1)
    % Find the brighest acceptor dots and eliminate dots from the same
    % channel within the threshold
    A_dist = pdist2(A_pos(dd, :), A_pos);
    A_dist(A_dist == 0) = nan;
    if any(A_dist < thr)
        V = max([Values(dd); Values(A_dist < thr)]);
        if find(V == Values) ~= dd
            continue
        end
    end
    
    % Find the closest donor dot to the acceptor
    [m, idx] = min(pdist2(A_pos(dd, :), position(D_idx, :)));
    if m < thr
    D(dd) = Values(D_idx(idx));
    end
    
    % Find the closest fret dot to the acceptor
    [m, idx] = min(pdist2(A_pos(dd, :), position(F_idx, :)));
    if m < thr
    F(dd) = Values(F_idx(idx));
    end
    
    % If the donor or acceptor dot were not found within the threshold, is
    % attributed non-specific value to donor and acceptor
    if isnan(D(dd)) || isnan(F(dd))
        D(dd) = nan;
        F(dd) = nan;
    end
end

% all the cases that donor and fret were not clearly identified, then the
% values are turn empty
D(isnan(D)) = [];
F(isnan(F)) = [];
end