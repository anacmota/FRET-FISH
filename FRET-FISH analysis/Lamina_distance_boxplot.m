clear all
close all
clc

gene = {'Atp2b3', 'Magix', 'Kdm5c'};
for g = 1:3
    AlleleG = [];
    geneID = [];
    %% Input information
    filenum = {};
        
    if strcmp(gene{g}, 'Kdm5c')
        filenum = {'iAM592'; 'iAM682'; 'iAM602'};
    elseif strcmp(gene{g}, 'Magix')
        filenum = {'iAM594'; 'iAM680'; 'iAM600'};
    elseif strcmp(gene{g}, 'Atp2b3')
        filenum = {'iAM618'; 'iAM684'; ''};
    end
    
    %% Extract data from file
    for file = 1:3
        
        if isempty(filenum) || isempty(filenum{file})
            disp('file does not exist')
            continue
        end
        
        try
            ID = readtable(['data/dlamin_fret_' filenum{file} '_autoPick.csv']);
        catch
            disp('file does not exist')
            continue
        end
        
        r = 0;
        D = [];
        FRET = [];
        Lam = [];
        for Fields = 1:ID.fov(end)
            
            NucleiAuto = unique(ID.original_nuclei(ID.fov == Fields));
            
            for n = 1:length(NucleiAuto)
                Indx = (ID.fov == Fields & ID.original_nuclei == NucleiAuto(n));
                
                r = r + 1;
                D(r, :) = unique(ID.Donor(Indx));
                FRET(r, :) = unique(ID.Fret(Indx));
                Lam(r, :) = unique(ID.dlamin_A(Indx));
            end
        end
        
        EfficiencyFRET = FRET./(FRET+D)*100;
        EfficiencyFRET = EfficiencyFRET(:);
        Lam = Lam(:);
        
        Div = quantile(Lam, [0 0.25 0.50 0.75 1]);
        G = [];
        grp = [];
        for q = 1:4
            % group the FRET efficiencies by the 4 concentric regions in
            % the nucleus
            G = [G; EfficiencyFRET(Lam >= Div(q) & Lam < Div(q+1))];
            grp = [grp;q*ones(length(EfficiencyFRET(Lam >= Div(q) & Lam < Div(q+1))), 1)];
        end
        
        AlleleG = [AlleleG; G];
        geneID = [geneID; grp];
    end

    figure('units','normalized','outerposition',[0 0 0.2 0.3])
    boxplot(AlleleG, geneID, 'labels', {'0 - 25%', '25 - 50%', '50 - 75%','75 - 100%'},'Symbol', '.k', 'Widths',0.7)
    hold on
    for n = 1:3
        % Calculate and plot the p-value
        plot(n:n+1, [50-n, 50-n], '-k', 'Linewidth', 1)
        A = AlleleG(geneID == n);
        B = AlleleG(geneID == n+1);
        text((n+n+1)/2, 55-n, ['P = ' num2str(round(ranksum(A(~isnan(A)),...
            B(~isnan(B))), 2,'significant'))],'HorizontalAlignment', 'center','FontSize', 8)
        text(n, 10, ['n = ' num2str(sum(~isnan(A)))],'HorizontalAlignment', 'center','FontSize', 8)
    end
    text(n+1, 10, ['n = ' num2str(sum(~isnan(B)))],'HorizontalAlignment', 'center','FontSize', 8)
    ylabel('FRET efficiency (%)')
    xlabel('Lamina distance quantiles (pixel)')
    title(gene{g})
end