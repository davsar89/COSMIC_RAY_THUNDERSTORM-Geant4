clearvars -except total nb_sets
close all
clc

% addpath('./maxdistcolor/')
% fun = @CIELab_to_DIN99o;
% rgb = maxdistcolor(9,fun);

plot_x_min = 0.026;
plot_x_max = 100;

fontsize_names_types = 15;

linew_backgrnd = 3;

RREA_ALT = [13 14 15];
RREA_THRES = [108.5 125.6 143.6];

figure('WindowStyle','normal','Position', [0 0 800 600*2])
tiledlayout(4,2, 'Padding', 'tight', 'TileSpacing', 'tight');

%%
sett = load_settings();

MIN_ENER_PHOT = 5; % keV
MAX_ENER_PHOT = 1000000;
MIN_ENER_ELEC = 5;
MAX_ENER_ELEC = 1000000;

WANTED_EFIELD_ALT_CENTER = 14;
WANTED_EFIELD_FULL_SIZE = 2;
WANTED_REC_POS = 1;

idx_list_plot_spec = [1 2 3 5 6 7 13 18 19 20 21 23 25];

% ! sh COPY_FINAL_DATAFILES_FROM_FRAM.sh

loaded = load([sett.base_path 'BIG_DATAFILE_all.mat']);
BIG_DATAFILE = loaded.BIG_DATAFILE;

[status,cmdout] = system(['date -r ' [sett.base_path 'BIG_DATAFILE_all.mat'] ' "+%m-%d-%Y %H:%M:%S"']);
disp(cmdout)

POTENTIAL_LIST = sett.POTENTIAL_LIST;
RECORD_POS_LIST = sett.RECORD_POS_LIST;
EFIELD_CENTER_list = sett.EFIELD_CENTER_list;
EFIELD_SIZE_list = sett.EFIELD_SIZE_list;

i_0MV = find_idx(POTENTIAL_LIST,0);

i_efield_c = find_idx(EFIELD_CENTER_list,WANTED_EFIELD_ALT_CENTER);
i_recPos = find_idx(RECORD_POS_LIST,WANTED_REC_POS);
i_efield_s = find_idx(EFIELD_SIZE_list,WANTED_EFIELD_FULL_SIZE);

glow_database = cell(...
    length(sett.record_type_PDG_NB),...
    length(POTENTIAL_LIST),...
    length(RECORD_POS_LIST),...
    length(EFIELD_CENTER_list),...
    length(EFIELD_SIZE_list));

%%
nb_files=[];

ALT = EFIELD_CENTER_list(i_efield_c);
RECORD_POS = RECORD_POS_LIST(i_recPos);
EFIELD_SIZE = EFIELD_SIZE_list(i_efield_s);

rec_alt = ALT+RECORD_POS*EFIELD_SIZE/2.0;

%%

for i_pot = 1:length(POTENTIAL_LIST)
    
    POTENTIAL = POTENTIAL_LIST(i_pot);
    
    if ~isempty(BIG_DATAFILE.photon{i_pot,i_recPos,i_efield_c,i_efield_s})
        
        nb_recorded(1,i_pot) = ...
            get_NB_RECORDED_in_energy_range(BIG_DATAFILE.photon{i_pot,i_recPos,i_efield_c,i_efield_s},...
            MIN_ENER_PHOT,...
            MAX_ENER_PHOT);
        nb_recorded(2,i_pot) = ...
            get_NB_RECORDED_in_energy_range(BIG_DATAFILE.electron{i_pot,i_recPos,i_efield_c,i_efield_s},...
            MIN_ENER_ELEC,...
            MAX_ENER_ELEC);
        nb_recorded(3,i_pot) = ...
            get_NB_RECORDED_in_energy_range(BIG_DATAFILE.positron{i_pot,i_recPos,i_efield_c,i_efield_s},...
            MIN_ENER_ELEC,...
            MAX_ENER_ELEC);
        
        nb_files(1,i_pot) = BIG_DATAFILE.photon{i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES;
        nb_files(2,i_pot) = BIG_DATAFILE.electron{i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES;
        nb_files(3,i_pot) = BIG_DATAFILE.positron{i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES;
    else
        
        nb_recorded(1,i_pot) = 0;
        nb_recorded(2,i_pot) = 0;
        nb_recorded(3,i_pot) = 0;
    end
end

nb_recorded_before = nb_recorded;

for ii=1:3
    
    ref = nb_recorded(ii,i_0MV);
    if ref==0
        ref = interp1([-10 10],nb_recorded_before(ii,[i_0MV-1 i_0MV+1]),0);
    end
    
    vals = nb_recorded(ii,:);
    
    nb_recorded(ii,:) = (vals) ./ ref ;
end

disp(' ')
disp(num2str(nb_recorded))

disp(' ')
disp(num2str(nb_files))

%%

nexttile([1 2])
% colororder(rgb)
% colororder([[0,0.45,0.74];...
%     [0.85,0.33,0.10];...
%     [0.93,0.69,0.13];...
%     [0.49,0.18,0.56];...
%     [0.30,0.75,0.93]...
%     ])
plot(POTENTIAL_LIST,nb_recorded(1,:),'linewidth',2,'color','k')
hold on
plot(POTENTIAL_LIST,nb_recorded(2,:),'linewidth',2,'color','b')
plot(POTENTIAL_LIST,nb_recorded(3,:),'linewidth',2,'color','r')

%%
NAMES_COL = {'Potential (MV)', 'photon increase factor', 'electron increase factor', 'positron increase factor'};
make_datafile_4vals(POTENTIAL_LIST, nb_recorded(1,:), nb_recorded(2,:), nb_recorded(3,:), NAMES_COL, 'relative_increase')
%%

legend('photon','electron','positron','location','north','fontsize',8)
ax = gca;
ax.YGrid = 'on';
ax.YMinorGrid = 'off';
ax.XGrid = 'on';
ax.XMinorGrid = 'off';
xlabel('Potential $\Delta U$ (MV)','interpreter','latex','fontsize',14)
ylabel('Multiplication factor $F_{inc}$','interpreter','latex','fontsize',12)
rec_alt = ALT+RECORD_POS*EFIELD_SIZE/2.0;
% title(['ALT = ' num2str(ALT) ' km; EFIELD SIZE = ' num2str(EFIELD_SIZE) ' km; RECORD POS = ' num2str(RECORD_POS) '   (=' num2str(rec_alt) ' km)'])
set(gca,'yscale','log')
set(gca, 'YLim', [0.6 6e3]);

set(gca, 'XLimSpec', 'Tight');

text(0.02,0.92,'a.','Units','normalized','fontweight','bold','FontSize',25)

plot([RREA_THRES(2) RREA_THRES(2)],get(gca,'ylim'),'k--','HandleVisibility','off')
plot([-RREA_THRES(2) -RREA_THRES(2)],get(gca,'ylim'),'k--','HandleVisibility','off')

plot([RREA_THRES(1) RREA_THRES(1)],get(gca,'ylim'),'--','color',[0.6,0.6,0.6],'HandleVisibility','off')
plot([RREA_THRES(3) RREA_THRES(3)],get(gca,'ylim'),'--','color',[0.6,0.6,0.6],'HandleVisibility','off')
plot([-RREA_THRES(1) -RREA_THRES(1)],get(gca,'ylim'),'--','color',[0.6,0.6,0.6],'HandleVisibility','off')
plot([-RREA_THRES(3) -RREA_THRES(3)],get(gca,'ylim'),'--','color',[0.6,0.6,0.6],'HandleVisibility','off')

myText=text(0.205,0.6,'RREA','Units','normalized','fontweight','bold','FontSize',10,'color',[0.3,0.3,0.3]);
set(myText,'Rotation',90);
myText=text(0.775,0.6,'RREA','Units','normalized','fontweight','bold','FontSize',10,'color',[0.3,0.3,0.3]);
set(myText,'Rotation',90);

yticks([1 10 100 1000 6000])
yticklabels({'10^0','10^1','10^2','10^3','6.10^3'})

hold off

%% plot glow photon NEGATIVE
nexttile
type_names = sett.record_names;
colors_type = {'r','g','b'};

for i_pot = idx_list_plot_spec
    
    if sett.POTENTIAL_LIST(i_pot)>0
        continue
    end
    
    failed=false;
    
    for i_t = 1:1
        try
            BINS = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID;
            ENERGY_HIST = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
            ENERGY_HIST_0MV = BIG_DATAFILE.(type_names{i_t}){i_0MV,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
            %                         pc_fact(i_t) = 1.0+nb_recorded(i_t,i_pot)/100.0;
            %                         ENERGY_HIST = ENERGY_HIST ./ sum(ENERGY_HIST) .* pc_fact(i_t);
            %                         ENERGY_HIST_0MV = ENERGY_HIST_0MV ./sum(ENERGY_HIST_0MV);
            hold on
            if sett.POTENTIAL_LIST(i_pot)==0
                hh = histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST./diff(BINS),'displayStyle','stairs','LineWidth',linew_backgrnd,'edgecolor','k','LineStyle','-.');
            else
                histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST./diff(BINS),'displayStyle','stairs','LineWidth',2)
            end
            
            %%
            NAMES_COL = {'Energy bin edges (MeV)', 'dn/de spectrum'};
            make_datafile_1vals(BINS, ENERGY_HIST./diff(BINS), NAMES_COL, ['photon_spectrun_' num2str(sett.POTENTIAL_LIST(i_pot)) 'MV' ])
            %%
            
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST = ENERGY_HIST;
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV = ENERGY_HIST_0MV;
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID = BINS;
        catch ME
            disp(['type: ' num2str(i_t) ' ; failed to plot'])
            failed=true;
        end
    end
    
    if ~failed
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        set(gca, 'YLimSpec', 'Tight');
        set(gca,'xlim',[plot_x_min plot_x_max])
        set(gca,'ylim',[1e-6 70])
        %         xlabel('Energy (MeV)')
        ylabel('n/de spectrum (a.u.)','interpreter','latex','fontsize',12)
        ax = gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'off';
        ax.XGrid = 'on';
        ax.XMinorGrid = 'off';
%         grid on
        
        %         legend('photon','electron','positron')
    end
end

legend_names={};
for i_pot = idx_list_plot_spec
    if sett.POTENTIAL_LIST(i_pot)>0
        continue
    end
    legend_names{end+1} = ['\DeltaU = ' num2str(sett.POTENTIAL_LIST(i_pot)) ' MV'];
end

legend(legend_names,'location','northoutside','NumColumns',3,'FontSize',7)
yticks([1.e-6 1.e-5 1.e-4 1.e-3 1.e-2 1.e-1 1.e-0 1.e1 1.e2])
yticklabels({'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}'})
text(0.3,0.90,'b.','Units','normalized','fontweight','bold','FontSize',25)
text(0.10,0.10,'Photon $\Delta U \leq$ 0','Units','normalized','FontSize',fontsize_names_types,'color',[0.3,0.3,0.3],'interpreter','latex','fontsize',14)
uistack(hh,'top');
box on;

%% plot glow photon POSITIVE
nexttile
type_names = sett.record_names;
colors_type = {'r','g','b'};

for i_pot = idx_list_plot_spec
    
    if sett.POTENTIAL_LIST(i_pot)<0
        continue
    end
    
    failed=false;
    
    for i_t = 1:1
        try
            BINS = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID;
            ENERGY_HIST = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
            ENERGY_HIST_0MV = BIG_DATAFILE.(type_names{i_t}){i_0MV,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
            %                         pc_fact(i_t) = 1.0+nb_recorded(i_t,i_pot)/100.0;
            %                         ENERGY_HIST = ENERGY_HIST ./ sum(ENERGY_HIST) .* pc_fact(i_t);
            %                         ENERGY_HIST_0MV = ENERGY_HIST_0MV ./sum(ENERGY_HIST_0MV);
            hold on
            if sett.POTENTIAL_LIST(i_pot)==0
                hh = histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST./diff(BINS),'displayStyle','stairs','LineWidth',linew_backgrnd,'edgecolor','k','LineStyle','-.');
            else
                histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST./diff(BINS),'displayStyle','stairs','LineWidth',2)
            end
            
            %%
            NAMES_COL = {'Energy bin edges (MeV)', 'dn/de spectrum'};
            make_datafile_1vals(BINS, ENERGY_HIST./diff(BINS), NAMES_COL, ['photon_spectrun_' num2str(sett.POTENTIAL_LIST(i_pot)) 'MV' ])
            %%
            
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST = ENERGY_HIST;
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV = ENERGY_HIST_0MV;
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID = BINS;
        catch ME
            disp(['type: ' num2str(i_t) ' ; failed to plot'])
            failed=true;
        end
    end
    
    if ~failed
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        set(gca, 'YLimSpec', 'Tight');
        set(gca,'xlim',[plot_x_min plot_x_max])
        set(gca,'ylim',[1e-6 70])
        %         xlabel('Energy (MeV)')
        %         ylabel('n/de spectrum (a.u.)','interpreter','latex','fontsize',12)
        ax=gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'off';
        ax.XGrid = 'on';
        ax.XMinorGrid = 'off';
        
        %         legend('photon','electron','positron')
    end
end

legend_names={};
for i_pot = idx_list_plot_spec
    if sett.POTENTIAL_LIST(i_pot)<0
        continue
    end
    legend_names{end+1} = ['\DeltaU = ' num2str(sett.POTENTIAL_LIST(i_pot)) ' MV'];
end

legend(legend_names,'location','northoutside','NumColumns',3,'FontSize',7)
yticks([1.e-6 1.e-5 1.e-4 1.e-3 1.e-2 1.e-1 1.e-0 1.e1 1.e2])
yticklabels({'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}'})
text(0.5,0.95,'c.','Units','normalized','fontweight','bold','FontSize',25)
text(0.10,0.10,'Photon $\Delta U \geq$ 0','Units','normalized','FontSize',fontsize_names_types,'color',[0.3,0.3,0.3],'interpreter','latex','fontsize',14)
uistack(hh,'top');
box on;

%% plot glow electron NEGATIVE
nexttile
type_names = sett.record_names;
colors_type = {'r','g','b'};

for i_pot = idx_list_plot_spec
    
    if sett.POTENTIAL_LIST(i_pot)>0
        continue
    end
    
    failed=false;
    
    for i_t = 2:2
        try
            BINS = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID;
            ENERGY_HIST = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
            ENERGY_HIST_0MV = BIG_DATAFILE.(type_names{i_t}){i_0MV,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
            %                         pc_fact(i_t) = 1.0+nb_recorded(i_t,i_pot)/100.0;
            %                         ENERGY_HIST = ENERGY_HIST ./ sum(ENERGY_HIST) .* pc_fact(i_t);
            %                         ENERGY_HIST_0MV = ENERGY_HIST_0MV ./sum(ENERGY_HIST_0MV);
            hold on
            if sett.POTENTIAL_LIST(i_pot)==0
                hh = histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST./diff(BINS),'displayStyle','stairs','LineWidth',linew_backgrnd,'edgecolor','k','LineStyle','-.');
            else
                histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST./diff(BINS),'displayStyle','stairs','LineWidth',2)
            end
            
            %%
            NAMES_COL = {'Energy bin edges (MeV)', 'dn/de spectrum'};
            make_datafile_1vals(BINS, ENERGY_HIST./diff(BINS), NAMES_COL, ['electron_spectrun_' num2str(sett.POTENTIAL_LIST(i_pot)) 'MV' ])
            %%
            
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST = ENERGY_HIST;
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV = ENERGY_HIST_0MV;
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID = BINS;
        catch ME
            disp(['type: ' num2str(i_t) ' ; failed to plot'])
            failed=true;
        end
    end
    
    if ~failed
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        set(gca, 'YLimSpec', 'Tight');
        set(gca,'xlim',[plot_x_min plot_x_max])
        set(gca,'ylim',[0.5e-6 0.15])
        %         xlabel('Energy (MeV)')
        ylabel('n/de spectrum (a.u.)','interpreter','latex','fontsize',12)
        ax=gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'off';
        ax.XGrid = 'on';
        ax.XMinorGrid = 'off';
        
        %         legend('photon','electron','positron')
    end
end

legend_names={};
for i_pot = idx_list_plot_spec
    if sett.POTENTIAL_LIST(i_pot)>0
        continue
    end
    legend_names{end+1} = ['\DeltaU = ' num2str(sett.POTENTIAL_LIST(i_pot)) ' MV'];
end

legend(legend_names,'location','northoutside','NumColumns',3,'FontSize',7)
yticks([1.e-6 1.e-5 1.e-4 1.e-3 1.e-2 1.e-1])
yticklabels({'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'})

text(0.02,0.90,'d.','Units','normalized','fontweight','bold','FontSize',25)
text(0.20,0.90,'Electron $\Delta U \leq$ 0','Units','normalized','FontSize',fontsize_names_types,'color',[0.3,0.3,0.3],'interpreter','latex','fontsize',14)
uistack(hh,'top');
box on;

%% plot glow electron POSITIVE
nexttile
type_names = sett.record_names;
colors_type = {'r','g','b'};

for i_pot = idx_list_plot_spec
    
    if sett.POTENTIAL_LIST(i_pot)<0
        continue
    end
    
    failed=false;
    
    for i_t = 2:2
        try
            BINS = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID;
            ENERGY_HIST = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
            ENERGY_HIST_0MV = BIG_DATAFILE.(type_names{i_t}){i_0MV,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
            %                         pc_fact(i_t) = 1.0+nb_recorded(i_t,i_pot)/100.0;
            %                         ENERGY_HIST = ENERGY_HIST ./ sum(ENERGY_HIST) .* pc_fact(i_t);
            %                         ENERGY_HIST_0MV = ENERGY_HIST_0MV ./sum(ENERGY_HIST_0MV);
            hold on
            if sett.POTENTIAL_LIST(i_pot)==0
                hh = histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST./diff(BINS),'displayStyle','stairs','LineWidth',linew_backgrnd,'edgecolor','k','LineStyle','-.');
            else
                histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST./diff(BINS),'displayStyle','stairs','LineWidth',2)
            end
            
            %%
            NAMES_COL = {'Energy bin edges (MeV)', 'dn/de spectrum'};
            make_datafile_1vals(BINS, ENERGY_HIST./diff(BINS), NAMES_COL, ['electron_spectrun_' num2str(sett.POTENTIAL_LIST(i_pot)) 'MV' ])
            %%
            
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST = ENERGY_HIST;
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV = ENERGY_HIST_0MV;
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID = BINS;
        catch ME
            disp(['type: ' num2str(i_t) ' ; failed to plot'])
            failed=true;
        end
    end
    
    if ~failed
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        set(gca, 'YLimSpec', 'Tight');
        set(gca,'xlim',[plot_x_min plot_x_max])
        set(gca,'ylim',[0.5e-6 0.15])
        %         xlabel('Energy (MeV)')
        %         ylabel('n/de spectrum (a.u.)','interpreter','latex','fontsize',12)
        ax=gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'off';
        ax.XGrid = 'on';
        ax.XMinorGrid = 'off';
        
        %         legend('photon','electron','positron')
    end
end

legend_names={};
for i_pot = idx_list_plot_spec
    if sett.POTENTIAL_LIST(i_pot)<0
        continue
    end
    legend_names{end+1} = ['\DeltaU = ' num2str(sett.POTENTIAL_LIST(i_pot)) ' MV'];
end

legend(legend_names,'location','northoutside','NumColumns',3,'FontSize',7)
yticks([1.e-6 1.e-5 1.e-4 1.e-3 1.e-2 1.e-1])
yticklabels({'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'})
text(0.90,0.95,'e.','Units','normalized','fontweight','bold','FontSize',25)
text(0.10,0.10,'Electron $\Delta U \geq$ 0','Units','normalized','FontSize',fontsize_names_types,'color',[0.3,0.3,0.3],'interpreter','latex','fontsize',14)
uistack(hh,'top');
box on;

%% plot glow positron NEGATIVE
nexttile
type_names = sett.record_names;
colors_type = {'r','g','b'};

for i_pot = idx_list_plot_spec
    
    if sett.POTENTIAL_LIST(i_pot)>0
        continue
    end
    
    failed=false;
    
    for i_t = 3:3
        try
            BINS = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID;
            ENERGY_HIST = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
            ENERGY_HIST_0MV = BIG_DATAFILE.(type_names{i_t}){i_0MV,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
            %                         pc_fact(i_t) = 1.0+nb_recorded(i_t,i_pot)/100.0;
            %                         ENERGY_HIST = ENERGY_HIST ./ sum(ENERGY_HIST) .* pc_fact(i_t);
            %                         ENERGY_HIST_0MV = ENERGY_HIST_0MV ./sum(ENERGY_HIST_0MV);
            hold on
            if sett.POTENTIAL_LIST(i_pot)==0
                hh = histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST./diff(BINS),'displayStyle','stairs','LineWidth',linew_backgrnd,'edgecolor','k','LineStyle','-.');
            else
                histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST./diff(BINS),'displayStyle','stairs','LineWidth',2)
            end
            
            %%
            NAMES_COL = {'Energy bin edges (MeV)', 'dn/de spectrum'};
            make_datafile_1vals(BINS, ENERGY_HIST./diff(BINS), NAMES_COL, ['positron_spectrun_' num2str(sett.POTENTIAL_LIST(i_pot)) 'MV' ])
            %%
            
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST = ENERGY_HIST;
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV = ENERGY_HIST_0MV;
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID = BINS;
        catch ME
            disp(['type: ' num2str(i_t) ' ; failed to plot'])
            failed=true;
        end
    end
    
    if ~failed
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        set(gca, 'YLimSpec', 'Tight');
        set(gca,'xlim',[plot_x_min plot_x_max])
        set(gca,'ylim',[0.25e-6 4.e-4])
        xlabel('Energy (MeV)','interpreter','latex','fontsize',14)
        ylabel('n/de spectrum (a.u.)','interpreter','latex','fontsize',12)
        ax=gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'off';
        ax.XGrid = 'on';
        ax.XMinorGrid = 'off';
        
        %         legend('photon','electron','positron')
    end
end

legend_names={};
for i_pot = idx_list_plot_spec
    if sett.POTENTIAL_LIST(i_pot)>0
        continue
    end
    legend_names{end+1} = ['\DeltaU = ' num2str(sett.POTENTIAL_LIST(i_pot)) ' MV'];
end

legend(legend_names,'location','northoutside','NumColumns',3,'FontSize',7)
yticks([1.e-6 1.e-5 1.e-4 1.e-3 1.e-2 1.e-1])
yticklabels({'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'})
text(0.5,0.15,'f.','Units','normalized','fontweight','bold','FontSize',25)
text(0.20,0.90,'Positron $\Delta U$ $\leq$ 0','Units','normalized','FontSize',fontsize_names_types,'color',[0.3,0.3,0.3],'interpreter','latex','fontsize',14)
uistack(hh,'top');
box on;

%% plot glow positron POSITIVE
nexttile
type_names = sett.record_names;
colors_type = {'r','g','b'};

for i_pot = idx_list_plot_spec
    
    if sett.POTENTIAL_LIST(i_pot)<0
        continue
    end
    
    failed=false;
    
    for i_t = 3:3
        try
            BINS = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID;
            ENERGY_HIST = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
            ENERGY_HIST_0MV = BIG_DATAFILE.(type_names{i_t}){i_0MV,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
            %                         pc_fact(i_t) = 1.0+nb_recorded(i_t,i_pot)/100.0;
            %                         ENERGY_HIST = ENERGY_HIST ./ sum(ENERGY_HIST) .* pc_fact(i_t);
            %                         ENERGY_HIST_0MV = ENERGY_HIST_0MV ./sum(ENERGY_HIST_0MV);
            hold on
            if sett.POTENTIAL_LIST(i_pot)==0
                hh = histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST./diff(BINS),'displayStyle','stairs','LineWidth',linew_backgrnd,'edgecolor','k','LineStyle','-.');
            else
                histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST./diff(BINS),'displayStyle','stairs','LineWidth',2)
            end
            
            %%
            NAMES_COL = {'Energy bin edges (MeV)', 'dn/de spectrum'};
            make_datafile_1vals(BINS, ENERGY_HIST./diff(BINS), NAMES_COL, ['positron_spectrun_' num2str(sett.POTENTIAL_LIST(i_pot)) 'MV' ])
            %%
            
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST = ENERGY_HIST;
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV = ENERGY_HIST_0MV;
            glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID = BINS;
        catch ME
            disp(['type: ' num2str(i_t) ' ; failed to plot'])
            failed=true;
        end
    end
    
    if ~failed
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        set(gca, 'YLimSpec', 'Tight');
        set(gca,'xlim',[plot_x_min plot_x_max])
        set(gca,'ylim',[0.25e-6 4.e-4])
        xlabel('Energy (MeV)','interpreter','latex','fontsize',14)
        %         ylabel('n/de spectrum (a.u.)','interpreter','latex','fontsize',12)
        ax=gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'off';
        ax.XGrid = 'on';
        ax.XMinorGrid = 'off';
        
        %         legend('photon','electron','positron')
    end
end

legend_names={};
for i_pot = idx_list_plot_spec
    if sett.POTENTIAL_LIST(i_pot)<0
        continue
    end
    legend_names{end+1} = ['\DeltaU = ' num2str(sett.POTENTIAL_LIST(i_pot)) ' MV'];
end

legend(legend_names,'location','northoutside','NumColumns',3,'FontSize',7)
yticks([1.e-6 1.e-5 1.e-4 1.e-3 1.e-2 1.e-1])
yticklabels({'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'})
text(0.90,0.6,'g.','Units','normalized','fontweight','bold','FontSize',25)
text(0.6,0.91,'Positron $\Delta U \geq$ 0','Units','normalized','FontSize',fontsize_names_types,'color',[0.3,0.3,0.3],'interpreter','latex','fontsize',14)
uistack(hh,'top');
box on;

%%
sgtitle(['$H_E$ = ' num2str(ALT) ' km; $\Delta H_E$ = ' num2str(EFIELD_SIZE) ' km; ' ...
    '$H_D$ = ' num2str(rec_alt) ' km, instr. resp. not included' ],'interpreter','latex','fontsize',12)

copygraphics(gcf)

saveas(gcf,'/Home/siv29/dsa030/Desktop/article_cr_modeling/figures/COMPA_FLUX_SPEC_EXAMPLE.eps','epsc')
saveas(gcf,'/Home/siv29/dsa030/Desktop/article_cr_modeling/figures/COMPA_FLUX_SPEC_EXAMPLE.fig')

%%
function NB_RECORDED = get_NB_RECORDED_in_energy_range(data_struct,min_ener,max_ener)

grid = data_struct.ENERGY_GRID*1000.0;

eh = data_struct.ENERGY_HIST;

NB_RECORDED = sum(eh(grid>min_ener & grid<max_ener))*1e5;

end


%% 