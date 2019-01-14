clearvars -except total nb_sets
close all
clc

cd(fileparts(which(mfilename)));

sett = load_settings();

MIN_ENER_PHOT = 5; % keV
MAX_ENER_PHOT = 100000;
MIN_ENER_ELEC = 5;
MAX_ENER_ELEC = 100000;

% ! sh COPY_FINAL_DATAFILES_FROM_FRAM.sh

loaded = load([sett.base_path 'BIG_DATAFILE_all.mat']);
BIG_DATAFILE = loaded.BIG_DATAFILE;

% [status,cmdout] = system(['date -r ' [sett.base_path 'BIG_DATAFILE_all.mat'] ' "+%m-%d-%Y %H:%M:%S"']);
% disp(cmdout)

% pause(2)

POTENTIAL_LIST = sett.POTENTIAL_LIST;
RECORD_POS_LIST = sett.RECORD_POS_LIST;
EFIELD_CENTER_list = sett.EFIELD_CENTER_list;
EFIELD_SIZE_list = sett.EFIELD_SIZE_list;

i_0MV = find_idx(POTENTIAL_LIST,0);

glow_database = cell(...
    length(sett.record_type_PDG_NB),...
    length(POTENTIAL_LIST),...
    length(RECORD_POS_LIST),...
    length(EFIELD_CENTER_list),...
    length(EFIELD_SIZE_list));

%%
for i_recPos = 1 : length(RECORD_POS_LIST)
    for i_efield_c = 1 : length(EFIELD_CENTER_list)
        for i_efield_s = 2
            
            nb_files=[];
            
            ALT = EFIELD_CENTER_list(i_efield_c);
            RECORD_POS = RECORD_POS_LIST(i_recPos);
            EFIELD_SIZE = EFIELD_SIZE_list(i_efield_s);
            
            rec_alt = ALT+RECORD_POS*EFIELD_SIZE/2.0;
            
            if ~ismember(rec_alt,sett.WANTED_RECORD_ALTS)
                continue;
            end
            
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
                
                nb_recorded(ii,:) = (vals) ./ ref;
            end
            
            disp(' ')
            disp(num2str(nb_recorded))
            
            disp(' ')
            disp(num2str(nb_files))
            
            %%
            
            subplot(1,2,1)
            plot(POTENTIAL_LIST,nb_recorded(1,:))
            hold on
            plot(POTENTIAL_LIST,nb_recorded(2,:))
            plot(POTENTIAL_LIST,nb_recorded(3,:))
            hold off
            legend('photon','electron','positron')
            grid on
            xlabel('potential (MV)')
            ylabel('Change multiplicaton factor')
            rec_alt = ALT+RECORD_POS*EFIELD_SIZE/2.0;
            title(['ALT = ' num2str(ALT) ' km; EFIELD SIZE = ' num2str(EFIELD_SIZE) ' km; RECORD POS = ' num2str(RECORD_POS) '   (=' num2str(rec_alt) ' km)'])
            set(gca,'yscale','log')
            
            %% build glow database
            
            type_names = sett.record_names;
            colors_type = {'r','g','b'};
            
            for i_pot = 1:length(POTENTIAL_LIST)
                subplot(1,2,2)
                hold off
                
                failed=false;
                
                for i_t = 1:3
                    try
                        BINS = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID;
                        ENERGY_HIST = BIG_DATAFILE.(type_names{i_t}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);
                        
                        ENERGY_HIST_0MV = BIG_DATAFILE.(type_names{i_t}){i_0MV,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(1:end-1);

                        histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST,'displayStyle','stairs','LineWidth',2,'edgecolor',colors_type{i_t})
                        hold on
                        histogram('BinEdges',BINS,'BinCounts',ENERGY_HIST_0MV,'displayStyle','stairs','LineWidth',1,'edgecolor','k','LineStyle','--','HandleVisibility','off')
                        
                        glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST = ENERGY_HIST(1:end-1);
                        glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV = ENERGY_HIST_0MV(1:end-1);
                        glow_database{i_t,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID = BINS(1:end-1).*1000; % MeV to keV
                    catch ME
                        disp(['type: ' num2str(i_t) ' ; failed to plot'])
                        failed=true;
                    end
                end
                
                if ~failed
                    set(gca,'xscale','log')
                    set(gca,'yscale','log')
                    set(gca,'xlim',[0.005 200])
                    xlabel('Energy (MeV)')
                    ylabel('Counts per bin (a.u.)')
                    grid on
                    title(['ALT = ' num2str(ALT) ' km; EFIELD SIZE = ' num2str(EFIELD_SIZE) ' km; ' ...
                        'RECORD POS = ' num2str(RECORD_POS) '   (=' num2str(rec_alt) ' km)  ; ' ...
                        'POTENTIAL = ' num2str(POTENTIAL_LIST(i_pot)) ])
                    
                    legend('photon','electron','positron')
                    
                    if usejava('jvm')
%                         pause(0.1)
                        filename = [num2str(ALT) '_' num2str(EFIELD_SIZE) ...
                            '_' num2str(RECORD_POS) '_' num2str(rec_alt) '_km_' num2str(POTENTIAL_LIST(i_pot)) '.png' ];
                        saveas(gcf,['./PLOTS/' filename])
                    end
                    
                end
            end
            
        end
    end
end


%%

save([sett.base_path 'glow_database.mat'],'glow_database');

disp(' ')
disp('DONE.')
close all



%%
function NB_RECORDED = get_NB_RECORDED_in_energy_range(data_struct,min_ener,max_ener)

sett=load_settings();

grid = data_struct.ENERGY_GRID*1000.0;

% data_struct.ENERGY_HIST(data_struct.ENERGY_HIST<10^-20) = sett.MIN_HIST_BCKGRND;
% data_struct.ENERGY_HIST(data_struct.ENERGY_HIST>10^20) = sett.MIN_HIST_BCKGRND;

eh = data_struct.ENERGY_HIST;

NB_RECORDED = sum(eh(grid>min_ener & grid<max_ener))*1e5;

end