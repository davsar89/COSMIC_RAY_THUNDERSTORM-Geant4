close all
clc
clear all

cd(fileparts(which(mfilename)));

sett=load_settings();

POTENTIAL_LIST = sett.POTENTIAL_LIST;
RECORD_POS_LIST = sett.RECORD_POS_LIST;
EFIELD_CENTER_list = sett.EFIELD_CENTER_list;
EFIELD_SIZE_list = sett.EFIELD_SIZE_list;

IS_FRAM = sett.IS_FRAM;
base_path = sett.base_path;

record_names={'photon','electron','positron'};

BIG_DATAFILE=[];

BIG_DATAFILE.photon = cell(length(POTENTIAL_LIST),length(RECORD_POS_LIST),...
    length(EFIELD_CENTER_list),length(EFIELD_SIZE_list));
BIG_DATAFILE.electron = cell(length(POTENTIAL_LIST),length(RECORD_POS_LIST),...
    length(EFIELD_CENTER_list),length(EFIELD_SIZE_list));
BIG_DATAFILE.positron = cell(length(POTENTIAL_LIST),length(RECORD_POS_LIST),...
    length(EFIELD_CENTER_list),length(EFIELD_SIZE_list));

%%
for i_pot = 1:length(POTENTIAL_LIST)
    for i_recPos = 1:length(RECORD_POS_LIST)
        for i_efield_c = 1:length(EFIELD_CENTER_list)
            for i_efield_s = 1:length(EFIELD_SIZE_list)
                
                POTENTIAL = POTENTIAL_LIST(i_pot);
                SELECTED_RECORD_POSITION = RECORD_POS_LIST(i_recPos);
                EFIELD_CENTER = EFIELD_CENTER_list(i_efield_c);
                EFIELD_SIZE = EFIELD_SIZE_list(i_efield_s);
                
                rec_alt = EFIELD_CENTER + SELECTED_RECORD_POSITION*EFIELD_SIZE/2.0;
                
                if ~ismember(rec_alt,sett.WANTED_RECORD_ALTS)
                    continue
                end
                
                folder_simu_results_for_param_set = [ base_path num2str(POTENTIAL) 'MV/' num2str(SELECTED_RECORD_POSITION) '/' ...
                    num2str(EFIELD_CENTER) 'km_' num2str(EFIELD_SIZE) 'km/'];
                
                if ~exist(folder_simu_results_for_param_set, 'dir')
                    continue
                end
                
                disp(['POTENTIAL: ' num2str(POTENTIAL)])
                disp(['REC POS: ' num2str(SELECTED_RECORD_POSITION)])
                disp(['EFIELD CENTER: ' num2str(EFIELD_CENTER)])
                disp(['EFIELD FULL SIZE: ' num2str(EFIELD_SIZE)])
                
                data_files_for_run = dir([folder_simu_results_for_param_set '*/ALL_ener_mom_dists*.out']);
                
                counter = 0;
                
                for i_file = 1:length(data_files_for_run)
                    
                    filename = [data_files_for_run(i_file).folder '/' data_files_for_run(i_file).name];
                    ok = false;
                    PARAMETER_SET = [POTENTIAL SELECTED_RECORD_POSITION EFIELD_SIZE EFIELD_CENTER];
                    [BIG_DATAFILE,ok] = add_to_BIG_DATAFILE(filename,BIG_DATAFILE,i_pot,i_recPos,i_efield_c,i_efield_s,PARAMETER_SET);
                    if ok
                        counter = counter +1;
                    end
                end
                disp(['NB FILES: ' num2str(counter)])
            end
        end
    end
end

% normalize with the number of files
for irt=1:length(record_names)
    for i_pot =1:length(POTENTIAL_LIST)
        for i_efield_s = 1:length(EFIELD_SIZE_list)
            for i_efield_c = 1:length(EFIELD_CENTER_list)
                for i_recPos = 1:length(RECORD_POS_LIST)
                    
                    if isfield(BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s},'NB_FILES')
                        
                        nb_f = double(BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES);
                        
                        BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.ANGLE_ENERGY_HIST = ...
                            BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.ANGLE_ENERGY_HIST./nb_f;
                        
                        BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.TOTAL_NB_RECORDED = ...
                            BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.TOTAL_NB_RECORDED./nb_f;
                        
                        %% add the energy hisogram (sum for all zenith angles energy histograms)
                        BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST = ...
                            sum(BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.ANGLE_ENERGY_HIST,1);
                        
%                         BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST<10^-20) = sett.MIN_HIST_BCKGRND;
%                         BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST(BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST>10^20) = sett.MIN_HIST_BCKGRND;
%                         
                        BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.EGRID_UNIT = 'Energy grid unit is MeV';
                        BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.comment1 = 'already taken into account that the number of files to build the data base may not be the same for all input parameters (by dividing by number of file for given parameter set).';
                        BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.comment2 = 'last energy histogram bin values is overflow (> 100 MeV)';
                        BIG_DATAFILE.(record_names{irt}){i_pot,i_recPos,i_efield_c,i_efield_s}.comment3 = 'histograms are in unit of counts per bin (NO division by bin size). Scaling is arbitrary but consistent.';
                    end
                    
                end
            end
        end
    end
end

%% save to .mat files

if IS_FRAM
    save('BIG_DATAFILE_all.mat','BIG_DATAFILE','-v7.3')
else
    save([sett.base_path 'BIG_DATAFILE_all.mat'],'BIG_DATAFILE','-v7.3')
end

%% plot for sanity check
% if ~IS_FRAM
%     BINS = load_energy_grid(BIG_DATAFILE.photon{1,6,7,2});
%     full_hist = sum(BIG_DATAFILE.photon{1,6,7,2}.ANGLE_ENERGY_HIST,1);
%     disp(['nb files used: ' num2str(BIG_DATAFILE.photon{1,6,7,2}.NB_FILES)])
%     histogram('BinEdges',BINS,'BinCounts',full_hist(1:end-1),'displayStyle','stairs','LineWidth',2)
%     set(gca,'xscale','log')
%     set(gca,'yscale','log')
%     xlabel('Energy (MeV)')
%     ylabel('Counts per bin (a.u.)')
%     grid on
%     hold on
%     full_hist2 = sum(BIG_DATAFILE.electron{1,6,7,2}.ANGLE_ENERGY_HIST,1);
%     histogram('BinEdges',BINS,'BinCounts',full_hist2(1:end-1),'displayStyle','stairs','LineWidth',2)
%     full_hist2 = sum(BIG_DATAFILE.positron{1,6,7,2}.ANGLE_ENERGY_HIST,1);
%     histogram('BinEdges',BINS,'BinCounts',full_hist2(1:end-1),'displayStyle','stairs','LineWidth',2)
%     
%     legend('photon','electron','positron')
% end





%%

function ener_grid = load_energy_grid(data)
ener_grid= double(data.ENERGY_GRID);
end

%%

function [out_hist,bins] = get_angle_energy_histogram(data)
bins = load_energy_grid(data);
out_hist = double(data.ANGLE_ENERGY_HIST)./double(data.SAMPLED_NB);
end

%%

function [nb] = get_number_recorded(data)
nb = double(data.NB_RECORDED)./double(data.SAMPLED_NB);
end

%%

function [BIG_DATAFILE,ok]=add_to_BIG_DATAFILE(filename,BIG_DATAFILE,i_pot,i_recPos,i_efield_c,i_efield_s,PARAMETER_SET)

ok = false;

record_names={'photon','electron','positron'};

sett=load_settings();

% try
data_file = parse_output_file(filename);
% catch ME
%     disp(['skipped file ' num2str(filename) '...'])
%     ok = false;
%     return
% end

for i_r = 1:length(record_names)
    
    [out_hist,~] = get_angle_energy_histogram(data_file.(record_names{i_r}));
    
    nb_rec = get_number_recorded(data_file.(record_names{i_r}));
    
    if ~isfield(BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s},'NB_FILES')
        BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES = 0;
    end
    
    if ~isfield(BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s},'ANGLE_ENERGY_HIST')
        BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.ANGLE_ENERGY_HIST = zeros(13,256);
    end
    
    if ~isfield(BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s},'TOTAL_NB_RECORDED')
        BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.TOTAL_NB_RECORDED = 0;
    end
    
    if ~isfield(BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s},'ENERGY_GRID')
        BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID = data_file.(record_names{i_r}).ENERGY_GRID;
    end
    
    if ~isfield(BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s},'ZENITH_ANGLE_GRID')
        BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.ZENITH_ANGLE_GRID = data_file.(record_names{i_r}).ZENITH_ANGLE_GRID;
    end
    
    if ~isfield(BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s},'PARAMETER_SET')
        BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.PARAMETER_SET = PARAMETER_SET;
    end
    
    if ~isfield(BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s},'PARAMETER_SET_DESCRIPTION')
        BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.PARAMETER_SET_DESCRIPTION = 'POTENTIAL SELECTED_RECORD_POSITION EFIELD_SIZE_ALONG_ALT EFIELD_CENTER_ALT';
    end
    
%     out_hist(out_hist>10^20)=sett.MIN_HIST_BCKGRND;
%     out_hist(out_hist<10^-20)=sett.MIN_HIST_BCKGRND;
    
    BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.ANGLE_ENERGY_HIST = ...
        BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.ANGLE_ENERGY_HIST + out_hist;
    
    BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.TOTAL_NB_RECORDED = ...
        BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.TOTAL_NB_RECORDED + nb_rec;
    
    BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES = ...
        BIG_DATAFILE.(record_names{i_r}){i_pot,i_recPos,i_efield_c,i_efield_s}.NB_FILES + 1;
    
end

ok = true;

end

%%

function dir_list = get_dir_list(the_folder)

dir_list = dir(the_folder);
dir_list = dir_list([dir_list.isdir]);
dir_list = {dir_list([dir_list.isdir]).name};
dir_list = dir_list(~ismember(dir_list ,{'.','..',...
    'initial_photon','initial_electron','initial_positron',...
    'initial_muonN','initial_muonP','initial_neutron','initial_proton'}));

end