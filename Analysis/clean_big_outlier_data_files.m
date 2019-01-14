function clean_big_outlier_data_files()

sett=load_settings();

FACTOR = 1000;
NB_FILES_THRESHOLD = 300;

IS_FRAM = sett.IS_FRAM;
base_path = sett.base_path;

sett=load_settings();
POTENTIAL_LIST = sett.POTENTIAL_LIST;
RECORD_POS_LIST = sett.RECORD_POS_LIST;
EFIELD_CENTER_list = sett.EFIELD_CENTER_list;
EFIELD_SIZE_list = sett.EFIELD_SIZE_list;
record_names = sett.record_names;

nb_recorded = cell(3,length(POTENTIAL_LIST),length(RECORD_POS_LIST),...
    length(EFIELD_CENTER_list),length(EFIELD_SIZE_list));

for i_pot =1:length(POTENTIAL_LIST)
    for i_recPos = 1:length(RECORD_POS_LIST)
        for i_efield_c = 1:length(EFIELD_CENTER_list)
            for i_efield_s = 1:length(EFIELD_SIZE_list)
                
                POTENTIAL = POTENTIAL_LIST(i_pot);
                SELECTED_RECORD_POSITION = RECORD_POS_LIST(i_recPos);
                EFIELD_CENTER = EFIELD_CENTER_list(i_efield_c);
                EFIELD_SIZE = EFIELD_SIZE_list(i_efield_s);
                
                rec_alt = EFIELD_CENTER + SELECTED_RECORD_POSITION*EFIELD_SIZE/2.0;
                
                if ~ismember(rec_alt, sett.WANTED_RECORD_ALTS)
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
                
                linear_indexing = 1:length(data_files_for_run);
                shuffled_indexing = linear_indexing(randperm(length(linear_indexing)));
                
                for i_file = shuffled_indexing
                    for i_r=1:length(record_names)
                        filename = [data_files_for_run(i_file).folder '/' data_files_for_run(i_file).name];
                        data_file = parse_output_file(filename);
                        this_nb_recorded = get_number_recorded(data_file.(record_names{i_r}));
                        
                        if ~isempty(nb_recorded{i_r,i_pot,i_recPos,i_efield_c,i_efield_s})
                            
                            mean_prev_nb_recorded = mean(nb_recorded{i_r,i_pot,i_recPos,i_efield_c,i_efield_s});
                            
                            nb_done_before = length(nb_recorded{i_r,i_pot,i_recPos,i_efield_c,i_efield_s});
                            
                            if (this_nb_recorded>FACTOR*mean_prev_nb_recorded && mean_prev_nb_recorded~=0 && nb_done_before>NB_FILES_THRESHOLD)
                                disp(num2str([this_nb_recorded mean_prev_nb_recorded this_nb_recorded/mean_prev_nb_recorded]))
                                pause(1)
                                disp(['deleting : ' filename ])
                                disp(num2str([mean_prev_nb_recorded this_nb_recorded nb_done_before]));
                                delete(filename)
                                %                                 disp('problem')
                                %                                 error('stop')
                                had_error = true;
                                clean_big_outlier_data_files();
                                return;
                            end
                        end
                        
                        nb_recorded{i_r,i_pot,i_recPos,i_efield_c,i_efield_s}(end+1) = this_nb_recorded;
                    end
                end
                
            end
        end
    end
end

end

%%

function [nb] = get_number_recorded(data)
nb = double(data.NB_RECORDED)./double(data.SAMPLED_NB).*1e3;
end