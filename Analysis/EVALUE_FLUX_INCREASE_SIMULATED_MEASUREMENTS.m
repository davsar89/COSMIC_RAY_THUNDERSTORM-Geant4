clear all
close all
clc

sett = load_settings();

glow_db_path = [sett.base_path 'glow_database.mat'];
responses_fodler = '/Home/siv29/dsa030/Desktop/article_cr_modeling/MATLAB_ANALYSIS/responses/';

mission = 'ALOFT1';

% load glow simualtion database
loaded = load(glow_db_path);
glow_database = loaded.glow_database;

% load instrument response
loaded = load([responses_fodler mission]);
response = loaded.response;

increase=[];

s_egrid=[];


%%
for i_recPos = 1 : length(sett.RECORD_POS_LIST)
    for i_efield_c = 1 : length(sett.EFIELD_CENTER_list)
        for i_efield_s = 2
            for i_pot = 1:length(sett.POTENTIAL_LIST)
                
                POTENTIAL = sett.POTENTIAL_LIST(i_pot);
                ALT = sett.EFIELD_CENTER_list(i_efield_c);
                RECORD_POS = sett.RECORD_POS_LIST(i_recPos);
                EFIELD_SIZE = sett.EFIELD_SIZE_list(i_efield_s);
                
                rec_alt = ALT+RECORD_POS*EFIELD_SIZE/2.0;
                
                if rec_alt<=0 || rec_alt>=20
                    continue;
                end
                
                if isempty(s_egrid)
                    s_egrid = glow_database{1,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID;
                    s_egrid = (s_egrid(1:end-1)+s_egrid(2:end))/2.0;
                end
                
                photon_spec=glow_database{1,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST;
                electron_spec=glow_database{2,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST;
                positron_spec=glow_database{3,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST;
                
                r_egrid = response.energy_grid;
                photon_resp=response.effective_area.photon;
                electron_resp=response.effective_area.electron;
                positron_resp=response.effective_area.positron;
                
                megagrid = sort(unique([s_egrid r_egrid]));
                
                %% reference(background)
                [detection_factor_bg]=get_background(glow_database,...
                    i_pot,i_recPos,i_efield_c,i_efield_s,megagrid,s_egrid,r_egrid,photon_resp,electron_resp,positron_resp);
                
                
                %% enhanced (glow)
                
                photon_spec = interp1(s_egrid,photon_spec,megagrid);
                photon_resp = interp1(r_egrid,photon_resp,megagrid);
                photon_spec(isnan(photon_spec))=0;
                photon_resp(isnan(photon_resp))=0;
                detection_factor(1) = sum(photon_spec.*photon_resp);
                
                electron_spec = interp1(s_egrid,electron_spec,megagrid);
                electron_resp = interp1(r_egrid,electron_resp,megagrid);
                electron_spec(isnan(electron_spec))=0;
                electron_resp(isnan(electron_resp))=0;
                detection_factor(2) = sum(electron_spec.*electron_resp);
                
                positron_spec = interp1(s_egrid,positron_spec,megagrid);
                positron_resp = interp1(r_egrid,positron_resp,megagrid);
                positron_spec(isnan(positron_spec))=0;
                positron_resp(isnan(positron_resp))=0;
                detection_factor(3) = sum(positron_spec.*positron_resp);
                
                increase(end+1) = sum(detection_factor)/sum(detection_factor_bg);
                
                increase_due_to_photons = detection_factor(1)/detection_factor_bg(1);
                increase_due_to_electrons = detection_factor(2)/detection_factor_bg(2);
                increase_due_to_positrons = detection_factor(3)/detection_factor_bg(3);
                
                constraints = rec_alt>=18 && ALT<13 && increase(end)>1 && increase(end)<1.5;
                
                if constraints
                    disp(['POTENTIAL: ' num2str(POTENTIAL)])
                    disp(['EFIELD magnitude: ' num2str(POTENTIAL/EFIELD_SIZE) ' MV/km'])
                    disp(['REC POS: ' num2str(RECORD_POS)])
                    disp(['EFIELD CENTER: ' num2str(ALT)])
                    disp(['EFIELD FULL SIZE: ' num2str(EFIELD_SIZE)])
                    disp(['EFIELD EXTEND: ' num2str([ALT-EFIELD_SIZE/2]) '<->' num2str([ALT+EFIELD_SIZE/2])])
                    disp(['Rec alt: ' num2str(rec_alt)])
                    disp(['increase in record: ' num2str(increase(end))])
                    disp(['increase in records due to photons: ' num2str(increase_due_to_photons)])
                    disp(['increase in records due to electrons: ' num2str(increase_due_to_electrons)])
                    disp(['increase in records due to positrons: ' num2str(increase_due_to_positrons)])
                    disp(' ')
                end
                
            end
        end
    end
end

increase(increase<0)=0;

plot(increase,'+')


%%



function [detection_factor_bg]=get_background(glow_database,...
    i_pot,i_recPos,i_efield_c,i_efield_s,megagrid,s_egrid,r_egrid,photon_resp,electron_resp,positron_resp)

photon_spec_bg=glow_database{1,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV;
electron_spec_bg=glow_database{2,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV;
positron_spec_bg=glow_database{3,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV;

photon_spec_bg = interp1(s_egrid,photon_spec_bg,megagrid);
photon_resp = interp1(r_egrid,photon_resp,megagrid);
photon_spec_bg(isnan(photon_spec_bg))=0;
photon_resp(isnan(photon_resp))=0;
detection_factor_bg(1) = sum(photon_spec_bg.*photon_resp);

electron_spec_bg = interp1(s_egrid,electron_spec_bg,megagrid);
electron_resp = interp1(r_egrid,electron_resp,megagrid);
electron_spec_bg(isnan(electron_spec_bg))=0;
electron_resp(isnan(electron_resp))=0;
detection_factor_bg(2) = sum(electron_spec_bg.*electron_resp);

positron_spec_bg = interp1(s_egrid,positron_spec_bg,megagrid);
positron_resp = interp1(r_egrid,positron_resp,megagrid);
positron_spec_bg(isnan(positron_spec_bg))=0;
positron_resp(isnan(positron_resp))=0;
detection_factor_bg(3) = sum(positron_spec_bg.*positron_resp);

end