clear all
close all
clc

% ! sh COPY_FINAL_DATAFILES_FROM_FRAM.sh

loaded = load('BIG_DATAFILE_all.mat');
BIG_DATAFILE = loaded.BIG_DATAFILE;

sett=load_settings();
POTENTIAL_LIST = sett.POTENTIAL_LIST;
RECORD_POS_LIST = sett.RECORD_POS_LIST;
EFIELD_CENTER_list = sett.EFIELD_CENTER_list;
EFIELD_SIZE_list = sett.EFIELD_SIZE_list;

POT = 50;
ALT = 16;
RECORD_POS = 0;
EFIELD_SIZE = 1;

i_pot = find_idx(POTENTIAL_LIST,POT);
i_recPos = find_idx(RECORD_POS_LIST,RECORD_POS);
i_efield_c = find_idx(EFIELD_CENTER_list,ALT);
i_efield_s = find_idx(EFIELD_SIZE_list,EFIELD_SIZE);

kept.photon = BIG_DATAFILE.photon{i_pot,i_recPos,i_efield_c,i_efield_s};
kept.electron = BIG_DATAFILE.electron{i_pot,i_recPos,i_efield_c,i_efield_s};
kept.positron = BIG_DATAFILE.positron{i_pot,i_recPos,i_efield_c,i_efield_s};

ZENITH_ANGLE_GRID = kept.photon.ZENITH_ANGLE_GRID;

ENERGY_GRID = kept.photon.ENERGY_GRID;

NB_ANGLES = length(kept.photon.ZENITH_ANGLE_GRID);

weight_type = [kept.photon.TOTAL_NB_RECORDED kept.electron.TOTAL_NB_RECORDED kept.positron.TOTAL_NB_RECORDED] ./ sum([kept.photon.TOTAL_NB_RECORDED kept.electron.TOTAL_NB_RECORDED kept.positron.TOTAL_NB_RECORDED])*100.0;

weight_angles(1,:) = sum(kept.photon.ANGLE_ENERGY_HIST,2)./ sum(sum(kept.photon.ANGLE_ENERGY_HIST,2))*100.0;
weight_angles(2,:) = sum(kept.electron.ANGLE_ENERGY_HIST,2)./ sum(sum(kept.electron.ANGLE_ENERGY_HIST,2))*100.0;
weight_angles(3,:) = sum(kept.positron.ANGLE_ENERGY_HIST,2)./ sum(sum(kept.positron.ANGLE_ENERGY_HIST,2))*100.0;

weight_energy_for_angles={};

for i_type=1:3
    for i_ang=1:NB_ANGLES
        weight_energy_for_angles{i_type}(i_ang,:) = kept.photon.ANGLE_ENERGY_HIST(i_ang,:)./sum( kept.photon.ANGLE_ENERGY_HIST(i_ang,:))*100.0;
    end
    weight_energy_for_angles{i_type}(isnan(weight_energy_for_angles{i_type}))=0.0;
end

%%
output_string{1,1} = [num2str(ZENITH_ANGLE_GRID,' %e')];
output_string{end+1,1} = [num2str(ENERGY_GRID,' %e')];

output_string{end+1,1} = [num2str(weight_type,' %e')];

for i_type=1:3
    output_string{end+1,1} = [num2str(weight_angles(i_type,:),' %e')];
end

for i_type=1:3
    output_string{end+1,1} = [num2str(9.999999,' %e')];
    for i_ang=1:NB_ANGLES
        output_string{end+1,1} = [num2str(weight_energy_for_angles{i_type}(i_ang,:),' %e')];
    end
end

%%
out_filename = ['input_' num2str(POT) '_' num2str(ALT) '_' num2str(RECORD_POS) '_' num2str(EFIELD_SIZE) '.in'];

fileID = fopen(out_filename,'w');
for is=1:length(output_string)
    fprintf(fileID,'%s \n',output_string{is});
end

fclose(fileID);
