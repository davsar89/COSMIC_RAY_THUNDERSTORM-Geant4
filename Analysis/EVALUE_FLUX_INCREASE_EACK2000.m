clear all
close all
clc

set(0, 'DefaultFigureRenderer', 'painters');
PLOT_GRAY_AREA=false;

global APPLY_RESP
global LOG_INTERP

APPLY_RESP = true;
LOG_INTERP = true;

sett = load_settings();

mission = 'EACK2000';
WANTED_REC_ALT = 14;
MAX_WANTED_INCREASE = 3;
MIN_WANTED_INCREASE = 1.6;

YLIM_PLOT = [10.0 16];
XLIM_PLOT = [-222 222];

SMOOTHING_LEVEL = sett.SMOOTHING_LEVEL;

%%
NB_SUPER_SAMPLING_PLOT = sett.NB_SUPER_SAMPLING_PLOT;

glow_db_path = [sett.base_path 'glow_database.mat'];
responses_fodler = sett.responses_fodler;

% load glow simualtion database
loaded = load(glow_db_path);
glow_database = loaded.glow_database;

% load instrument response
loaded = load([responses_fodler mission]);
response = loaded.response;

increase_list = [];

s_egrid_centers = [];
s_egrid = [];

max_rec_alt = 0;

increase_for_level_curve_plot=zeros(length(sett.EFIELD_CENTER_list),length(sett.POTENTIAL_LIST));

POTENTIAL_grid=zeros(length(sett.EFIELD_CENTER_list),length(sett.POTENTIAL_LIST));

%%
for i_recPos = 1 : length(sett.RECORD_POS_LIST)
    for i_efield_c = 1 : length(sett.EFIELD_CENTER_list)
        for i_efield_s = 2
            for i_pot = 1:length(sett.POTENTIAL_LIST)
                
                POTENTIAL = sett.POTENTIAL_LIST(i_pot);
                ALT = sett.EFIELD_CENTER_list(i_efield_c);
                RECORD_POS = sett.RECORD_POS_LIST(i_recPos);
                EFIELD_SIZE = sett.EFIELD_SIZE_list(i_efield_s);
                
                rec_alt = ALT + RECORD_POS*EFIELD_SIZE/2.0;
                
                if rec_alt~=WANTED_REC_ALT
                    continue;
                end
                
                if isempty(glow_database{1,i_pot,i_recPos,i_efield_c,i_efield_s})
                    continue
                end
                
                if isempty(s_egrid)
                    s_egrid = glow_database{1,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_GRID;
                    s_egrid_centers = (s_egrid(1:end-1)+s_egrid(2:end))/2.0;
                end
                
                r_egrid = response.energy_grid;
                
                photon_resp = response.effective_area.photon;
                electron_resp = response.effective_area.electron;
                positron_resp = response.effective_area.positron;
                
                all_egrid_vals = sort(unique([s_egrid_centers r_egrid]));
                
                % megagrid = linspace(min(all_egrid_vals),max(all_egrid_vals),4096*100);
                megagrid = logspace(log10(min(all_egrid_vals)),log10(max(all_egrid_vals)),4096/4);
                
                %% reference(background)
                [detection_factor_bg]=get_background(...
                    glow_database,...
                    i_recPos,i_efield_c,i_efield_s,...
                    megagrid,s_egrid_centers,r_egrid,...
                    photon_resp,electron_resp,positron_resp);
                
                %% enhanced (glow)
                [detection_factor]=get_glow(...
                    glow_database,...
                    i_pot,i_recPos,i_efield_c,i_efield_s,...
                    megagrid,s_egrid_centers,r_egrid,...
                    photon_resp,electron_resp,positron_resp);
                
                %%
                
                increase_due_to_photons = detection_factor(1)/detection_factor_bg(1);
                increase_due_to_electrons = detection_factor(2)/detection_factor_bg(2);
                increase_due_to_positrons = detection_factor(3)/detection_factor_bg(3);
                
                if APPLY_RESP
                    increase_list(end+1) = sum(detection_factor.*sett.particles_extra_weight)/sum(detection_factor_bg.*sett.particles_extra_weight);
                else
                    increase_list(end+1) = increase_due_to_photons;
                end
                
                constraints = rec_alt==WANTED_REC_ALT && increase_list(end)>=MIN_WANTED_INCREASE && increase_list(end)<=MAX_WANTED_INCREASE;
                
                if constraints
                    disp(['POTENTIAL: ' num2str(POTENTIAL)])
                    disp(['EFIELD magnitude: ' num2str(POTENTIAL/EFIELD_SIZE) ' MV/km'])
                    disp(['REC POS: ' num2str(RECORD_POS)])
                    disp(['EFIELD CENTER: ' num2str(ALT)])
                    disp(['EFIELD FULL SIZE: ' num2str(EFIELD_SIZE)])
                    disp(['EFIELD EXTEND: ' num2str([ALT-EFIELD_SIZE/2]) '<->' num2str([ALT+EFIELD_SIZE/2])])
                    disp(['Rec alt: ' num2str(rec_alt)])
                    disp(['increase in record: ' num2str(increase_list(end))])
                    disp(['increase in records due to photons: ' num2str(increase_due_to_photons)])
                    disp(['increase in records due to electrons: ' num2str(increase_due_to_electrons)])
                    disp(['increase in records due to positrons: ' num2str(increase_due_to_positrons)])
                    disp(' ')
                end
                
                if rec_alt==WANTED_REC_ALT
                    increase_for_level_curve_plot(i_efield_c,i_pot) = increase_list(end);
                    
                    POTENTIAL_grid(i_efield_c,i_pot) = POTENTIAL;
                end
                
            end
        end
    end
end

%%
increase_for_level_curve_plot(increase_for_level_curve_plot<=0)=min(increase_for_level_curve_plot(increase_for_level_curve_plot~=0));
increase_for_level_curve_plot(isnan(increase_for_level_curve_plot))=1.0;
increase_for_level_curve_plot(isinf(increase_for_level_curve_plot))=1.0;
increase_for_level_curve_plot(increase_for_level_curve_plot<1.0) = 1.0;
increase_list(increase_list<1.0) = 1.0;

%%

figure('WindowStyle','normal','Position', [0 0 1200 600].*0.9)
hold on

[X,Y] = meshgrid(sett.POTENTIAL_LIST,sett.EFIELD_CENTER_list);

pot_megagrid = linspace(min(sett.POTENTIAL_LIST),max(sett.POTENTIAL_LIST),NB_SUPER_SAMPLING_PLOT);
pot_megagrid = unique(sort([sett.POTENTIAL_LIST pot_megagrid]));
alt_megagrid = linspace(min(sett.EFIELD_CENTER_list),max(sett.EFIELD_CENTER_list),NB_SUPER_SAMPLING_PLOT);
alt_megagrid = unique(sort([sett.EFIELD_CENTER_list alt_megagrid]));

increase_for_level_curve_plot_interpolated = zeros(NB_SUPER_SAMPLING_PLOT,NB_SUPER_SAMPLING_PLOT+1);

if SMOOTHING_LEVEL==2
    for ii=1:size(increase_for_level_curve_plot,2)
        increase_for_level_curve_plot(:,ii)=smooth(increase_for_level_curve_plot(:,ii));
    end
end

if SMOOTHING_LEVEL==3
    for ii=1:size(increase_for_level_curve_plot,1)
        increase_for_level_curve_plot(ii,:)=smooth(increase_for_level_curve_plot(ii,:));
    end
end

if SMOOTHING_LEVEL==4
    for ii=1:size(increase_for_level_curve_plot,2)
        increase_for_level_curve_plot(:,ii)=smooth(increase_for_level_curve_plot(:,ii));
    end
    for ii=1:size(increase_for_level_curve_plot,1)
        increase_for_level_curve_plot(ii,:)=smooth(increase_for_level_curve_plot(ii,:));
    end
end

pot_list_positive=[];
alt_list_positive=[];
pot_list_negative=[];
alt_list_negative=[];

for ii=1:length(pot_megagrid)
    for jj=1:length(alt_megagrid)
        
        pot = pot_megagrid(ii);
        alt = alt_megagrid(jj);
        
        val = interp2(X,Y,log(increase_for_level_curve_plot),pot,alt,'linear');
        val = exp(val);
        
        if abs(pot)<=1
            val=1;
        end
        
        increase_for_level_curve_plot_interpolated(jj,ii) = val;
        
        if val>=MIN_WANTED_INCREASE && val<=MAX_WANTED_INCREASE
%         if val>=MIN_WANTED_INCREASE
            if pot>=0
                pot_list_positive(end+1)=pot;
                alt_list_positive(end+1)=alt;
            else
                pot_list_negative(end+1)=pot;
                alt_list_negative(end+1)=alt; 
            end
        end
    end
end

%%

scatter(pot_list_positive,alt_list_positive,'g','filled')
scatter(pot_list_negative,alt_list_negative,'g','filled')


% k = boundary(pot_list_positive(:),alt_list_positive(:));
% fill(pot_list_positive(k),alt_list_positive(k),'g','LineStyle','none')
% 
% k2 = boundary(pot_list_negative(:),alt_list_negative(:));
% fill(pot_list_negative(k2),alt_list_negative(k2),'g','LineStyle','none')


%%

min_val=min(increase_for_level_curve_plot(:));
max_val=max(increase_for_level_curve_plot(:));
grid_val = round(logspace(log10(min_val),log10(max_val),20).*10.0)/10.0;
grid_val = unique(sort([1.05 1.1 1.2 1.3 1.4 1.5 1.6 1.7 grid_val]));
[X2,Y2] = meshgrid(pot_megagrid,alt_megagrid);

[C,h] = contour(X2,Y2,increase_for_level_curve_plot_interpolated,grid_val);
% [C,h] = contour(X,Y,increase_for_level_curve_plot,linspace(log10(0.2),log10(10),20));

clabel(C,h);

set(gca,'ylim',YLIM_PLOT)
set(gca,'xlim',XLIM_PLOT)

if PLOT_GRAY_AREA
    ybars = [14.05 20];
    pp = patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1), ybars(2) ybars(2)], [0.8 0.8 0.8],'LineStyle','none');
    alpha(0.8)
end

plot([-220 220],[WANTED_REC_ALT WANTED_REC_ALT],'-r','linewidth',8)

%% plot RREA theshold as function of altitude

RREA_thres = 2.84e5; % V/m
RREA_thres = RREA_thres / 1e6 * 1.e3; % MV/km
RREA_thres_MV = RREA_thres * sett.EFIELD_SIZE_list(2);
profile = get_density_profile();
RREA_thres_MV = RREA_thres_MV.*profile.vals;

plot(RREA_thres_MV,profile.altitude,'b','linewidth',8)
plot(-RREA_thres_MV,profile.altitude,'b','linewidth',8)

% title(mission,'interpreter','latex','fontsize',sett.figure_font_size)
xlabel('Potential $\Delta U$ (MV)','interpreter','latex','fontsize',sett.figure_font_size)
ylabel('Electric field center altitude $H_E$ (km)','interpreter','latex','fontsize',sett.figure_font_size)

grid on
xticks(unique(sort([0 -220:20:220])));
h = breakxaxis([0 80]);

path_to_folder = '/Home/siv29/dsa030/Desktop/article_cr_modeling/figures/';
saveas(gcf,[path_to_folder mission '.eps'],'epsc')
saveas(gcf,[path_to_folder mission '.png'],'png')

path_to_folder = '/Home/siv29/dsa030/Desktop/article_cr_modeling/DATA_REPO/ARTICLE_DATA/Figure_5/';
saveas(gcf,[path_to_folder mission '.eps'],'epsc')
saveas(gcf,[path_to_folder mission '.png'],'png')
savefig(gcf,[path_to_folder mission '.fig'])

%%

plot_command = {...
    'min_val = min(increase_factor(:));',...
    'max_val = max(increase_factor(:));',...
    'grid_val = round(logspace(log10(min_val),log10(max_val),20).*10.0)/10.0;',...
    'grid_val = unique(sort([1.05 1.1 1.2 1.3 1.4 1.5 1.6 1.7 grid_val]));',...
    '[X2,Y2] = meshgrid(potential_grid,altitude_grid);',...
    '[C,h] = contour(X2,Y2,increase_for_level_curve_plot_interpolated,grid_val);'};

[C,h] = contour(X2,Y2,increase_for_level_curve_plot_interpolated,grid_val);

datafile=[];
datafile.name=mission;
datafile.potential_grid = pot_megagrid;
datafile.altitude_grid = alt_megagrid;
datafile.potential_unit='MV';
datafile.altitude_unit='km';
datafile.increase_factor = increase_for_level_curve_plot_interpolated;
datafile.plot_command = plot_command;
save([path_to_folder mission '.mat'],'datafile');




%%
function [detection_factor_bg]=get_background(glow_database,...
    i_recPos,i_efield_c,i_efield_s,...
    megagrid,s_egrid_cen,r_egrid,...
    photon_resp,electron_resp,positron_resp)
global APPLY_RESP;

supersum = sum([photon_resp(~isnan(photon_resp)) electron_resp(~isnan(electron_resp)) positron_resp(~isnan(positron_resp))]);
photon_resp = photon_resp./supersum;
electron_resp = electron_resp./supersum;
positron_resp = positron_resp./supersum;

sett = load_settings();

i_0MV = sett.i_0MV;

photon_spec_bg=glow_database{1,i_0MV,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV;
electron_spec_bg=glow_database{2,i_0MV,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV;
positron_spec_bg=glow_database{3,i_0MV,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST_0MV;

photon_spec_bg = interp1_log(s_egrid_cen,photon_spec_bg,megagrid);
photon_resp = interp1_log(r_egrid,photon_resp,megagrid);
photon_spec_bg(isnan(photon_spec_bg))=0;
photon_resp(isnan(photon_resp))=0;
if APPLY_RESP
        detection_factor_bg(1) = trapz(photon_spec_bg.*photon_resp,megagrid);
%     detection_factor_bg(1) = sum(photon_spec_bg.*photon_resp);
else
    detection_factor_bg(1) = trapz(photon_spec_bg,megagrid);
end

electron_spec_bg = interp1_log(s_egrid_cen,electron_spec_bg,megagrid);
electron_resp = interp1_log(r_egrid,electron_resp,megagrid);
electron_spec_bg(isnan(electron_spec_bg))=0;
electron_resp(isnan(electron_resp))=0;
if APPLY_RESP
        detection_factor_bg(2) = trapz(electron_spec_bg.*electron_resp,megagrid);
%     detection_factor_bg(2) = sum(electron_spec_bg.*electron_resp);
else
    detection_factor_bg(2) = trapz(electron_spec_bg,megagrid);
end

positron_spec_bg = interp1_log(s_egrid_cen,positron_spec_bg,megagrid);
positron_resp = interp1_log(r_egrid,positron_resp,megagrid);
positron_spec_bg(isnan(positron_spec_bg))=0;
positron_resp(isnan(positron_resp))=0;
if APPLY_RESP
        detection_factor_bg(3) = trapz(positron_spec_bg.*positron_resp,megagrid);
%     detection_factor_bg(3) = sum(positron_spec_bg.*positron_resp);
else
    detection_factor_bg(3) = trapz(positron_spec_bg,megagrid);
end

end

%%

function [detection_factor]=get_glow(glow_database,...
    i_pot,i_recPos,i_efield_c,i_efield_s,...
    megagrid,s_egrid_cen,r_egrid,...
    photon_resp,electron_resp,positron_resp)
global APPLY_RESP

supersum = sum([photon_resp(~isnan(photon_resp)) electron_resp(~isnan(electron_resp)) positron_resp(~isnan(positron_resp))]);
photon_resp = photon_resp./supersum;
electron_resp = electron_resp./supersum;
positron_resp = positron_resp./supersum;

photon_spec = glow_database{1,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST;
electron_spec = glow_database{2,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST;
positron_spec = glow_database{3,i_pot,i_recPos,i_efield_c,i_efield_s}.ENERGY_HIST;

photon_spec = interp1_log(s_egrid_cen,photon_spec,megagrid);
photon_resp = interp1_log(r_egrid,photon_resp,megagrid);
photon_spec(isnan(photon_spec)) = 0;
photon_resp(isnan(photon_resp)) = 0;
if APPLY_RESP
        detection_factor(1) = trapz(photon_spec.*photon_resp,megagrid);
%     detection_factor(1) = sum(photon_spec.*photon_resp);
else
    detection_factor(1) = trapz(photon_spec,megagrid);
end

electron_spec = interp1_log(s_egrid_cen,electron_spec,megagrid);
electron_resp = interp1_log(r_egrid,electron_resp,megagrid);
electron_spec(isnan(electron_spec)) = 0;
electron_resp(isnan(electron_resp)) = 0;
if APPLY_RESP
        detection_factor(2) = trapz(electron_spec.*electron_resp,megagrid);
%     detection_factor(2) = sum(electron_spec.*electron_resp);
else
    detection_factor(2) = trapz(electron_spec,megagrid);
end

positron_spec = interp1_log(s_egrid_cen,positron_spec,megagrid);
positron_resp = interp1_log(r_egrid,positron_resp,megagrid);
positron_spec(isnan(positron_spec)) = 0;
positron_resp(isnan(positron_resp)) = 0;
if APPLY_RESP
        detection_factor(3) = trapz(positron_spec.*positron_resp,megagrid);
%     detection_factor(3) = sum(positron_spec.*positron_resp);
else
    detection_factor(3) = trapz(positron_spec,megagrid);
end

end

%%

function val = interp1_log(x,y,xv)
global LOG_INTERP

if LOG_INTERP
    if sum(y)~=0
        y(y==0)=min(y(y~=0))/100.0;
    end
    
    val = interp1(log10(x),log10(y),log10(xv));
    val = 10.^(val);
else
    val = interp1(x,y,xv);
end

end

%%