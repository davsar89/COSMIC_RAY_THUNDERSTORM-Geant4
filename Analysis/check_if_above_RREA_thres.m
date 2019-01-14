function true_or_false = check_if_above_RREA_thres(pot,alt)
sett = load_settings();

RREA_thres = 2.84e5; % V/m
RREA_thres = RREA_thres / 1e6 * 1.e3; % MV/km
RREA_thres_MV = RREA_thres * sett.EFIELD_SIZE_list(2);
profile = get_density_profile();
RREA_thres_MV = RREA_thres_MV.*profile.vals;

true_or_false = abs(pot)>abs(interp1(profile.altitude,RREA_thres_MV,alt));
end