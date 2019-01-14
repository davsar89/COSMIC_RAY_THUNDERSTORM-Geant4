#pdg_phot = 22
#pdg_elec = 11
#pdg_posi = -11
#pdg_muN = 13
#pdg_muP = -13
#pdg_neut = 2112
#pdg_prot = 2212
# potential nb_to_shoot type_to_shoot sample_altitude record_altitude
timeout 1000s ./mos_test 0 20000 22 18 16 16 1&
timeout 1000s ./mos_test 0 20000 11 18 16 16 1&
timeout 1000s ./mos_test 0 20000 -11 18 16 16 1&
timeout 1000s ./mos_test 0 20000 13 18 16 16 1&
timeout 1000s ./mos_test 0 20000 -13 18 16 16 1&
timeout 1000s ./mos_test 0 20000 2112 18 16 16 1&
timeout 1000s ./mos_test 0 20000 2212 18 16 16 1&
