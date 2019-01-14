* At the end of each simulation run, 3 output files are filled in the folder `./output/...`
    * `photon_***` for photons, `electron_***` for electrons, `positron_***` for positrons.
* For each file:
    * First line has 10 values :
        * random number seed
        * number of events  = number of events started (i.e. initial cosmic particles)
        * initialy sampled PDG number type
        * record PDG type
        * center position of the electric field in km
        * full length of the electric field (along altitude) in km
        * applied potential in MV over the electric field region
        * relative weight of the initial sampled PDG number type
        * record altitude
        * cosmic ray generation altitude
    * next line is empty
    * next line has the number of records for the given particle type
    * next line is empty
    * next line is the zenith angle bins (degree) (= angle grid)
    * next line is empty
    * next line is the energy bins (MeV) (= energy grid)
    * next non-empty lines are the counts per energy bin for the several zenith-angular bins, last value is the overflow bin (>100 MeV)
