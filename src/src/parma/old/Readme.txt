PARMA4 can calculate the atmospheric cosmic-ray fluxes & dose rates.
Details of the model is described in
T.Sato, Analytical model for estimating terrestrial cosmic-ray fluxes nearly anytime and anywhere in the world; Extension of PARMA/EXPACS, 10(12): e0144679 (2015)
T.Sato, Analytical model for estimating the zenith angle dependence of terrestrial cosmic ray fluxes, PLOS ONE, 11(8): e0160390 (2016)
This program is designed to be implemented in some other systems such as route-dose calculation systems. 
If you would like to calculate cosmic-ray fluxes or doses at certain conditions, 
it is much easier to use EXPACS.
You can download EXPACS from the webiste below
http://phits.jaea.go.jp/expacs

Table of contents
  Readme.txt:      This file
  main.f90:        Main program of PARMA for calculating cosmic-ray fluxes & doses
  main-simple.f90: Main program for calculating flux at a certain condition (easy to start up)
  main-generation.f90: Main program for generating energy and direction of cosmic-rays at a certain condition (for Monte Carlo simulation)
  subroutines.f90: Subroutines of PARMA
  CheckGeneration.out: Sample output file from main-generator.f90
  /AngOut:         Sample output files for angular distribution calculation
  /condition:      Sample input files for flux & dose calculation
  /dcc:            Database for fluence to dose (or other quantity) conversion coefficients
  /Doseout:        Sample output files for dose calculation
  /input:          Databases used in PARMA
  /SpecOut:        Sample output files for flux calculation

Q1. How to use?
A1. You have to compile the programs using a Fortran compiler.
    Intel Fortran or Gfortran are recommended to be used.

Q2. How to compile and execute?
A2. To start up, it is easier to use main-simple.f90 as the control routine.
    In that case, you have to compile main-simple.f90 and subroutines.f90 at once.
    The followings are an example of commands that should be typed (Gfortran case)
    > gfortran main-simple.f90 subroutines.f90
    > a.exe

Q3. How to change the condition?
A3. You can change the condition such as time and location by changing
    the parameters written in main-simple.f90.

Q4. How to use this program as an event generator?
A4. You have to use main-generator.f90 instead of main-simple.f90 as the control routine.
    You have to specify the energy and angle ranges for generation, as well as the
    number of particles to be generated. 
    The variables e,u,v,w are the energy & direction vector of the generated particle.
    Total flux (/cm2/s) is output in "CheckGeneration.out" file.

Q5. How to calculate dose?
A5. You have to use main.f90 instead of main-simple.f90 as the control routine.
    To specify the calculation condition, you have to make your own input file
    in "condition" folder by referring sample input files such as "Tokyo-Smin.inp".
    At the 1st line of the input file, you have to specify isout and istyle parameters.
    isout=0: No output for flux data for each condition, =1: Output Flux data in "SpecOut" folder
    istyle= 0:Direct input (s****-r***-d****-g***)
          = 1:year.month, latitude, longitude, atmospheric detph(g/cm^2), g(direct input)
          = 2:year.month, latitude, longitude, altitude (m), g(direct input) 
          = 3:year.month, latitude, longitude, altitude (ft), g(direct input)
          = 4:year.month, cutoff rigidity, altitude (m), g(direct input)
    After 2nd lines, you have to specify the conditions

Q6. How to change the surrounding conditions, such as ground or aircraft
A6. You have to change "g" parameter.
    If you want to calculate the cosmic-ray spectra in the ideal atmosphere
    (i.e. without considering the surrounding effect), you should set g=10.0
    If you want to calculate the ground level cosmic-ray spectra, you should
    set 0 =< g =< 1. In this case, "g" means the water fraction in ground.
    Ground level muon correction is also considered in this mode.
    If you do not know the water fraction at your specified location, 
    I recommended to use 0.15 for "g".
    If you want to calculate the cosmic-ray spectra in aircraft, you should
    set -10 < g < 0. In this case, the absolute value of "g" indicates
    the mass of the aircraft in the unit of 100 ton. 
    It should be noted that only neutron spectrum is influenced by 
    the surrounding condition.
    If you want to calculate the angular differential cosmic-ray fluxes without
    the albedo from the Earth (i.e. black hole mode), you should set g=100.0.

Q7. How to calculate other quantity such as count rates of neutron monitor?
A7. You can select the evaluation quantity by changing "dccname" parameter in "main.f90".
    Four quantities can be calculated, which are the effective dose for isotropic
    irradiation, H*(10), dose rate in air, and count rate of 6 tube neutron monitor.

Q8. What is the conditions to use this program?
A8. For non-commercial use, you should refer the following manuscripts in any published use of this program,
    T.Sato, Analytical model for estimating terrestrial cosmic-ray fluxes nearly anytime and anywhere in the world; Extension of PARMA/EXPACS, 10(12): e0144679 (2015)
    T.Sato, Analytical model for estimating the zenith angle dependence of terrestrial cosmic ray fluxes, PLOS ONE, 11(8): e0160390 (2016)
    The commercial use of this program is NOT allowed with a prior agreement with JAEA

Acknowledgement
  This program uses a random generator based on Mersenne Twister method,
  which is developed by Prof. M. Matsumoto of Hiroshima University.

Contact
  If you have any questions or requests on this program, 
  please E-mail to nsed-expacs@jaea.go.jp
