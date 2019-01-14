Using Docker, deploy the GEANT4 code rapidly on any computer (Linux, Mac Windows) and run simulations (but not for GUI and easy developement).

Requires `docker` and `docker-compose` installed.

Commands to run :
* `docker-compose up`

Output files produced by GEANT4 are written in the `./output` folder.

Will grab an Ubuntu 22.04 image with the full install Geant4 10.7.3 and dependencies, data. Will then compile and run the Geant4 project.

Modify the bash command and the end of the file `docker-compose.yml` to adapt it to another Geant4 code.
