Using Docker, deploy the GEANT4 code rapidly on any computer (Linux, Mac Windows) and run simulations (but not for GUI and easy developement).

Requires `docker` and `docker-compose` installed.

Commands to run :
* `docker-compose pull`
* `docker-compose build`
* `docker-compose up`

The folders `interface` and `output` are in common between this computer and the Docker Virtual Machine. `interface` is used for the installation process. Output files produced by GEANT4 are written in the `output` folder.

Will grab an Ubuntu 18.04 image, install updates, required packages, download compile and install Geant4 and dependencies. Will then compile and run the Geant4 project.

Modify the bash command and the end of the file `docker-compose.yml` to adapt it to another Geant4 code.
