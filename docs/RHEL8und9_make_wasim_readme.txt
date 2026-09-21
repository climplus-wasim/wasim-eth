How to build WaSiM in a VM running RHEL 9: 

mount the windows directory called E containing the sources (if not stored locally):
E must be configured as a shared directory (e.g. in virtualBox)

in a terminal:
su (enter passwort) (or prefix all commands with sudo)
umount E
mount -t vboxsf E /media/sf_E
cd /media/sf_E/<windows source code directory relative to shared directory E>


use of g++ version 11 and higher
using makefile_rhel89_openmpi
optional, depending on system/g++ version: before running the make command, the LD_LIRARY_PATH must be set in order to to use the latest versions of libstdc++ and openMPI.

then: load modul openmpi 
then: run the make command:
make -f makefile_rhel9_openmpi

after success, try out the compiled model:

LD_LIBRARY_PATH="$PATH:/usr/lib64/openmpi/lib"
module load openmpi-x86_64
export OMP_NUM_THREADS=1
cd /media/sf_E/<local path to sample wasim model setup>
mpiexec -n 18 --allow-run-as-root /media/sf_E/<local path to compiled model>/WaSiM-11.01.01/wasimuzr-11-01-01-x86-64-pc-linux-rhel9-openmpi sample_controlfile.ctl

using VSCode:
1. Install VSCode 
2. Install C++ extension
3. in makefile: insert -g -O0 
4. remove -fopenmp if one single thread is ok
5. copy latest sources to /home/<user>/wasim/xx.xx.xx
6. in lanuch.json: set name of the program to start with the correct version xx.xx.xx and set controlfile as parameter
7. open a terminal
   (when running in a VBox-environment and the directory to be mapped has the name E)
	sudo umount E
	sudo mount -tvboxsf E /medis/sf_E/
	module load mpi/openmpi-x86_64
	code &
the last command starts VS Code. There, run 
	tasks->clean
then
	tasks->build
then start debug (but make sure to set some brakpoints before)


RHEL93:
sample commands to run the model:

LD_LIBRARY_PATH="$PATH:/usr/lib64/openmpi/lib"
module load mpi/openmpi-x86_64
export OMP_NUM_THREADS=1
cd /media/sf_E/BirsErgolz/Calibration/Control
mpiexec -n 18 --allow-run-as-root /media/sf_E/m500/source/WaSiM-11.01.01/wasimuzr-11-01-01-x86-64-pc-linux-rhel9-openmpi BirsErg_500_0019_1991.ctl

