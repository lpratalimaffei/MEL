# MEL
Master Equation based Lumping code for integrating single PESs into global kinetic schemes

Requirements:

- Windows operating system (only tested on Windows for now)
- OpenSMOKE++ (available at https://www.opensmokepp.polimi.it/menu-download)
- python version > 3.6 ; version 3.9 might have problems with sklearn package

Installation: 

from MEL repository, run
python setup.py install

Testing:

All examples in the examples directory were already run. 
For each case, to run from scratch:
copy just the folder "inp" and the file input_lumping.txt

from the selected folder, run
MEL input_lumping.txt

Cite this work:

Pratali Maffei, L., Pelucchi, M., Cavallotti, C., Bertolino, A. and Faravelli, T., "Master equation lumping for multi-well potential energy surfaces: a bridge between ab initio based rate constants calculations and large kinetic mechanisms", Chemical Engineering Journal, 2021, just accepted


Contacts:
luna.pratali@polimi.it
