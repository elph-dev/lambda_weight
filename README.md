Weighted average method for calculating electron-phonon couplings
========
lambda_weight is a program to calculate electron-phonon couplings using weighted average method.
For details, see

  T. Koretsune and R. Arita  in preparation.

Currently, the program is based on the quantum-ESPRESSO v6.0.

# Install ######################################
	1. Download or clone the github repository, e.g.
		> git clone https://github.com/elph-dev/lambda_weight

	2. Apply a patch to the qe-6.0/PHonon/PH directory.
		> cd PATH/TO/QEDIR/qe-6.0/PHonon
		> patch -p1 < LAMBDA_WEIGHT_DIR/qe-6.0/lambda_weight.patch

	3. Compile
		> cd PHelph
		> make ( lambda_weight.x )

# Usage
	You can find an example in lambda_weight/qe-6.0/examples/Nb.

