Weighted average method for lambda
========
lambda_weight is a program to calculate electron-phonon couplings using weighted average method.
For details, see

  T. Koretsune and R. Arita, Computer Physics Communications 220 239-242 (2017) (https://doi.org/10.1016/j.cpc.2017.07.011) 
  ( [arxiv:1610.09441](http://arxiv.org/abs/1610.09441) )

Currently, the program is based on the quantum-ESPRESSO v6.0.

# Install ######################################
	1. Download or clone the github repository, e.g.
		> git clone https://github.com/elph-dev/lambda_weight

	2. Apply a patch to the qe-6.0/PHonon/PH directory.
		> cd PATH/TO/QEDIR/qe-6.0/PHonon
		> cp -a PH PHelph
		> cd PHelph
		> patch -p1 < LAMBDA_WEIGHT_DIR/qe-6.0/lambda_weight.patch

	3. Compile
		> cd PHelph
		> make ( lambda_weight.x )

# Usage
	You can find an example in qe-6.0/examples/Nb. See README in this directory for more details.

