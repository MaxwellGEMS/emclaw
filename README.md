emclaw
======

Solvers for Maxwell's electrodynamic equations with time-space varying coefficients based on
[PyClaw]([https://alals](http://www.clawpack.org/pyclaw/index.html))

## Setup

**Tested with Ubuntu 22.10, Python 3.9.**


### Clone this repository to your home directory

    git clone https://github.com/shredEngineer/emclaw.git /home/$USER/emclaw

### Create Anaconda environment

You need to have [Anaconda](https://docs.anaconda.com/anaconda/install/) installed. 

    conda create -n py39 python=3.9

### Install dependencies


    sudo apt install gfortran libblas-dev liblapack-dev
	conda install -c conda-forge numpy scipy matplotlib petsc4py clawpack

### Build and setup the `emclaw` package

	export EMCLAW_PATH=/home/$USER/emclaw
	cd $EMCLAW_PATH/emclaw/riemann

	conda activate py39
	python setup.py build_ext -i
	conda develop -n py39 $EMCLAW_PATH

### Quick Test 

	python -c "import emclaw; print(emclaw.__path__[0])"

## Examples

Run the Maxwell 2D example:

	cd /home/$USER/emclaw/examples/maxwell_vc_2d
	python maxwell_2d.py

The console output should look like this:

	averaging
	2022-12-30 10:27:00,924 INFO CLAW: Solution 0 computed for time t=0.000000
	2022-12-30 10:27:01,734 INFO CLAW: Solution 1 computed for time t=10.000000
	2022-12-30 10:27:02,526 INFO CLAW: Solution 2 computed for time t=20.000000
	2022-12-30 10:27:03,304 INFO CLAW: Solution 3 computed for time t=30.000000
	2022-12-30 10:27:04,080 INFO CLAW: Solution 4 computed for time t=40.000000
	2022-12-30 10:27:04,855 INFO CLAW: Solution 5 computed for time t=50.000000
	2022-12-30 10:27:05,636 INFO CLAW: Solution 6 computed for time t=60.000000
	2022-12-30 10:27:06,413 INFO CLAW: Solution 7 computed for time t=70.000000
	2022-12-30 10:27:07,193 INFO CLAW: Solution 8 computed for time t=80.000000
	2022-12-30 10:27:07,974 INFO CLAW: Solution 9 computed for time t=90.000000
	2022-12-30 10:27:08,750 INFO CLAW: Solution 10 computed for time t=100.000000

This should have generated a new subfolder named `_output` where the results are stored.

**TODO:** Generate plots.

## Support
[damian.sanroman@kaust.edu.sa](mailto:damian.sanroman@kaust.edu.sa)
