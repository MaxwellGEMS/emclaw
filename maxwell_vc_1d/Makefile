all:
	make maxwell_1d_rp.so
	make maxwell_1d_tfluct.so

maxwell_1d_rp.so: maxwell_1d_nl_fwave.f90
	f2py -m $(basename $(notdir $@)) -c $^

maxwell_1d_tfluct.so: maxwell_1d_nl_tfluct.f90
	f2py -m $(basename $(notdir $@)) -c $^