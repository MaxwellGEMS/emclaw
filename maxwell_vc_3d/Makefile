all:
	make maxwell_3d_rp.so
	make maxwell_3d_tfluct.so

maxwell_3d_rp.so: maxwell_3d_nl_fwave.f90
	f2py -m $(basename $(notdir $@)) -c $^

maxwell_3d_tfluct.so: maxwell_3d_nl_tfluct.f90
	f2py -m $(basename $(notdir $@)) -c $^