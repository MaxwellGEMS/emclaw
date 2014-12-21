all:
	make maxwell_1d_nc_rp.so
	make maxwell_1d_nc_tfluct.so
	make maxwell_1d_rp.so
	make maxwell_1d_tfluct.so

maxwell_1d_nc_rp.so: maxwell_1d_nc_fwave.f90
	f2py -m $(basename $(notdir $@)) -c $^

maxwell_1d_nc_tfluct.so: maxwell_1d_nc_tfluct.f90
	f2py -m $(basename $(notdir $@)) -c $^

maxwell_1d_rp.so: maxwell_1d_fwave.f90
	f2py -m $(basename $(notdir $@)) -c $^

maxwell_1d_tfluct.so: maxwell_1d_tfluct.f90
	f2py -m $(basename $(notdir $@)) -c $^