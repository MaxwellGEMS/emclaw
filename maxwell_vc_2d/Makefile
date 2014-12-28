all:
	make maxwell_2d_nc_rp.so
	make maxwell_2d_nc_tfluct.so
	make maxwell_2d_rp.so
	make maxwell_2d_tfluct.so

maxwell_2d_nc_rp.so: maxwell_2d_nc_fwave.f90
	f2py -m $(basename $(notdir $@)) -c $^

maxwell_2d_nc_tfluct.so: maxwell_2d_nc_tfluct.f90
	f2py -m $(basename $(notdir $@)) -c $^

maxwell_2d_rp.so: maxwell_2d_fwave.f90
	f2py -m $(basename $(notdir $@)) -c $^

maxwell_2d_tfluct.so: maxwell_2d_tfluct.f90
	f2py -m $(basename $(notdir $@)) -c $^