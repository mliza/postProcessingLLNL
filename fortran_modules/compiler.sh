module load python/3.8.2;
f2py3.8 -m f_mapping -c mapping_m.f90;
f2py3.8 -m f_gridReader -c gridReader_m.f90;
f2py3.8 -m f_scalarReader -c scalarReader_m.f90;
f2py3.8 -m f_vectorReader -c vectorReader_m.f90;
