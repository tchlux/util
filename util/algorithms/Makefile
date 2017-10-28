all:
	-make -C ./fmodpy_VTdelaunay && make -C ./fmodpy_VTdelaunay clean
	-make -C ./fmodpy_linear_shepard && make -C ./fmodpy_linear_shepard clean
	-make -C ./fmodpy_marspack && make -C ./fmodpy_marspack clean
	-make -C ./fmodpy_max_box_mod && make -C ./fmodpy_max_box_mod clean
	-make -C ./fmodpy_fit_box_mod && make -C ./fmodpy_fit_box_mod clean

clean:
	-rm oldCode.txt
	-rm */oldCode.txt
	-rm temp-plot.html
	-rm -r __pycache__
	-rm -r */temp-plot.html
	-rm -r */__pycache__
