all: pyagcm.so

pyagcm.so:
	python setup.py build_ext --inplace

clean:
	rm -f *.so *.o pyatom.cpp
	rm -rf build/
