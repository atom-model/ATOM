# FIXME: just for Ian's Mac, which has a mystery /usr/local/bin/clang which is broken
export CC=/usr/bin/clang

all: pyagcm.so

#atmosphere:
	#${CC} 

pyagcm.so:
	python setup.py build_ext --inplace

clean:
	rm -f *.so *.o pyatom.cpp
	rm -rf build/
