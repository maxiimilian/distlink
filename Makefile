.PHONY: build clean

CC=gcc
CXX=g++

distlink.py: distlink.i
	swig -python -c++ distlink.i

build: distlink.py setup.py
	CC=${CC} CXX=${CXX} python setup.py build_ext --inplace

clean:
	python setup.py clean --all