.PHONY: build clean

CC=gcc
CXX=g++

# Generate Python interface for module from C++ code using swig interface file (*.i)
distlink.py: distlink.i
	swig -python -c++ distlink.i

# Build extension
build: distlink.py setup.py
	CC=${CC} CXX=${CXX} python setup.py build_ext --inplace

# Clean temporary build files
clean:
	python setup.py clean --all

# Compile C++ test
test.out: test.cpp
	${CXX} test.cpp ./distlink-2.0/distlink.cpp -o test.out