CFLAGS_OPENCV=`pkg-config opencv --cflags --libs`
CFLAGS_EIGEN=`pkg-config eigen3 --cflags --libs`
CXX ?= g++
CC ?= gcc
CFLAGS = -Wl,--as-needed -fPIC
CXXFLAGS = -Wl,--as-needed -fPIC
MEX = mex
MEX_OPTION = CC="$(CC)" CFLAGS="$(CFLAGS)" CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)"\
						 -O -outdir . 
MEX_EXT = $(shell mexext)
SRC_DIR = .

.PHONY: all
all: undistort dapriltag

dapriltag: $(SRC_DIR)/dapriltag.cpp
	$(MEX) $(MEX_OPTION) $< -L./ -lapriltags -I./  $(CFLAGS_EIGEN) $(CFLAGS_OPENCV)

undistort: $(SRC_DIR)/undistort.cpp
	$(MEX) $(MEX_OPTION) $< $(CFLAGS_OPENCV)

.PHONY: clean
clean:
	rm -f *.$(MEX_EXT)
