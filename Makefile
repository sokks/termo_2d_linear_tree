COMPILER=mpicxx
BASE_GRID_SIZE=64
TIME_STEPS=200
N_PROCS=2

all: bin


bin: bin/test

translate: bin/translate
	./translate.sh

visualize: vis_2d.py
	python3 vis_2d.py $(BASE_GRID_SIZE) $(BASE_GRID_SIZE) data/temp temp.mp4
	python3 vis_2d.py $(BASE_GRID_SIZE) $(BASE_GRID_SIZE) data/refine refine.mp4 gist_ncar

test: bin/test
	rm -rf data/temp/*
	rm -rf data/refine/*
	bin/test $(BASE_GRID_SIZE) $(TIME_STEPS)

run_mpi: bin/test
	rm -rf data/temp/*
	rm -rf data/refine/*
	mpiexec -np $(N_PROCS) bin/test $(BASE_GRID_SIZE) $(TIME_STEPS)

bin/test: build/main.o build/area.o build/proc.o Makefile
	mkdir -p bin
	$(COMPILER) -std=c++11 -o $@ build/main.o build/area.o build/proc.o


bin/translate: build/translate.o build/proc.o build/area.o Makefile
	$(COMPILER) -std=c++11 -o $@ build/translate.o build/proc.o build/area.o

build/translate.o: src/translate.cpp Makefile
	mkdir -p build
	$(COMPILER) -std=c++11 -o $@ -c src/translate.cpp

build/main.o: src/main.cpp Makefile
	mkdir -p build
	$(COMPILER) -std=c++11 -o $@ -c src/main.cpp

build/proc.o: src/proc/proc.h src/proc/proc.cpp Makefile
	mkdir -p build
	$(COMPILER) -std=c++11 -o $@ -c src/proc/proc.cpp 

build/area.o: src/area/area.h src/area/area.cpp Makefile
	mkdir -p build
	$(COMPILER) -std=c++11 -o $@ -c src/area/area.cpp


clean:
	rm -rf build bin lib
	find . -name \*~ -delete


