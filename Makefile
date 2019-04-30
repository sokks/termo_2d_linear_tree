# COMPILER=mpixlC
# COMPILER=mpicxx
COMPILER=mpixlcxx_r
OPTS=-O0\ -g
#OPTS=-O0\ -qlanglvl=extended0x
PYTHON_I=python3
# PYTHON_I=python

BASE_LVL=7
MAX_LVL=10

N_PROCS=8

TIME_STEPS=100

all: bin

bin: bin/test



translate: bin/translate
	./translate.sh



test: bin/test
	rm -rf data/temp/*
	bin/test $(BASE_LVL) $(MAX_LVL) $(TIME_STEPS)

gen_grid: bin/gen_grid
	rm -rf data/refine/*
	mkdir -p data/refine
	bin/gen_grid $(BASE_LVL) $(MAX_LVL) data/refine/base_grid.dat $(N_PROCS) data/refine/offsets_$(N_PROCS).dat

polus_job_gen_grid: bin/gen_grid
	rm -rf data/refine/*
	mkdir -p data/refine
	mpisubmit.pl -p 1 bin/gen_grid -- $(BASE_LVL) $(MAX_LVL) data/refine/base_grid.dat $(N_PROCS) data/refine/offsets_$(N_PROCS).dat

bg_job_gen_grid: bin/gen_grid
	rm -rf data/refine/*
	mkdir -p data/refine
	mpisubmit.bg -n 1 -m smp bin/gen_grid -- $(BASE_LVL) $(MAX_LVL) data/refine/base_grid.dat $(N_PROCS) data/refine/offsets_$(N_PROCS).dat

vis_base_grid: update_txt
	$(PYTHON_I) vis_2d_nonuniform.py $(MAX_LVL) data/refine/base_grid.txt data/pics/grid_levels_$(BASE_LVL).png Greys lvls
	
vis_decomposition: update_txt
	$(PYTHON_I) vis_2d_nonuniform.py $(MAX_LVL) data/refine/base_grid.txt data/pics/grid_decomposition_$(BASE_LVL)_$(N_PROCS).png prism procs $(N_PROCS)

vis_temps: translate
	# python3 vis_2d_nonuniform.py $(MAX_LVL) data/refine/base_grid.txt data/pics/start_temp.png coolwarm temp
	# python3 vis_2d_nonuniform.py $(MAX_LVL) data/temp/000200.out.txt data/pics/000200.png coolwarm temp
	mkdir -p data/temp/img
	./plot_temps.sh $(MAX_LVL)

vis_start_temp: update_txt
	$(PYTHON_I) vis_2d_nonuniform.py $(MAX_LVL) data/refine/base_grid.txt data/pics/start_temp.png coolwarm temp

update_txt: bin/translate
	bin/translate data/refine/base_grid.dat data/refine/base_grid.txt



run_mpi: bin/test
	rm -rf data/temp/*
	mkdir -p data/temp
	mpiexec -np $(N_PROCS) bin/test $(BASE_LVL) $(MAX_LVL) data/refine/offsets_$(N_PROCS).dat data/refine/base_grid.dat $(TIME_STEPS)

polus_job_run_mpi: bin/test
	rm -rf data/temp/*
	mkdir -p data/temp
	mpisubmit.pl -p $(N_PROCS) bin/test -- $(BASE_LVL) $(MAX_LVL) data/refine/offsets_$(N_PROCS).dat data/refine/base_grid.dat $(TIME_STEPS)

bg_job_run_mpi: bin/test
	rm -rf data/temp/*
	mkdir -p data/temp
	mpisubmit.bg -n $(N_PROCS) -m smp bin/test -- $(BASE_LVL) $(MAX_LVL) data/refine/offsets_$(N_PROCS).dat data/refine/base_grid.dat $(TIME_STEPS)


bin/test: build/main.o build/area.o build/grid.o build/proc.o Makefile
	mkdir -p bin
	$(COMPILER) -std=c++11 -o $@ build/main.o build/area.o build/grid.o build/proc.o

build/main.o: src/main.cpp Makefile
	mkdir -p build
	$(COMPILER) -std=c++11 $(OPTS) -o $@ -c src/main.cpp


bin/gen_grid: build/gen_grid.o build/grid.o build/area.o Makefile
	mkdir -p bin
	$(COMPILER) -std=c++11 $(OPTS) -o $@ build/gen_grid.o build/grid.o build/area.o

build/gen_grid.o: src/gen_grid.cpp Makefile
	mkdir -p build
	$(COMPILER) -std=c++11 $(OPTS) -o $@ -c src/gen_grid.cpp


bin/translate: build/translate.o Makefile
	$(COMPILER) -std=c++11 $(OPTS) -o $@ build/translate.o

build/translate.o: src/translate.cpp Makefile
	mkdir -p build
	$(COMPILER) -std=c++11 $(OPTS) -o $@ -c src/translate.cpp


build/grid.o: src/grid/grid.h src/grid/grid.cpp Makefile
	mkdir -p build
	$(COMPILER) -std=c++11 $(OPTS) -o $@ -c src/grid/grid.cpp 

build/proc.o: src/proc/proc.h src/proc/proc.cpp Makefile
	mkdir -p build
	$(COMPILER) -std=c++11 $(OPTS) -o $@ -c src/proc/proc.cpp 

build/area.o: src/area/area.h src/area/area.cpp Makefile
	mkdir -p build
	$(COMPILER) -std=c++11 $(OPTS) -o $@ -c src/area/area.cpp


clean:
	rm -rf build bin lib
	find . -name \*~ -delete


