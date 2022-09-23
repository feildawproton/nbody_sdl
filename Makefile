basecase: src/basecase.c
	#gcc src/basecase.c -lSDL2 -lSDL2main -pg -o bin/basecase
	gcc src/basecase.c -lSDL2 -lSDL2main -O2 -o bin/basecase
	./bin/basecase
	#gprof bin/basecase gmon.out > prof/prof_basecase

basecase_omp: src/basecase_omp.c
	#gcc src/basecase_omp.c -lSDL2 -lSDL2main -fopenmp -pg -o bin/basecase_omp
	gcc src/basecase_omp.c -lSDL2 -lSDL2main -fopenmp -o bin/basecase_omp
	./bin/basecase_omp
	#gprof bin/basecase_omp gmon.out > prof/prof_basecase_omp

basecase_omp2: src/basecase_omp2.c
	gcc src/basecase_omp2.c -lSDL2 -lSDL2main -fopenmp -O2 -o bin/basecase
	./bin/basecase

