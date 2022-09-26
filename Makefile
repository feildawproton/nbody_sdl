basecase: src/basecase.c
	#gcc src/basecase.c -lSDL2 -lSDL2main -pg -o bin/basecase
	gcc src/basecase.c -lSDL2 -lSDL2main -O2 -o bin/basecase
	./bin/basecase
	#gprof bin/basecase gmon.out > prof/prof_basecase

basecase_omp: src/basecase_omp.c
	gcc src/basecase_omp.c -lSDL2 -lSDL2main -fopenmp -O2 -o bin/basecase_omp
	./bin/basecase_omp

field: src/field.c
	gcc src/field.c -lSDL2 -lSDL2main -O2 -o bin/field
	./bin/field

field_omp: src/field_omp.c
	gcc src/field_omp.c -lSDL2 -lSDL2main -fopenmp -O2 -o bin/field_omp
	./bin/field_omp

pull_omp: src/pull_omp.c
	gcc src/pull_omp.c -lSDL2 -lSDL2main -fopenmp -O2 -o bin/pull_omp
	./bin/pull_omp

pull2: src/pull2.c
	gcc src/pull2.c -lSDL2 -lSDL2main -fopenmp -O2 -o bin/pull2
	./bin/pull2


