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

field_omp2: src/field_omp2.c
	gcc src/field_omp2.c -lSDL2 -lSDL2main -fopenmp -O2 -o bin/field_omp2
	./bin/field_omp2

move_omp: src/move_omp.c
	gcc src/move_omp.c -lSDL2 -lSDL2main -fopenmp -O2 -o bin/move_omp
	./bin/move_omp


