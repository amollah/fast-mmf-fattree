export OMP_NESTED=TRUE
export OMP_NUM_THREADS=1

all:

#	gcc -O3 -o mmf_nlp_omp.x  -fopenmp -mcmodel=large \
#	xgft_tr.c  \
#	xgft_utils.c  \
#	model_utils.c \
#	mmf_engine_nonlp.c  \
#	cplex_driver_mmf.c
#	mmf_driver_nonlp.c

	gcc -O3 -o fast_mmf.out -fopenmp  -mcmodel=large \
	xgft_tr.c \
	model_utils.c \
	xgft_utils.c  \
	cplex_engine_mmf.c  \
	mmf_engine_nonlp.c  \
	cplex_driver_mmf.c

small:
#       gcc -O3 -o cplex_gen.x  -D_SMALL  -mcmodel=large -I../topology ../topology/dragonfly_tr.c ../topology/xg$
	gcc -O3 -o cplex_gen.x  -D_SMALL  -mcmodel=large -I../topology ../topology/xgft_tr.c cplex_engine_mmf.c cplex_driver_mmf.c

clean:
	rm -f *~
	rm -f *.o
	rm -f \#*\#
	rm -f *.out
	rm -f *.x
	rm -f *.log
	rm -f *.lp
