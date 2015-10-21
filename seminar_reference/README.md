# TODO List #

* [Add] CG Solver. This may reduce number of cycles needed and / or may increase the coarsest grid level -> short testing with some parameters did not result in better performance, maybe think of full-multigrid cycle or F cycle for first iteration in order to improve initial guess quickly.
* [Add] SIMD by hand or improve loops in order to get more loops autovectorized by compiler
* [Add] Blocking and therefore maybe padding to matrix-vector operations
* [Add] grid.fill(): replace for loops with memset


# DONE LIST #

* [Fix] Adapt allocated grid space to actually used space -> DONE