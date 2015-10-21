// SET SIZE OF COARSEST GRID: 0 = 3x3, 1=5x5 ...
#define COARSEST 1 

// SET PRECISION OF CALCULATION 
#define DOUBLE_PRECISION 0

// SET FILEOUTPUT: 0 = NO PRINT, 1 = PRINT FOR GNUPLOT, 2 = PRINT FOR TOOL
#define PRINTMODE 2 

// SET NUMBER OF PRESMOOTHING RBGS-CYCLES
#define PRE 2 

// SET NUMBER OF POSTSMOOTHING RBGS-CYCLES
#define POST 1

// SET NUMBER OF MULTIGRID-CYCLES
#define FLOATCYCLES 5
#define MAXCYCLES 7

// use CG AS SOLVER FOR COARSEST GRID
#define USE_CG 0 

// SET NUMBER OF CG-ITERATIONS FOR SOLVING COARSEST GRID
#define CGSTEPS 10 

// PRINT TOTAL EXECUTION TIME OF SUBROUTINES
#define MEASURE_TIME_PER_FUNCTION 0

// PRINT L2NORM OF RESIDUAL OF EVERY CYCLE
#define MEASUREL2NORM 0

// PRINT ERROR OF CALCULATION IN L2-NORM
#define CALCERROR 1 

#define BLOCKSIZE_RESIDUAL 1337 
#define MINHALFSIZE 2048

//USE UNROLLED RBGS FOR SOLVING
//#define USEUNROLLED

