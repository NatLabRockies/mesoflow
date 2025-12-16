MFLO_HOME ?= ../../../
PRECISION  = DOUBLE
PROFILE    = FALSE
SUNDIALS_DIR?=/Users/hsitaram/projects/sundials

DEBUG      = TRUE
DEBUG      = FALSE

DIM       = 3

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE 
USE_CUDA   = FALSE
USE_SUNDIALS = FALSE

#USE_SEC_ORDER_FLUX = TRUE

Bpack   := ./Make.package 
Blocs   := . 

ifeq ($(USE_SUNDIALS), TRUE)
    INCLUDE_LOCATIONS += ${SUNDIALS_DIR}/include
    LIBRARIES += -L$(SUNDIALS_DIR)/lib64 -lsundials_cvode -lsundials_arkode
    LIBRARIES += -L$(SUNDIALS_DIR)/lib64 -lsundials_core -lsundials_nvecserial -lsundials_nvecmanyvector
endif

include ${MFLO_HOME}/tests/Make.mflo
