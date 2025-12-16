MFLO_HOME ?= ../../
PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM        = 3

COMP       = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA   = FALSE

Bpack   := ./Make.package
Blocs   := .

USE_SUNDIALS = FALSE
SUNDIALS_DIR?=/projects/co2rue/sundials

ifeq ($(USE_SUNDIALS),TRUE)
    INCLUDE_LOCATIONS += ${SUNDIALS_DIR}/include
    LIBRARIES += -L$(SUNDIALS_DIR)/lib64 -lsundials_cvode -lsundials_arkode
    LIBRARIES += -L$(SUNDIALS_DIR)/lib64 -lsundials_core -lsundials_nvecserial -lsundials_nvecmanyvector
endif

include ${MFLO_HOME}/models/Make.mflo
