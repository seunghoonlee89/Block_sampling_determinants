#Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
#Copyright (c) 2012, Garnet K.-L. Chan                                                                                                                     
#This program is integrated in Molpro with the permission of 
#Sandeep Sharma and Garnet K.-L. Chan


######### GENERAL OPTIONS FOR USER #########
#specify boost include file
BOOSTINCLUDE = /usr/local/apps/gcc-5.4.0/boost_1_55_0/include

#specify boost and lapack-blas library locations
BOOSTLIB = -L/usr/local/apps/gcc-5.4.0/boost_1_55_0/lib  -lboost_system-mt -lboost_filesystem-mt -lboost_serialization-mt
#BOOSTLIB = -lboost_serialization -lboost_system -lboost_filesystem

# set if using BOOST version >= 1.56.0
USE_BOOST56 = no
ifeq ($(USE_BOOST56), yes)
	B56 = -DBOOST_1_56_0
endif

#use these variable to set if we will use mpi or not 
USE_MPI = yes
USE_MKL = yes

# change to icpc for Intel
CXX =  g++
MPICXX = mpic++
export CXX
export MPICXX

# use this variable to set if we will use integer size of 8 or not.
# molpro compilation requires I8, since their integers are long
I8_OPT = no

MOLPRO = no
DOPROF = no

# add Molcas interface to libqcdmrg.so
# molcas compilation w/ -64 option requires I8 as well
MOLCAS = no

OPENMP = no

ifeq ($(USE_MKL), yes)
MKLROOT = /usr/local/apps/intel/compilers_and_libraries_2017.0.098/linux/mkl
ifeq ($(OPENMP), yes)
LAPACKBLAS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
else
LAPACKBLAS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
endif
MKLFLAGS = ${MKLROOT}/include/
MKLOPT = -D_HAS_INTEL_MKL
else
MKLFLAGS = .
endif

RUN_UNITTEST=no
ifeq ($(RUN_UNITTEST), yes)
UNITTEST = -DUNITTEST
endif

AR=ar
ARFLAGS=-qs
RANLIB=ranlib

ifeq ($(I8_OPT), yes)
	I8 = -D_FORTINT_64
endif

EXECUTABLE = block.spin_adapted

BLOCKHOME = .
HOME = .
NEWMATINCLUDE = $(BLOCKHOME)/newmat10/
INCLUDE1 = $(BLOCKHOME)/include/
INCLUDE2 = $(BLOCKHOME)/
NEWMATLIB = $(BLOCKHOME)/newmat10/
BTAS = $(BLOCKHOME)/btas
.SUFFIXES: .C .cpp

ifeq ($(MOLPRO), yes)
   MOLPRO_BLOCK= -DMOLPRO
endif

ifeq ($(MOLCAS), yes)
   MOLCAS_BLOCK= -DMOLCAS
endif

FLAGS =  -I${MKLFLAGS} -I$(INCLUDE1) -I$(INCLUDE2) -I$(NEWMATINCLUDE) -I$(BOOSTINCLUDE) \
         -I$(HOME)/modules/generate_blocks/ -I$(HOME)/modules/onepdm -I$(HOME)/modules/twopdm/ \
         -I$(HOME)/modules/npdm -I$(HOME)/modules/two_index_ops -I$(HOME)/modules/three_index_ops -I$(HOME)/modules/four_index_ops -std=c++0x \
	 -I$(HOME)/modules/ResponseTheory -I$(HOME)/modules/nevpt2 -I$(HOME)/molcas -I$(HOME)/modules/mps_nevpt \
         -I$(HOME)/modules/ds1_onepdm -I$(HOME)/modules/ds0_onepdm

LIBS +=  -L$(NEWMATLIB) -lnewmat $(BOOSTLIB) $(LAPACKBLAS) 
MPI_OPT = -DSERIAL

ifeq (icpc, $(CXX))
   ifeq ($(OPENMP), yes)
      OPENMP_FLAGS= -openmp -D_OPENMP 
   endif
# Intel compiler
   OPT = -DNDEBUG -O3 -funroll-loops -Werror
#  OPT = -g -fPIC
   ifeq ($(USE_MPI), no) 
      CXX = icc
   endif
endif

# GNU compiler
ifeq (g++, $(CXX))
   ifeq ($(OPENMP), yes)
      OPENMP_FLAGS= -fopenmp #-D_OPENMP 
   endif
   OPT = -DNDEBUG -O3 -g -funroll-loops 
#   OPT = -g -fPIC
endif

ifeq (clang++, $(CXX))
   ifeq ($(OPENMP), yes)
      OPENMP_FLAGS= -fopenmp #-D_OPENMP 
   endif

   OPT = -DNDEBUG -g -Werror
endif

ifeq ($(DOPROF),yes)
     OPT += -pg
endif

ifeq ($(USE_MPI), yes)
     MPI_OPT = 
     MPI_LIB = -lboost_mpi -L/opt/intel/compilers_and_libraries_2016.0.109/linux/mpi/intel64/lib -lmpi 
     LIBS += $(MPI_LIB)
     CXX = $(MPICXX)
endif


OPT	+= $(OPENMP_FLAGS) -DBLAS -DUSELAPACK $(MPI_OPT) $(I8) $(B56) $(MOLPRO_BLOCK) $(MOLCAS_BLOCK) -DFAST_MTP -D_HAS_CBLAS -D_HAS_INTEL_MKL ${MKLOPT} ${UNITTEST} -fPIC

SRC_genetic = genetic/CrossOver.C genetic/Evaluate.C genetic/GAInput.C genetic/GAOptimize.C genetic/Generation.C genetic/Mutation.C genetic/RandomGenerator.C genetic/ReadIntegral.C

SRC_ds  = modules/ds1_onepdm/ds1_sweep.C modules/ds1_onepdm/ds1_onepdm.C modules/ds0_onepdm/ds0_sweep.C modules/ds0_onepdm/ds0_onepdm.C

SRC_npdm = modules/npdm/npdm.C modules/npdm/npdm_driver.C modules/npdm/npdm_patterns.C modules/npdm/npdm_expectations.C modules/npdm/npdm_expectations_engine.C  \
           modules/npdm/npdm_permutations.C modules/npdm/npdm_spin_adaptation.C modules/npdm/npdm_operator_selector.C modules/npdm/npdm_spin_ops.C \
           modules/npdm/onepdm_container.C modules/npdm/twopdm_container.C modules/npdm/threepdm_container.C modules/npdm/fourpdm_container.C  \
           modules/two_index_ops/two_index_wrappers.C modules/three_index_ops/three_index_wrappers.C modules/four_index_ops/four_index_wrappers.C  \
           modules/three_index_ops/three_index_compound_ops.C modules/four_index_ops/four_index_compound_ops.C  \
           modules/three_index_ops/three_index_op_components.C modules/four_index_ops/four_index_op_components.C  \
           modules/three_index_ops/three_index_ops.C modules/four_index_ops/four_index_ops.C  \
           modules/three_index_ops/build_3index_ops.C modules/four_index_ops/build_4index_ops.C modules/npdm/pairpdm_container.C \
           modules/npdm/nevpt2_npdm_driver.C modules/npdm/nevpt2_container.C modules/npdm/externalsort.C

#SRC_mps_nevpt = module/mps_nevpt/nevpt_Operators.C module/mps_nevpt/nevpt_op_component.C module/mps_nevpt/type1.C module/mps_nevpt/mps_nevpt.C
SRC_mps_nevpt = modules/mps_nevpt/nevpt_op_component.C modules/mps_nevpt/type1.C modules/mps_nevpt/mps_nevpt.C

SRC_spin_adapted =  modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C main.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C hfOccGenerator.C Schedule.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2) $(SRC_mps_nevpt)  $(SRC_ds)

SRC_OH = modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C wrapper.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C Schedule.C input.C hfOccGenerator.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2) $(SRC_mps_nevpt)  $(SRC_ds)
OBJ_OH= OverlapHelement.o

SRC_STOPT = modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C wrapper.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C Schedule.C input.C hfOccGenerator.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2) $(SRC_mps_nevpt)  $(SRC_ds) heatbath.C nonspinmps.C #Determinants.C
OBJ_STOPT= stochasticpt_new.o

SRC_MPS2CI_det = modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C wrapper.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C Schedule.C input.C hfOccGenerator.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2) $(SRC_mps_nevpt)  $(SRC_ds) heatbath.C nonspinmps.C #Determinants.C
OBJ_MPS2CI_det = MPS2CI_det.o

SRC_MPS2CI_ecCC = modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C wrapper.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C Schedule.C input.C hfOccGenerator.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2) $(SRC_mps_nevpt)  $(SRC_ds) heatbath.C nonspinmps.C #Determinants.C
OBJ_MPS2CI_ecCC = MPS2CI_ecCC.o

SRC_CITRIE = modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C wrapper.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C Schedule.C input.C hfOccGenerator.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2) $(SRC_mps_nevpt)  $(SRC_ds) heatbath.C nonspinmps.C #Determinants.C
OBJ_CITRIE = CIcoeff_from_mps_TRIE.o

SRC_MPS2CIdbg = modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C wrapper.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C Schedule.C input.C hfOccGenerator.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2) $(SRC_mps_nevpt)  $(SRC_ds) heatbath.C nonspinmps.C #Determinants.C
OBJ_MPS2CIdbg = MPS2CIdbg.o

SRC_CSFOH = CSFOverlapHelement.C modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C wrapper.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C Schedule.C input.C hfOccGenerator.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2)  $(SRC_mps_nevpt)  $(SRC_ds)

SRC_COEF= modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C wrapper.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C Schedule.C input.C hfOccGenerator.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2) $(SRC_mps_nevpt)  $(SRC_ds)
OBJ_COEF= Coefficients.o

SRC_spin_library =  modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C hfOccGenerator.C Schedule.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2)  $(SRC_mps_nevpt)  $(SRC_ds)

SRC_nevpt2 = modules/nevpt2/nevpt2.C modules/nevpt2/nevpt2_info.C modules/nevpt2/nevpt2_mpi.C \
             modules/nevpt2/nevpt2_opconstruct.C modules/nevpt2/nevpt2_operators.C \
             modules/nevpt2/nevpt2_pal.C modules/nevpt2/nevpt2_renormalize.C \
             modules/nevpt2/nevpt2_util.C modules/nevpt2/ripdm.C modules/nevpt2/sweep_gen_nevpt2.C \
             modules/nevpt2/sweep_nevpt2.C

SRC_molcas = molcas/block_calldmrg.C molcas/molpro_fcidump.C molcas/loadNpdm.C molcas/sortNpdm.C molcas/tranNpdm.C


OBJ_OH+=$(SRC_OH:.C=.o)
OBJ_STOPT+=$(SRC_STOPT:.C=.o)
OBJ_MPS2CI_det+=$(SRC_MPS2CI_det:.C=.o)
OBJ_MPS2CI_ecCC+=$(SRC_MPS2CI_ecCC:.C=.o)
OBJ_CITRIE+=$(SRC_CITRIE:.C=.o)
OBJ_MPS2CIdbg+=$(SRC_MPS2CIdbg:.C=.o)
OBJ_CSFOH+=$(SRC_CSFOH:.C=.o)
OBJ_COEF+=$(SRC_COEF:.C=.o)
OBJ_spin_adapted=$(SRC_spin_adapted:.C=.o)
OBJ_spin_library=$(SRC_spin_library:.C=.o)
OBJ_nevpt2=$(SRC_nevpt2:.C=.o)
OBJ_mps_nevpt=$(SRC_mps_nevpt:.C=.o)
ifeq ($(MOLCAS), yes)
	OBJ_spin_library+=$(SRC_molcas:.C=.o)
endif

.C.o :
	$(CXX)  $(FLAGS) $(OPT) -c $< -o $@
.cpp.o :
	$(CXX) $(FLAGS) $(OPT) -c $< -o $@

all	: library $(EXECUTABLE) OH COEF CSFOH STOPT MPS2CI_det MPS2CI_ecCC MPS2CIdbg CITRIE 
#all	: library $(EXECUTABLE) OH COEF CSFOH STOPT CITRIE

ifeq ($(MOLCAS), yes)
library : libqcdmrg.a $(NEWMATLIB)/libnewmat.a libqcdmrg.so
else                                           
library : libqcdmrg.a $(NEWMATLIB)/libnewmat.a 
endif

libqcdmrg.a : $(OBJ_spin_library) $(OBJ_molcas) $(OBJ_mps_nevpt) $(NEWMATLIB)/libnewmat.a
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

libqcdmrg.so : $(OBJ_spin_library) $(OBJ_molcas) $(OBJ_mps_nevpt)

	$(CXX) -shared -o $@ $^ $(LIBS)

$(EXECUTABLE) : $(OBJ_spin_adapted) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  $(EXECUTABLE) $(OBJ_spin_adapted) $(LIBS)

OH : $(OBJ_OH) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  OH $(OBJ_OH) $(LIBS)

STOPT : $(OBJ_STOPT) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  STOPT $(OBJ_STOPT) $(LIBS)

CITRIE : $(OBJ_CITRIE) $(NEWMATLIB)/libnewmat.a
	$(CXX) -O3  $(FLAGS) $(OPT) -o  CITRIE $(OBJ_CITRIE) $(LIBS)
#	$(CXX) -g -O3  $(FLAGS) $(OPT) -o  CITRIE $(OBJ_CITRIE) $(LIBS)

MPS2CI_det : $(OBJ_MPS2CI_det) $(NEWMATLIB)/libnewmat.a
	$(CXX) -O3  $(FLAGS) $(OPT) -o  MPS2CI_det $(OBJ_MPS2CI_det) $(LIBS)
#	$(CXX) -g -O3  $(FLAGS) $(OPT) -o  MPS2CI_det $(OBJ_MPS2CI_det) $(LIBS)

MPS2CI_ecCC : $(OBJ_MPS2CI_ecCC) $(NEWMATLIB)/libnewmat.a
	$(CXX) -O3  $(FLAGS) $(OPT) -o  MPS2CI_ecCC $(OBJ_MPS2CI_ecCC) $(LIBS)
#	$(CXX) -g -O3  $(FLAGS) $(OPT) -o  MPS2CI $(OBJ_MPS2CI) $(LIBS)

MPS2CIdbg : $(OBJ_MPS2CIdbg) $(NEWMATLIB)/libnewmat.a
	$(CXX)  $(FLAGS) $(OPT) -o  MPS2CIdbg $(OBJ_MPS2CIdbg) $(LIBS)

CSFOH : $(OBJ_CSFOH) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  CSFOH $(OBJ_CSFOH) $(LIBS)

COEF : $(OBJ_COEF) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  COEF $(OBJ_COEF) $(LIBS)

$(NEWMATLIB)/libnewmat.a :
	cd $(NEWMATLIB) &&  $(MAKE) -f makefile libnewmat.a

clean:
	rm *.o include/*.o modules/generate_blocks/*.o modules/onepdm/*.o modules/twopdm/*.o modules/npdm/*.o $(NEWMATLIB)*.o libqcdmrg.a libqcdmrg.so $(EXECUTABLE) $(NEWMATLIB)/libnewmat.a genetic/gaopt genetic/*.o btas/lib/*.o btas/lib/libbtas.a modules/two_index_ops/*.o modules/three_index_ops/*.o modules/four_index_ops/*.o modules/ResponseTheory/*.o modules/nevpt2/*.o molcas/*.o modules/mps_nevpt/*o OH CSFOH COEF modules/ds0_onepdm/*.o modules/ds1_onepdm/*.o MPS2CI_det MPS2CI_ecCC MPS2CIdbg CITRIE STOPT
	#rm *.o include/*.o modules/generate_blocks/*.o modules/onepdm/*.o modules/twopdm/*.o modules/npdm/*.o $(NEWMATLIB)*.o libqcdmrg.a libqcdmrg.so $(EXECUTABLE) $(NEWMATLIB)/libnewmat.a genetic/gaopt genetic/*.o btas/lib/*.o btas/lib/libbtas.a modules/two_index_ops/*.o modules/three_index_ops/*.o modules/four_index_ops/*.o modules/ResponseTheory/*.o modules/nevpt2/*.o molcas/*.o modules/mps_nevpt/*o OH CSFOH COEF modules/ds0_onepdm/*.o modules/ds1_onepdm/*.o CITRIE STOPT

# DO NOT DELETE


