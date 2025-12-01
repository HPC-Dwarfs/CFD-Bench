# Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
# All rights reserved.
# Use of this source code is governed by a MIT-style
# license that can be found in the LICENSE file.

#CONFIGURE BUILD SYSTEM
TARGET	   = NusifSolver-$(TOOLCHAIN)
BUILD_DIR  = ./build/$(TOOLCHAIN)
SRC_DIR    = ./src
MAKE_DIR   = ./mk
Q         ?= @

#DO NOT EDIT BELOW
include config.mk
include $(MAKE_DIR)/include_$(TOOLCHAIN).mk
INCLUDES  += -I$(SRC_DIR) -I$(BUILD_DIR)

VPATH     = $(SRC_DIR)
OBJ       = $(filter-out $(BUILD_DIR)/vtkWriter-%.o $(BUILD_DIR)/solver-%.o, $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(wildcard $(SRC_DIR)/*.c)))
OBJ      += $(BUILD_DIR)/vtkWriter-$(VTK_OUTPUT_FMT).o
OBJ      += $(BUILD_DIR)/solver-$(SOLVER).o
ASM       = $(patsubst $(BUILD_DIR)/%.o, $(BUILD_DIR)/%.s, $(OBJ))

ifeq ($(VTK_OUTPUT_FMT),mpi)
DEFINES  += -D_VTK_WRITER_MPI
endif
SRC       =  $(wildcard $(SRC_DIR)/*.h $(SRC_DIR)/*.c)
CPPFLAGS := $(CPPFLAGS) $(DEFINES) $(OPTIONS) $(INCLUDES)
c := ,
clist = $(subst $(eval) ,$c,$(strip $1))

define CLANGD_TEMPLATE
CompileFlags:
  Add: [$(call clist,$(CPPFLAGS)), $(call clist,$(CFLAGS)), -xc]
  Compiler: clang
endef

${TARGET}: sanity-checks $(BUILD_DIR) .clangd $(OBJ)
	$(info ===>  LINKING  $(TARGET))
	$(Q)${LD} ${LFLAGS} -o $(TARGET) $(OBJ) $(LIBS)

$(BUILD_DIR)/%.o:  %.c $(MAKE_DIR)/include_$(TOOLCHAIN).mk config.mk
	$(info ===>  COMPILE  $@)
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(Q)$(CC) $(CPPFLAGS) -MT $(@:.d=.o) -MM  $< > $(BUILD_DIR)/$*.d

$(BUILD_DIR)/%.s:  %.c
	$(info ===>  GENERATE ASM  $@)
	$(CC) -S $(CPPFLAGS) $(CFLAGS) $< -o $@

.PHONY: clean distclean info asm format plot

plot:
	$(info ===>  GENERATE PLOT)
	@gnuplot -e "filename='residual.dat'" ./residual.plot

clean:
	$(info ===>  CLEAN)
	@rm -rf $(BUILD_DIR)

distclean: clean
	$(info ===>  DIST CLEAN)
	@rm -rf build
	@rm -rf .cache
	@rm -f $(TARGET)
	@rm -f tags .clangd compile_commands.json
	@rm -f profile-*.txt
	@rm -f *.dat
	@rm -f *.vtk
	@rm -f *.png

info:
	$(info Nusif-Solver v1.0)
	$(info Flags: $(CFLAGS))
	$(Q)$(CC) $(VERSION)

asm:  $(BUILD_DIR) $(ASM)

format:
	@for src in $(SRC) ; do \
		echo "Formatting $$src" ; \
		clang-format -i $$src ; \
	done
	@echo "Done"

sanity-checks:
ifeq ($(VTK_OUTPUT_FMT),mpi)
ifeq ($(ENABLE_MPI),false)
	$(error VTK_OUTPUT_FMT mpi only supported for ENABLE_MPI true!)
endif
endif

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

.clangd:
	$(file > .clangd,$(CLANGD_TEMPLATE))

-include $(OBJ:.o=.d)
