#=======================================================================================
# Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
# All rights reserved.
# Use of this source code is governed by a MIT-style
# license that can be found in the LICENSE file.
#=======================================================================================

#CONFIGURE BUILD SYSTEM
TARGET	   = exe-$(TAG)
BUILD_DIR  = ./$(TAG)
SRC_DIR    = ./src
MAKE_DIR   = ./
Q         ?= @

#DO NOT EDIT BELOW
include $(MAKE_DIR)/config.mk
include $(MAKE_DIR)/include_$(TAG).mk
INCLUDES  += -I$(SRC_DIR) -I$(BUILD_DIR)

VPATH     = $(SRC_DIR)
SRC       = $(filter-out $(wildcard $(SRC_DIR)/*-*.c),$(wildcard $(SRC_DIR)/*.c))
ASM       = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.s, $(SRC))
OBJ       = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRC))
OBJ      += $(BUILD_DIR)/vtkWriter-$(VTK_OUTPUT_FMT).o
OBJ      += $(BUILD_DIR)/solver-$(SOLVER).o
SOURCES   = $(SRC) $(wildcard $(SRC_DIR)/*.h)
ifeq ($(VTK_OUTPUT_FMT),mpi)
DEFINES  += -D_VTK_WRITER_MPI
endif

CPPFLAGS := $(CPPFLAGS) $(DEFINES) $(OPTIONS) $(INCLUDES)

${TARGET}: sanity-checks $(BUILD_DIR) $(OBJ)
	$(info ===>  LINKING  $(TARGET))
	$(Q)${LINKER} ${LFLAGS} -o $(TARGET) $(OBJ) $(LIBS)

$(BUILD_DIR)/%.o:  %.c $(MAKE_DIR)/include_$(TAG).mk $(MAKE_DIR)/config.mk
	$(info ===>  COMPILE  $@)
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(Q)$(GCC) $(CPPFLAGS) -MT $(@:.d=.o) -MM  $< > $(BUILD_DIR)/$*.d

$(BUILD_DIR)/%.s:  %.c
	$(info ===>  GENERATE ASM  $@)
	$(CC) -S $(CPPFLAGS) $(CFLAGS) $< -o $@

.PHONY: clean distclean vis vis_clean tags info asm format

vis:
	$(info ===>  GENERATE VISUALIZATION)
	@gnuplot -e "filename='residual.dat'" ./residual.plot

vis_clean:
	$(info ===>  CLEAN VISUALIZATION)
	@rm -f *.dat
	@rm -f *.vtk
	@rm -f *.png

clean: vis_clean
	$(info ===>  CLEAN)
	@rm -rf $(BUILD_DIR)
	@rm -f tags

distclean: clean
	$(info ===>  DIST CLEAN)
	@rm -f $(TARGET)

info:
	$(info $(CFLAGS))
	$(Q)$(CC) $(VERSION)

asm:  $(BUILD_DIR) $(ASM)

tags:
	$(info ===>  GENERATE TAGS)
	$(Q)ctags -R

format:
	@for src in $(SOURCES) ; do \
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
	@mkdir $(BUILD_DIR)

-include $(OBJ:.o=.d)
