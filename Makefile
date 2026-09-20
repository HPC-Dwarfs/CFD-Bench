# Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
# All rights reserved. This file is part of CFD-Bench.
# Use of this source code is governed by a MIT-style
# license that can be found in the LICENSE file.

#CONFIGURE BUILD SYSTEM
TARGET	   = CFD-Bench-$(TOOLCHAIN)
BUILD_DIR  = ./build/$(TOOLCHAIN)
SRC_DIR    = ./src
MAKE_DIR   = ./mk
Q         ?= @

#DO NOT EDIT BELOW
ifeq (,$(wildcard config.mk))
$(info )
$(info ====================================================================)
$(info config.mk does not exist!)
$(info Creating config.mk from ./mk/config-default.mk)
$(info Please adapt config.mk to your needs and run make again.)
$(info ====================================================================)
$(info )
$(shell cp ./mk/config-default.mk config.mk)
$(error Stopping after creating config.mk - please review and run make again)
endif
include config.mk
include $(MAKE_DIR)/include_$(TOOLCHAIN).mk
INCLUDES  += -I$(SRC_DIR) -I$(BUILD_DIR)

# Which solver variant this build directory holds.
#
# The Makefile links exactly one solver-$(SOLVER).o and compiles everything with
# -DSOLVER_$(SOLVER), but none of that is visible to make as a timestamp. So
# switching back to a variant whose object file is already up to date used to
# leave the previously linked binaries in place -- tests/run-all.sh loops over
# every solver and so hits this on its second run, silently reporting one
# variant's results under another's name. Rewriting this file whenever the
# selection changes makes the choice an ordinary prerequisite.
$(shell mkdir -p $(BUILD_DIR))
$(shell [ "$$(cat $(BUILD_DIR)/solver.sel 2>/dev/null)" = "$(SOLVER)" ] || printf '%s' "$(SOLVER)" > $(BUILD_DIR)/solver.sel)
SOLVER_SEL = $(BUILD_DIR)/solver.sel

VPATH     = $(SRC_DIR)
OBJ       = $(filter-out $(BUILD_DIR)/vtkWriter-%.o $(BUILD_DIR)/solver-%.o, $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(wildcard $(SRC_DIR)/*.c)))
OBJ      += $(BUILD_DIR)/vtkWriter-$(VTK_OUTPUT_FMT).o
OBJ      += $(BUILD_DIR)/solver-$(SOLVER).o
ASM       = $(patsubst $(BUILD_DIR)/%.o, $(BUILD_DIR)/%.s, $(OBJ))

# Test build path. The same sources are compiled a second time with -DTEST into
# their own object directory, so the normal build is never affected. That build
# produces the solver with the field dump enabled, plus one binary per check
# driver in $(CHECK_DIR).
TEST_BUILD_DIR  = $(BUILD_DIR)/test
TEST_TARGET     = $(TARGET)-test
CHECK_DIR       = ./tests/checks
CHECK_SRC       = $(wildcard $(CHECK_DIR)/*.c)
CHECK_BIN       = $(patsubst $(CHECK_DIR)/%.c, $(BUILD_DIR)/check-%, $(CHECK_SRC))
TEST_OBJ        = $(patsubst $(BUILD_DIR)/%.o, $(TEST_BUILD_DIR)/%.o, $(OBJ))
TEST_OBJ_NOMAIN = $(filter-out $(TEST_BUILD_DIR)/main.o, $(TEST_OBJ))
TOOL_BIN        = tools/fieldcmp

ifeq ($(VTK_OUTPUT_FMT),mpi)
DEFINES  += -D_VTK_WRITER_MPI
endif
# Lets a check driver tell which solver variant it was linked against, so a
# driver that exercises one solver's internals can skip under the others.
DEFINES  += -DSOLVER_$(SOLVER)
SRC       =  $(wildcard $(SRC_DIR)/*.h $(SRC_DIR)/*.c)
CPPFLAGS := $(CPPFLAGS) $(DEFINES) $(OPTIONS) $(INCLUDES)
c := ,
clist = $(subst $(eval) ,$c,$(strip $1))

define CLANGD_TEMPLATE
CompileFlags:
  Add: [$(call clist,$(CPPFLAGS)), $(call clist,$(CFLAGS)), -xc]
  Compiler: clang
endef

${TARGET}: sanity-checks $(BUILD_DIR) .clangd $(SOLVER_SEL) $(OBJ)
	$(info ===>  LINKING  $(TARGET))
	$(Q)${LD} ${LFLAGS} -o $(TARGET) $(OBJ) $(LIBS)

$(BUILD_DIR)/%.o:  %.c $(MAKE_DIR)/include_$(TOOLCHAIN).mk config.mk
	$(info ===>  COMPILE  $@)
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(Q)$(CC) $(CPPFLAGS) -MT $(@:.d=.o) -MM  $< > $(BUILD_DIR)/$*.d

$(BUILD_DIR)/%.s:  %.c
	$(info ===>  GENERATE ASM  $@)
	$(CC) -S $(CPPFLAGS) $(CFLAGS) $< -o $@

$(TEST_BUILD_DIR)/%.o:  %.c $(MAKE_DIR)/include_$(TOOLCHAIN).mk config.mk
	$(info ===>  COMPILE (test)  $@)
	$(Q)$(CC) -c $(CPPFLAGS) -DTEST $(CFLAGS) $< -o $@
	$(Q)$(CC) $(CPPFLAGS) -DTEST -MT $@ -MM $< > $(TEST_BUILD_DIR)/$*.d

$(TEST_TARGET): sanity-checks $(TEST_BUILD_DIR) $(SOLVER_SEL) $(TEST_OBJ)
	$(info ===>  LINKING  $(TEST_TARGET))
	$(Q)${LD} ${LFLAGS} -o $(TEST_TARGET) $(TEST_OBJ) $(LIBS)

$(BUILD_DIR)/check-%: $(CHECK_DIR)/%.c $(SOLVER_SEL) $(TEST_OBJ_NOMAIN)
	$(info ===>  LINKING  $@)
	$(Q)$(CC) $(CPPFLAGS) -I$(CHECK_DIR) -DTEST $(CFLAGS) -o $@ $< $(TEST_OBJ_NOMAIN) ${LFLAGS} $(LIBS)

tools/fieldcmp: tools/fieldcmp.c
	$(info ===>  LINKING  $@)
	$(Q)$(CC) $(CFLAGS) -o $@ $< -lm

tests: $(TEST_TARGET) $(CHECK_BIN) $(TOOL_BIN)

.PHONY: clean distclean info asm format plot tests

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
	@rm -f $(TARGET) $(TEST_TARGET) $(TOOL_BIN)
	@rm -f tags .clangd compile_commands.json
	@rm -f profile-*.txt
	@rm -f *.dat
	@rm -f *.vtk
	@rm -f *.png

info:
	$(info CFD-Bench v1.0)
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

$(TEST_BUILD_DIR):
	@mkdir -p $(TEST_BUILD_DIR)

.clangd:
	$(file > .clangd,$(CLANGD_TEMPLATE))

-include $(OBJ:.o=.d)
-include $(TEST_OBJ:.o=.d)
