#  /******************************************************************************
#  *                        ETSI TS 103 634 V1.1.1                               *
#  *              Low Complexity Communication Codec Plus (LC3plus)              *
#  *                                                                             *
#  * Copyright licence is solely granted through ETSI Intellectual Property      *
#  * Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
#  * estoppel or otherwise.                                                      *
#  ******************************************************************************/

# GNU Makefile

# Options
AFL         = 0
CLANG       = 0
GCOV        = 0
KISSFFT     = 0
NO_POST_REL = 0
OPTIM       = 0
PLC         = 1
SUBSET      = 
RELEASE     = PLUS
WMOPS       = 1

# Paths
VPATH  = . basic_op
BUILD  = build
CC     = gcc
LINK   = $(CC)

# Binary Name
NAME_LC3   = LC3plus

# Default tool settings
RM        = rm -f

ifndef VERBOSE
QUIET_CC  = @echo '   ' Compiling $<;
QUIET_LINK= @echo '   ' Linking $@;
QUIET     = @
endif

CFLAGS += -std=c99

# C compiler flags
# Preprocessor(-I/-D) / Compiler / Linker flags
CPPFLAGS += -Ibasic_op -DLC3_$(RELEASE) -DSUBSET_$(SUBSET)
CFLAGS   += -pedantic -Wcast-qual -Wall -W -Wextra -Wno-long-long     \
            -Wpointer-arith -Wstrict-prototypes -Wmissing-prototypes  \
            -Werror-implicit-function-declaration

ifneq "$(DEBUG)" "0"
CFLAGS   += -g3
LDFLAGS  += -g3
endif

LDFLAGS += -lm

DEPFLAGS = -MT $@ -MMD -MP -MF $(BUILD)/$*.Td

ifeq "$(GCOV)" "1"
CFLAGS  += -fprofile-arcs -ftest-coverage
LDFLAGS += -fprofile-arcs -ftest-coverage
endif

OPTIM    ?= 0
CFLAGS   += -O$(OPTIM)

CFLAGS   += $(foreach DIR,$(SRC_DIRS),-I$(DIR))

ifeq "$(NO_POST_REL_CHANGES_TEST)" "1"
CFLAGS   += -DNO_POST_REL_CHANGES
endif

ifeq "$(SUBSET_SQ)" "1"
CFLAGS += -DSUBSET_SQ
endif

ifeq "$(SUBSET_HQ)" "1"
CFLAGS += -DSUBSET_HQ
endif

ifeq "$(SUBSET_SWB)" "1"
CFLAGS += -DSUBSET_SWB
endif

ifeq "$(SUBSET_FB)" "1"
CFLAGS += -DSUBSET_FB
endif

ifeq "$(PLC)" "0"
CFLAGS   += -DDISABLE_PLC
endif

ifeq "$(PLC)" "2"
CFLAGS += -DDISABLE_ADVANCED_PLC
endif

# disable wmops instrumentation
ifeq "$(WMOPS)" "0"
    CPPFLAGS += -DWMOPS=0 -DDONT_COUNT_MEM
endif

# dependency magic
CC_FLAGS    = '$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS)'
POSTCOMPILE = mv -f $(BUILD)/$*.Td $(BUILD)/$*.d && touch $@

###############################################################################

SRCS := $(notdir $(foreach DIR, $(VPATH), $(wildcard $(DIR)/*.c)))
OBJS := $(addprefix $(BUILD)/, $(SRCS:.c=.o))

###############################################################################

.PHONY: all clean help force

all: $(NAME_LC3)

help:
	@echo 'Syntax: make [OPTION=VALUE ...]'
	@echo 'Build options:'
	@echo '    KISSFFT     $(KISSFFT) [0,1]'
	@echo '    NO_POST_REL $(NO_POST_REL) [0,1]'
	@echo '    OPTIM       $(OPTIM) [0-3]'
	@echo '    PLC         $(PLC) [0-2]'
	@echo '    SUBSET      $(SUBSET) [SQ,HQ,SWB,FB]'
	@echo '    WMOPS       $(WMOPS) [0,1]'
	@echo 'Debug options:'
	@echo '    AFL         $(AFL) [0,1]'
	@echo '    CLANG       $(CLANG) [0-3]'
	@echo '    GCOV        $(GCOV) [0,1]'

$(NAME_LC3): $(OBJS)
	@echo 'Linking' $@
	$(QUIET) $(LINK)  $(OBJS) -o $@ $(LDFLAGS)

clean:
	$(QUIET) rm -rf $(NAME_LC3) $(BUILD)

$(BUILD)/%.o : %.c $(BUILD)/cc_flags
	@echo 'Compiling' $<
	$(QUIET) $(CC) $(DEPFLAGS) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<
	$(QUIET) $(POSTCOMPILE)

# force rebuild if compilation flags changed
$(BUILD)/cc_flags: force
	$(QUIET) mkdir -p $(BUILD)
	$(QUIET) echo $(CC_FLAGS) | cmp -s - $@ || echo $(CC_FLAGS) > $@

# force rebuild if include dependency changed
$(BUILD)/%.d: ;
include $(wildcard $(patsubst %, $(BUILD)/%.d, $(basename $(SRCS))))
