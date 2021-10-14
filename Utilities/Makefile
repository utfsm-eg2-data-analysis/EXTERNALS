.DELETE_ON_ERROR:

ifndef ROOTSYS
    $(error "ERROR: please, set the variable ROOTSYS")
endif

ifndef WORKDIR
    $(error "ERROR: please, set the variable WORKDIR")
endif

MKDIR_P := mkdir -p

BINDIR := ./bin
SRCDIR := ./src

ROOTCONFIG := root-config
ROOTCFLAGS := $(shell $(ROOTCONFIG) --cxx --cflags)
ROOTGLIBS := $(shell $(ROOTCONFIG) --glibs)

INCLUDES := -I./include

PROG := GetBinning GetCentroids
LIST := $(addprefix $(BINDIR)/, $(PROG))

.PHONY: all clean

all: $(LIST)

$(BINDIR)/%: $(SRCDIR)/%.cxx
	@echo "Compiling" $<
	$(MKDIR_P) $(BINDIR)
	$(ROOTCFLAGS) $(INCLUDES) -o $@ $< $(ROOTGLIBS)

clean:
	rm -f $(LIST:%=%.o) $(LIST)
