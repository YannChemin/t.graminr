MODULE_TOPDIR = ../..

PGM = t.graminr

include $(MODULE_TOPDIR)/include/Make/Script.make
include $(MODULE_TOPDIR)/include/Make/Html.make

default: script html etcfiles

# t.graminr is a thin GRASS front end over a GRASS-free, numpy-vectorised
# model package. Both the package and the default parameter sets are whole
# directories rather than single files, so ETCFILES cannot carry them;
# they are copied into $(ETC)/$(PGM)/, which Script.make's own install rule
# then takes to $(INST_DIR)/etc/$(PGM) as a unit. t.graminr.py finds them
# there relative to its own installed location.
ETCDIR = $(ETC)/$(PGM)

etcfiles: $(ETCDIR)/.stamp

$(ETCDIR)/.stamp: $(wildcard graminr_lib/*.py) $(wildcard params/*.csv)
	@mkdir -p $(ETCDIR)/graminr_lib $(ETCDIR)/params
	$(INSTALL_DATA) graminr_lib/*.py $(ETCDIR)/graminr_lib/
	$(INSTALL_DATA) params/*.csv $(ETCDIR)/params/
	touch $@

.PHONY: etcfiles
