GMX ?= ~/projects/gromacs/build/bin/gmx
GMXFLAGS ?= 
SOLVATE_FLAGS ?= -cs spc216.gro
MINIMIZATION_MDRUN_FLAGS ?= -nb cpu -v
EQUILIBRATION_MDRUN_FLAGS ?= -nb gpu -bonded cpu -update gpu -v
PRODUCTION_MDRUN_FLAGS ?= -nb gpu -bonded cpu -update gpu -v
ENERGY_TERM ?= Total-Energy

X = 100
Y = 3
Z = 3

TOP = topol.top
TEMPLATE_TOP = topol.template.top
CONF = conf.gro

EM_MDP = em.mdp
EM_DEFFNM = em
EM_TPR = $(EM_DEFFNM).tpr
EM_CONF = $(EM_DEFFNM).gro

NVT_MDP = nvt.mdp
NVT_DEFFNM = nvt
NVT_TPR = $(NVT_DEFFNM).tpr
NVT_CONF = $(NVT_DEFFNM).gro

MDP = md.mdp
DEFFNM = md
TPR = topol.tpr
PBC_TRAJ = $(DEFFNM)_pbc.xtc
ENERGY_XVG = total_energy.xvg
ENERGY_TXT = total_energy.dat
DRIFT_TXT = energy_drift.txt

.PHONY: all info solvate minimize equilibrate prepare run pbc energy clean distclean

all: energy

info:
	@echo "Build = candidate"
	@echo "GMX = $(GMX)"
	@echo "Box = $(X) $(Y) $(Z) nm"
	@echo "PRODUCTION_MDRUN_FLAGS = $(PRODUCTION_MDRUN_FLAGS)"

solvate:
	cp $(TEMPLATE_TOP) $(TOP)
	$(GMX) solvate $(GMXFLAGS) $(SOLVATE_FLAGS) \
		-box $(X) $(Y) $(Z) \
		-o $(CONF) \
		-p $(TOP)

minimize: solvate
	$(GMX) grompp $(GMXFLAGS) \
		-f $(EM_MDP) \
		-c $(CONF) \
		-p $(TOP) \
		-o $(EM_TPR)
	$(GMX) mdrun $(GMXFLAGS) \
		-s $(EM_TPR) \
		-deffnm $(EM_DEFFNM) \
		$(MINIMIZATION_MDRUN_FLAGS)

equilibrate: minimize
	$(GMX) grompp $(GMXFLAGS) \
		-f $(NVT_MDP) \
		-c $(EM_CONF) \
		-p $(TOP) \
		-o $(NVT_TPR)
	$(GMX) mdrun $(GMXFLAGS) \
		-s $(NVT_TPR) \
		-deffnm $(NVT_DEFFNM) \
		$(EQUILIBRATION_MDRUN_FLAGS)

prepare: equilibrate
	$(GMX) grompp $(GMXFLAGS) \
		-f $(MDP) \
		-c $(NVT_CONF) \
		-p $(TOP) \
		-o $(TPR)

run: prepare
	$(GMX) mdrun $(GMXFLAGS) \
		-s $(TPR) \
		-deffnm $(DEFFNM) \
		$(PRODUCTION_MDRUN_FLAGS)

pbc: run
	printf "System\nSystem\n" | $(GMX) trjconv $(GMXFLAGS) \
		-s $(TPR) \
		-f $(DEFFNM).trr \
		-o $(PBC_TRAJ) \
		-pbc mol \
		-ur compact \
		-center

energy: pbc
	printf "$(ENERGY_TERM)\n0\n" | $(GMX) energy $(GMXFLAGS) \
		-f $(DEFFNM).edr \
		-o $(ENERGY_XVG) \
		-xvg none
	awk 'BEGIN { first_set=0 } /^[@#]/ { next } NF >= 2 { if (!first_set) { first_t=$$1; first_e=$$2; first_set=1 } last_t=$$1; last_e=$$2 } END { if (!first_set) { print "No energy samples found" > "/dev/stderr"; exit 1 } drift=last_e-first_e; duration=last_t-first_t; printf "first_time_ps %.8g\nfirst_energy_kj_mol %.12g\nlast_time_ps %.8g\nlast_energy_kj_mol %.12g\ndrift_kj_mol %.12g\ndrift_per_ps_kj_mol %.12g\n", first_t, first_e, last_t, last_e, drift, (duration == 0 ? 0 : drift / duration) }' $(ENERGY_XVG) > $(DRIFT_TXT)
	cp $(ENERGY_XVG) $(ENERGY_TXT)

clean:
	rm -f $(CONF) $(TOP)
	rm -f $(EM_TPR) $(EM_CONF) em.log em.edr em.trr em.xtc em.cpt em_prev.cpt
	rm -f $(NVT_TPR) $(NVT_CONF) nvt.log nvt.edr nvt.trr nvt.xtc nvt.cpt nvt_prev.cpt
	rm -f $(TPR)
	rm -f md.log md.edr md.trr md.xtc md.gro md.cpt md_prev.cpt $(PBC_TRAJ)
	rm -f $(ENERGY_XVG) $(ENERGY_TXT) $(DRIFT_TXT)

distclean: clean
	rm -f \#*
