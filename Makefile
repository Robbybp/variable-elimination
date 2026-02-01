# TODO: Variable for results dir and/or suffix
DIR=var_elim/scripts
COMDIR=commands
RESDIR=results
IMDIR=images

all:
	echo "WARNING: Reproducing all results in serial via `make all` is time-consuming."
	echo "We recommend reproducing these results in an HPC environment with a resource"
	echo "manager like Slurm, if available."
	# Generate results for Table 3, Table 4, Figure 5, Figure 6, and Figure 7
	python $(DIR)/analyze_structure.py
	# Table 3
	python $(DIR)/write_latex_table.py results/structure.csv
	# Table 4
	python $(DIR)/write_latex_table.py results/structure.csv --which=matching-bounds
	# Figures 5 and 7
	python $(DIR)/plot_structure_bargraphs.py results/structure.csv
	# TODO: Figure 6
	#
	# Generate results for Figure 8 and Table 5
	python $(DIR)/analyze_solvetime.py
	# Table 5
	python $(DIR)/write_latex_table.py results/solvetime.csv
	# Figure 8
	python $(DIR)/plot_timing_bargraphs.py results/solvetime.csv
	# Generate results for Table 6 and Figures 9, 10, 11, and 12
	python $(DIR)/run_param_sweep.py --model=distill
	python $(DIR)/run_param_sweep.py --model=mb-steady
	python $(DIR)/run_param_sweep.py --model=pipeline
	# Table 6
	python $(DIR)/summarize_sweep_results.py --model=distill
	python $(DIR)/summarize_sweep_results.py --model=mb-steady
	python $(DIR)/summarize_sweep_results.py --model=pipeline
	#
	# Each parameter sweep plot is produced individually...
	# Figure 9
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/distill-no-elim-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/distill-d1-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/distill-trivial-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/distill-linear-d2-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/distill-d2-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/distill-ampl-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/distill-matching-sweep.csv
	# Figure 10
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/mb-steady-no-elim-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/mb-steady-d1-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/mb-steady-trivial-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/mb-steady-linear-d2-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/mb-steady-d2-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/mb-steady-ampl-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/mb-steady-matching-sweep.csv
	# Figure 11
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/pipeline-no-elim-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/pipeline-d1-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/pipeline-trivial-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/pipeline-linear-d2-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/pipeline-d2-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/pipeline-ampl-sweep.csv
	python $(DIR)/plot_sweep_results.py $(DIR)/results/sweep/pipeline-matching-sweep.csv

structure-parallel:
	mkdir -p $(RESDIR)/structure # Letting this get created by the parallel scripts below leads to a race condition
	python $(DIR)/write_command_lines.py structure --results-dir=$(RESDIR) --commands-dir=$(COMDIR)
	parallel -a $(COMDIR)/structure-commands.txt
	python $(DIR)/collect_results.py structure --results-dir=$(RESDIR)
	python $(DIR)/write_latex_table.py $(RESDIR)/structure.csv --results-dir=$(RESDIR)
	python $(DIR)/write_latex_table.py $(RESDIR)/structure.csv --results-dir=$(RESDIR) --which=matching-bounds
	python $(DIR)/plot_structure_bargraphs.py $(RESDIR)/structure.csv --image-dir=$(IMDIR)

solvetime-batch:
	mkdir -p $(RESDIR)/solvetime
	python $(DIR)/write_command_lines.py solvetime --results-dir=$(RESDIR)/solvetime --commands-dir=$(COMDIR)
	# This is a custom command I use to submit batch jobs on multiple HPC nodes
	submit-batch.sh $(COMDIR)/solvetime-commands.txt

solvetime-collect:
	python $(DIR)/collect_results.py solvetime --results-dir=$(RESDIR)
	python $(DIR)/write_latex_table.py $(RESDIR)/solvetime.csv --results-dir=$(RESDIR)
	python $(DIR)/plot_timing_bargraphs.py $(RESDIR)/solvetime.csv --image-dir=$(IMDIR)
