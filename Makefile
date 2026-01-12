# TODO: Variable for results dir and/or suffix

all:
	echo "WARNING: Reproducing all results in serial via `make all` is time-consuming."
	echo "We recommend reproducing these results in an HPC environment with a resource"
	echo "manager like Slurm, if available."
	# Generate results for Table 3, Table 4, Figure 5, Figure 6, and Figure 7
	python var_elim/scripts/analyze_structure.py
	# Table 3
	python var_elim/scripts/write_latex_table.py results/structure.csv
	# Table 4
	python var_elim/scripts/write_latex_table.py results/structure.csv --which=matching-bounds
	# Figures 5 and 7
	python var_elim/scripts/plot_structure_bargraphs.py results/structure.csv
	# TODO: Figure 6
	#
	# Generate results for Figure 8 and Table 5
	python var_elim/scripts/analyze_solvetime.py
	# Table 5
	python var_elim/scripts/write_latex_table.py results/solvetime.csv
	# Figure 8
	python var_elim/scripts/plot_timing_bargraphs.py results/solvetime.csv
	# Generate results for Table 6 and Figures 9, 10, 11, and 12
	python var_elim/scripts/run_param_sweep.py --model=distill
	python var_elim/scripts/run_param_sweep.py --model=mb-steady
	python var_elim/scripts/run_param_sweep.py --model=pipeline
	# Table 6
	python var_elim/scripts/summarize_sweep_results.py --model=distill
	python var_elim/scripts/summarize_sweep_results.py --model=mb-steady
	python var_elim/scripts/summarize_sweep_results.py --model=pipeline
	#
	# Each parameter sweep plot is produced individually...
	# Figure 9
	python plot_sweep_results.py var_elim/scripts/results/sweep/distill-no-elim-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/distill-d1-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/distill-trivial-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/distill-linear-d2-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/distill-d2-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/distill-ampl-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/distill-matching-sweep.csv
	# Figure 10
	python plot_sweep_results.py var_elim/scripts/results/sweep/mb-steady-no-elim-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/mb-steady-d1-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/mb-steady-trivial-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/mb-steady-linear-d2-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/mb-steady-d2-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/mb-steady-ampl-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/mb-steady-matching-sweep.csv
	# Figure 11
	python plot_sweep_results.py var_elim/scripts/results/sweep/pipeline-no-elim-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/pipeline-d1-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/pipeline-trivial-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/pipeline-linear-d2-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/pipeline-d2-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/pipeline-ampl-sweep.csv
	python plot_sweep_results.py var_elim/scripts/results/sweep/pipeline-matching-sweep.csv

structure-parallel:
	python var_elim/scripts/write_command_lines.py structure
	parallel -a structure-commands.txt
	python var_elim/scripts/collect_results.py structure
	python var_elim/scripts/write_latex_table.py results/structure.csv
	python var_elim/scripts/write_latex_table.py results/structure.csv --which=matching-bounds
	python var_elim/scripts/plot_structure_bargraphs results/structure.csv
