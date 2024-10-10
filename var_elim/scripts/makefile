structure-local: # Run structural results in parallel, locally, with GNU parallel
	python write_command_lines.py structure           # Writes commands to structure-commands.txt
	parallel -a structure-commands.txt                # Writes results to results/structure
	python collect_results.py structure               # Collects results into results/structure.csv
	python write_latex_table.py results/structure.csv # Writes table to results/structure.txt

# I don't think it makes sense to run solvetime results locally with multi-process
# parallelism. There will be too much noise. So I won't include a make recipe for it.

# Running all this at once might be a lot. It may be useful to just run each
# sweep individually.
sweep-local:
	python write_sweep_command_lines.py
	parallel -a commands/sweep-commands-distill-no-elim.txt
	parallel -a commands/sweep-commands-distill-d1.txt
	parallel -a commands/sweep-commands-distill-trivial.txt
	parallel -a commands/sweep-commands-distill-linear-d2.txt
	parallel -a commands/sweep-commands-distill-d2.txt
	parallel -a commands/sweep-commands-distill-ampl.txt
	parallel -a commands/sweep-commands-distill-matching.txt
	parallel -a commands/sweep-commands-mb-steady-no-elim.txt
	parallel -a commands/sweep-commands-mb-steady-d1.txt
	parallel -a commands/sweep-commands-mb-steady-trivial.txt
	parallel -a commands/sweep-commands-mb-steady-linear-d2.txt
	parallel -a commands/sweep-commands-mb-steady-d2.txt
	parallel -a commands/sweep-commands-mb-steady-ampl.txt
	parallel -a commands/sweep-commands-mb-steady-matching.txt
	parallel -a commands/sweep-commands-pipeline-no-elim.txt
	parallel -a commands/sweep-commands-pipeline-d1.txt
	parallel -a commands/sweep-commands-pipeline-trivial.txt
	parallel -a commands/sweep-commands-pipeline-linear-d2.txt
	parallel -a commands/sweep-commands-pipeline-d2.txt
	parallel -a commands/sweep-commands-pipeline-ampl.txt
	parallel -a commands/sweep-commands-pipeline-matching.txt
	parallel -a commands/collect-sweep-commands.txt # Collects results of each sweep into its own CSV file
	parallel -a plot-sweep-commands.txt
