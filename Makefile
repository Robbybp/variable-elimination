DIR=var_elim/scripts
COMDIR=commands
RESDIR=results
IMDIR=images

structure-parallel:
	mkdir -p $(RESDIR)/structure # Letting this get created by the parallel scripts below leads to a race condition
	python $(DIR)/write_command_lines.py structure --results-dir=$(RESDIR) --commands-dir=$(COMDIR)
	parallel -a $(COMDIR)/structure-commands.txt
	python $(DIR)/collect_results.py structure --results-dir=$(RESDIR)
	python $(DIR)/write_latex_table.py $(RESDIR)/structure.csv --results-dir=$(RESDIR)
	python $(DIR)/write_latex_table.py $(RESDIR)/structure.csv --results-dir=$(RESDIR) --which=matching-bounds
	python $(DIR)/plot_structure_bargraphs.py $(RESDIR)/structure.csv --image-dir=$(IMDIR)
	python $(DIR)/plot_sparsity.py --model=mb-steady --image-dir=$(IMDIR)

solvetime-batch:
	mkdir -p $(RESDIR)/solvetime
	python $(DIR)/write_command_lines.py solvetime --results-dir=$(RESDIR)/solvetime --commands-dir=$(COMDIR)
	# This is a custom command I use to submit batch jobs on multiple HPC nodes
	submit-batch.sh $(COMDIR)/solvetime-commands.txt

solvetime-collect:
	python $(DIR)/collect_results.py solvetime --results-dir=$(RESDIR)
	python $(DIR)/write_latex_table.py $(RESDIR)/solvetime.csv --results-dir=$(RESDIR)
	python $(DIR)/plot_timing_bargraphs.py $(RESDIR)/solvetime.csv --image-dir=$(IMDIR)

sweep-batch:
	mkdir -p $(COMDIR)
	mkdir -p $(RESDIR)/sweep
	python $(DIR)/write_sweep_command_lines.py --results-dir=$(RESDIR)/sweep --commands-dir=$(COMDIR) --image-dir=$(IMDIR)
	# This is a custom command I use to submit batch jobs on multiple HPC nodes
	submit-batch.sh $(COMDIR)/parallel-sweep-commands.txt

sweep-collect:
	mkdir -p $(IMDIR)
	parallel -a $(COMDIR)/collect-sweep-commands.txt
	parallel -a $(COMDIR)/plot-sweep-commands.txt
	python $(DIR)/summarize_sweep_results.py --results-dir=$(RESDIR)/sweep --model=distill
	python $(DIR)/summarize_sweep_results.py --results-dir=$(RESDIR)/sweep --model=mb-steady
	python $(DIR)/summarize_sweep_results.py --results-dir=$(RESDIR)/sweep --model=pipeline
	python $(DIR)/summarize_sweep_results.py --results-dir=$(RESDIR)/sweep
	python $(DIR)/write_latex_table.py $(RESDIR)/sweep-summary.csv --results-dir=$(RESDIR)
