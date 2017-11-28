echo ">>>> RUN SIMULATION <<<<"
model_distances_per_transcipt.pl \
	> simulation_table.tab

do_plots.R \
	-i simulation_table.tab \
	-f simulation_table.pdf
 
