# PDC_Analysis_2023
Analysis work on Sherbrooke PDC from MIEL measurements

The PDC_DarkAnalysis.C Macro takes in 5 data runs (I used a VOV range of 1-5 for higher temperature, 
and 4-8 for lower temperature) from root_output_files and creates time difference distributions, 
fits, and graphs to show the dark count rate vs VOV. 

Specific fit parameters are given for each temperature based on how the distribution behaves (fit ranges are done by eye). 

Run PDC_DarkAnalysis.C from the root command line (requires a root install on your computer!) with

root -l

.x PDC_DarkAnalysis.C(run1,run2,run3,run4,run5,anode,temp)

where anode specifies the PDC anode and temp is the temperature (both ints). 

The PDC_HoldOff.C Macro runs on the same body as the dark analysis, but is less developed. 
It was used to do a specific study of hold off (digital pulse width) vs dark count rate.

TGraph_Stacker_PDC.C takes any number of outputs from PDC_DarkAnalysis and creates a stacked graph of the outputs.
It is currently written to make plots that show each run with a legend for each temp, and a projection along 4VoV to show the temperature dependence of dark count rate.

bulkProcessing.py is useful if you want to do a bulk run of the DarkAnalysis. Edit it with the run numbers you want to run and the temps and anodes that correspond 
and it will automate the macro for any number of iterations.
