# MFWS_RapidFloristics_2021

Start with "TubeMFWS_SummaryFigures.R" which imports "TuberousPlants_MFWS_Data.xslx" to plot the raw data as spatial maps (lat/long location of points coloured by area and sized by response variable. The data file also includes a metadata tab that explains each variable in the dataframe. Any questions on methodology etc. please get in contact.

Second is "TubeMFWS_2_FirstDataAnalysis" which imports "MFWS-RapidFloristicData-2021.csv" (same observations as previous data file but with additional mapped explanatory variables joined), "MFWS-VegStructData-2021.csv" (the matched vegetation structure monitoring that occurred at most (88/112) plots), and "MULLonly-RapidFloristicData-2021.csv" (same observations as previous data file but only > Area == "MULL" observations and a few new columns that represent the experimental treatments within MULL in a few different ways

SEPTEMBER 2022 UPDATE
All analyses have been re-run for the manuscript that is in preparation. 
All summary and analyses results are now reproducable from the "TubeMFWS_4_FINALANALYSIS.R" code. No need to revist older code unless wanting to reproduce the preliminary results as presented and shared in early-2022
These results have been / will be publicly presented at:
- "Ecofocus" at Gold Creek Station, 15-09-2022
- "Science at the Sanctuary" at Wildbark, 06-10-2022
- ESA-SCOB 2022 Conference at Wollongong, 28-11-2022
