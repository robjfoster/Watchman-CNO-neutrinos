# watchman-CNO-neutrinos

Analysis of watchman sensitivity to CNO solar neutrinos

compare_and_plot_spectra.py: Compare the data for solar neutrino spectra for both the raw data and the rat-pac formatted data by plotting both. All data is found in the data file, the rat-pac formatted data has been copied into the script rather than read in.

comparespectra.py: Compares the CNO spectrum created to the one created by a summer student.

plotsolarspectra.py: Plots the raw data for the components of the solar neutrino spectrum.

spectra4ratpacX.py: Converts the raw data for component X to a format that can be read for the rat-pac simulation software. All data can be found in the data directory, with the .ratdb files being written there.

CNOsolarNeutrinoBins.py: Compares how the number of bins used in creating CNO spectrum from its components effected the flux. This was used to solve an issue in normalising the spectrum.

activities.py, rates.py, significance.py, masses.py all taken from a summer studemt's code. These contain funtionality from watchmakers.
