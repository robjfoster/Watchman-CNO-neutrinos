# Set the oscillation parameters for solar neutrinos. Patrignani et al. 2016 "Review of Particle Physics"
t12 = np.radians(33)
t13 = np.radians(8)
# Integrated fluxes from N Vingoles
fi_n = 4.93e8
fi_o = 5.46e8
fi_f = 2.05e6


# Convert the fluxes as is appropriate for the chosen oscillation scenario. Only vacuum oscillations are considered, as the relevant energies are below 2 MeV.
# See "The MSW effect and solar neutrinos" - Smirnov [arXiv:hep-ph/0305106]
# For the survival probabilities see "Numerical Methods and Approximations for the Calculation of Solar Neutrino Fluxes in Three Flavours Focused on Results for the Sudbury Neutrino Observatory" 
# - Ryan Martin. Section 2.1.2 where we are here using the average over the oscillitory terms. We can justify this because the distance to the Sun and the volume of production in the Sun is large 
# compared to the wavelength of oscillation for these energies.
P = 1 - 0.5 * np.cos(t13) ** 4 * np.sin(2 * t12) ** 2 - 0.5 * np.sin(2 * t13) ** 2 
fbinned = fbinned * P