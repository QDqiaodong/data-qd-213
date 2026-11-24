# Spectral wavelength calibration contract

Laboratory data processors supply 3..1000 identified reference lines ordered by strictly increasing detector pixel and wavelength, plus an increasing list of pixels to calibrate. Pixel and wavelength values are bounded decimal strings and are handled as exact rational numbers. The service fits wavelength = c0 + c1*pixel + c2*pixel² by exact least squares using the normal equations.

Every reference retains its signed residual. The full calibration is accepted only when the maximum absolute residual does not exceed max_residual_nm. Query pixels must remain inside the reference range; the service refuses extrapolation. It does not discover peaks, identify chemical species, correct instrument drift or estimate statistical uncertainty.

The output is a calibration data product for downstream spectra, not a visualization page or generic regression utility. Future exposure records can reference the coefficient set and residual evidence through a versioned calibration identity.

