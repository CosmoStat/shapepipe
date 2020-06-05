[Home](./shapepipe.md)

# Validation Checklist

- [x] PSF: per-tile ellipticity vs position plots (model, residual)  
- [x] PSF: compute diagnostics on star residuals. Make plots of residual PSF pattern.
- [ ] PSF: use above diagnostics to compute "Paulin" diagnostics (àla DES Y1 eqs (3.4-3.8)
- [ ] PSF: run diagnostics on 5% of test stars kept unused for the PSF modelling
- [ ] PSF: compare diagnostics for `PSF_SAMP` and `BASIS_TYPE` PSFEx parameters
- [ ] B-modes consistent with zero
- [x] Histogram of ellipticity components: 0 average, symmetric and same for both components
- [ ] Cross-correlation of PSF shape and other outputs (kappa maps, galaxy shapes)
- [ ] Angular correlation function of galaxies (to check photometric homogeneity between tiles, and masks)
- [x] Maps: Number counts map correlated with the data quality of the regions
- [x] Maps: Number counts map uncorrelated with: e1, e2, kE, kB, PSFellip, PSFsize (other maps to be added)
- [x] Maps: ellipticity maps uncorrelated with: num counts, the other e component, kE, kB, PSF ellip, PSF size (other maps to be added)
- [x] Maps: convergence (kE) maps uncorrelated with: num counts, ellipticity, kB, PSF ellip, PSF size (other maps to be added)
- [x] Maps: B-mode (kB) maps uncorrelated with: num counts, ellipticity, kE, PSF ellip, PSF size (other maps to be added)
- [x] Maps: PSF ellipticity maps uncorrelated with: num counts, ellipticity, kE, kB, the other ellip component, PSF size (other maps to be added)
- [x] Maps: PSF size maps uncorrelated with: num counts, ellipticity, kE, kB, PSF ellipticity (other maps to be added)
- [x] Calibrated ellipticities uncorrelated with: magnitude, flux, size
- [x] Calibrated ellipticities uncorrelated with flux PSF properties: ellipticity, size
