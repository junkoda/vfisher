vfisher
=======

Fisher matrix forecast code for peculiar velocity surveys (cosmology)

## Compile

```
  % cd src
  % make
```
GSL and Boost libraries are required. See `src/Makefile` if you have them in non-standard paths.

## Run

Check available options with,

```
  % vfisher --help
```

Typical run is:

```
  % src/vfihser --noise=uniform_ngal --ngal=1.0e-2 \
  --steradian=4pi --verr=0.2 --pk-dir=../data/power_spectra \
  --planck-prior=../data/planck_cov.txt --b=1.0 --sigma_u=13.4 \
  --pk=rpt --zmax=0.1 --kmax=0.2 k
```

- uniform galaxy density is assumed `--noise=uniforme_ngal`
- galaxy number density is 1.0e-2 (1/h Mpc)^-3 ``--ngal=1.0e-2`
- field of view is full sky ``--steradian=4pi`
- galaxy linear bias `--b=1.0`
- `k` mode prints the result as a function of kmax,
- an alternative is `z` mode, as a function of zmax,
- or `final` that prints only at --kmax and --zmax

## Output

See header of the output.
In `k` mode, column 1 is kmax,

```
# 1 k
```

Following means that free parameters fitted against data would be beta, fsigma8, Omega[b] h^2, Omega[CDM] h^2, h, and ns, and columns 26-31 show the fractional constraints on these parameters with the information k < kmax. Note that Planck prior is used on the latter 4 parameters with the option --planck-prior=../data/planck_cov.txt.

```
# --free: beta fsigma8 omegabh2 omegach2 h ns
# 26: dbeta/beta
# 27: dfsigma8/fsigma8
# 28: domegabh2/omegabh2
# 29: domegach2/omegach2
# 30: dh/h
# 31: dns/ns
```

# References

Please cite Koda et al (2014), MNRAS, **445**, 4267 [[1]], if you use this code.

If you would like the read the code, I also recommend looking at a simpler code by Martin White [[2]], [[3]]

1. Koda et al. (2014),
*Are peculiar velocity surveys competitive as a cosmological probe?* MNRAS, **445**, 4267

2. 	White, Song, and Percival (2009)
*Forecasting cosmological constraints from redshift surveys*
MNRAS, **397**, 1348

[1]: http://adsabs.harvard.edu/abs/2014MNRAS.445.4267K
[2]: http://mwhite.berkeley.edu/Redshift/
[3]: http://adsabs.harvard.edu/abs/2009MNRAS.397.1348W







  