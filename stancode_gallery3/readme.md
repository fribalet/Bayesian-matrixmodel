# Model description

## Notes
 * Code compiles and intial tests are *again* underway. After testing is done, it is recommended to use gallery 3 code for the new data files.

## The model versions contained in gallery 3:
`m1`:`matrixmodel_mlmultinom_estinilnorm2_freedelta_normparam_trackgrowthvol_xval2.stan`
`m2`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_normparam_trackgrowthvol_xval2.stan`
`m3`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_gammaiv6_normparam_trackgrowthvol_xval2.stan`
`m4`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_respv1_normparam_trackgrowthvol_xval2.stan`
`m5`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_respv2_normparam_trackgrowthvol_xval2.stan`
`m6`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_respiv6_normparam_trackgrowthvol_xval2.stan`
`m7`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_respiv7_normparam_trackgrowthvol_xval2.stan`
`m8`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_respiv6_normparam_trackgrowthvol_xval2.stan`
`m9`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_respiv7_normparam_trackgrowthvol_xval2.stan`
`m10`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_respv2_normparam_trackgrowthvol_xval2.stan`
`m11`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_resp_gammaiv6_normparam_trackgrowthvol_xval2.stan`
`m12`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_resp_gammaiv7_normparam_trackgrowthvol_xval2.stan`
`m13`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_resp_gammaiv6_normparam_trackgrowthvol_xval2.stan`
`m14`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_resp_gammaiv7_normparam_trackgrowthvol_xval2.stan`
`m3u`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_gammaiv8_normparam_trackgrowthvol_xval2.stan`
`m6u`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_respiv8_normparam_trackgrowthvol_xval2.stan`
`m7u`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_respiv9_normparam_trackgrowthvol_xval2.stan`
`m8u`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_respiv8_normparam_trackgrowthvol_xval2.stan`
`m9u`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_respiv9_normparam_trackgrowthvol_xval2.stan`
`m11u`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_resp_gammaiv8_normparam_trackgrowthvol_xval2.stan`
`m12u`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_resp_gammaiv9_normparam_trackgrowthvol_xval2.stan`
`m13u`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_resp_gammaiv8_normparam_trackgrowthvol_xval2.stan`
`m14u`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_resp_gammaiv9_normparam_trackgrowthvol_xval2.stan`
`m15`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_resp_freegamma_normparam_trackgrowthvol_xval2.stan`
`m16`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_freeresp_freegamma_normparam_trackgrowthvol_xval2.stan`
`m4s6`:`../g3_timedep/matrixmodel_mlmultinom_estinilnorm2_monodelta2-timespline_respv1_normparam_trackgrowthvol_xval2.stan` 
`m12us6`:`../g3_timedep/matrixmodel_mlmultinom_estinilnorm2_monodelta2-timespline_resp_gammaiv9_normparam_trackgrowthvol_xval2.stan` 
`m15s6`:`../g3_timedep/matrixmodel_mlmultinom_estinilnorm2_monodelta2-timespline_resp_freegamma_normparam_trackgrowthvol_xval2.stan` 
`m16s6`:`../g3_timedep/matrixmodel_mlmultinom_estinilnorm2_monodelta2-timespline_freeresp_freegamma_normparam_trackgrowthvol_xval2.stan`

## Changes compared to gallery 2

The main change only affects model with either size-dependent growth or respiration. Gallery 2 only contained exponential size relationships, the updated code in gallery 3 (denoted by a `u` in the name) contains power law relationships.
```
if (exponent > 0){
    sizelim = (v_mid[i]^exponent)/(v_mid[m]^exponent);
} else {
    sizelim = (v_mid[i]^exponent)/(v_mid[1]^exponent);
}

```
The old exponential relationships have also been updated: while size limits were previously computed based on the absolute size difference between the size classes:
```
// compute size-dependent gamma and rho
if (xi > 0){
    sizelim = exp(xi*(v_mid[i]-v_mid[m]));
} else {
    sizelim = exp(xi*(v_mid[i]-v_mid[1]));
}
```
it is now based on the _relative_ difference:
```
// compute size-dependent gamma and rho
if (xi > 0){
    sizelim = exp(xi*(v_mid[i]-v_mid[m])/(v_mid[m]-v_mid[1]));
} else {
    sizelim = exp(xi*(v_mid[i]-v_mid[1])/(v_mid[m]-v_mid[1]));
}
```
This change makes the prior of `xi` (and `xir` for respiration) less dependent on changes in `v`. Note that the updated code now also uses `v_mid` for size-dependence and not `v`. `v_mid` is defined as:
```
for (i in 1:m){
    // using geometric mean
    v_mid[i] = sqrt(v[i] * v[i+1]);
}
```

**Note:** Size limits are also pre-computed more efficiently now.


## What all models have in common:
 * All models estimate the initial conditions the same way.
 * All models are currently using the same, updated multinomial to fit the data.
 * All models are using an updated priors compared to the code in [gallery 1](/stancode_gallery1) but the same as in [gallery 2](/stancode_gallery2).

## Model differences:

| name     | version | core model | old core model <sup>[\[1\]](#corefootnote) | `delta_max` | using respiration | size-dep respiration | size-dep growth | light-dep division | using net growth <sup>[\[2\]](#netfootnote) | using time-dep division | growth/respiration version <sup>[\[3\]](#versionfootnote) |
| -------  | ------- | --- | --- | ----------  | --- | --- | --- | --- | --- | --- | -------------------------- |
| `m_bfx`  |`m1`     |     |     | free        |     |     |     |     |     |     | basic                      |
| `m_bmx`  |`m2`     | ✓   | ✓   | monotonic   |     |     |     |     |     |     | basic                      |
| `m_emx`  |`m3`     |     |     | monotonic   |     |     | ✓   |     |     |     | `gammaiv6`                 |
| `m_bmb`  |`m4`     | ✓   | ✓   | monotonic   | ✓   |     |     |     |     |     | `respv1`                   |
| `m_bmb-` |`m5`     |     | ✓   | monotonic   | ✓   |     |     |     | ✓   |     | `respv2`                   |
| `m_eme-` |`m6`     |     |     | monotonic   | ✓   | ✓   | ✓   |     | ✓   |     | `respiv6`                  |
| `m_eme`  |`m7`     |     |     | monotonic   | ✓   | ✓   | ✓   |     |     |     | `respiv7`                  |
| `m_ele-` |`m8`     |     |     | monotonic   | ✓   | ✓   | ✓   | ✓   | ✓   |     | `respiv6`                  |
| `m_ele`  |`m9`     |     |     | monotonic   | ✓   | ✓   | ✓   | ✓   |     |     | `respiv7`                  |
| `m_blb`  |`m10`    |     |     | monotonic   | ✓   |     |     | ✓   | ✓   |     | `respv2`                   |
| `m_emb-` |`m11`    |     |     | monotonic   | ✓   |     | ✓   |     | ✓   |     | `resp_gammaiv6`            |
| `m_emb`  |`m12`    |     | ✓   | monotonic   | ✓   |     | ✓   |     |     |     | `resp_gammaiv7`            |
| `m_elb-` |`m13`    |     |     | monotonic   | ✓   |     | ✓   | ✓   | ✓   |     | `resp_gammaiv6`            |
| `m_elb`  |`m14`    |     | ✓   | monotonic   | ✓   |     | ✓   | ✓   |     |     | `resp_gammaiv7`            |
| `m_pmx`  |`m3u`    |     |     | monotonic   |     |     | ✓   |     |     |     | `gammaiv8`                 |
| `m_pmp-` |`m6u`    |     |     | monotonic   | ✓   | ✓   | ✓   |     | ✓   |     | `respiv8`                  |
| `m_pmp`  |`m7u`    |     |     | monotonic   | ✓   | ✓   | ✓   |     |     |     | `respiv9`                  |
| `m_plp-` |`m8u`    |     |     | monotonic   | ✓   | ✓   | ✓   | ✓   | ✓   |     | `respiv8`                  |
| `m_plp`  |`m9u`    |     |     | monotonic   | ✓   | ✓   | ✓   | ✓   |     |     | `respiv9`                  |
| `m_pmb-` |`m11u`   |     |     | monotonic   | ✓   |     | ✓   |     | ✓   |     | `resp_gammaiv8`            |
| `m_pmb`  |`m12u`   | ✓   |     | monotonic   | ✓   |     | ✓   |     |     |     | `resp_gammaiv9`            |
| `m_plb-` |`m13u`   |     |     | monotonic   | ✓   |     | ✓   | ✓   | ✓   |     | `resp_gammaiv8`            |
| `m_plb`  |`m14u`   |     |     | monotonic   | ✓   |     | ✓   | ✓   |     |     | `resp_gammaiv9`            |
| `m_fmb`  |`m15`    | ✓   |     | monotonic   | ✓   |     | ✓   |     |     |     | free growth                |
| `m_fmf`  |`m16`    | ✓   |     | monotonic   | ✓   | ✓   | ✓   |     |     |     | free growth, free division |
| `m_btb`  |`m4s6`   | ✓   |     | monotonic   | ✓   |     |     |     |     | ✓   | `respv1`                   |
| `m_ptb`  |`m12us6` | ✓   |     | monotonic   | ✓   |     | ✓   |     |     | ✓   | `resp_gammaiv9`            |
| `m_ftb`  |`m15s6`  | ✓   |     | monotonic   | ✓   |     | ✓   |     |     | ✓   | free growth                |
| `m_ftf`  |`m16s6`  | ✓   |     | monotonic   | ✓   | ✓   | ✓   |     |     | ✓   | free growth, free division |

<a name="corefootnote">[1]</a> Based on discussion on 2020-09-01.

<a name="netfootnote">[2]</a> Growth and respiration cannot occur at the same time.

<a name="versionfootnote">[3]</a> See [this notebook](/sizedep_formulations.ipynb).

## Initial tests

 * Initial tests are underway.
 
## Proposed naming schemes

[see wiki](https://github.com/fribalet/Bayesian-matrixmodel/wiki/Model-names)




