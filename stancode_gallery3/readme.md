# Model description

## Notes
 * Code compiles and intial tests are underway. After testing is done, it is recommended to use gallery 3 code for the new data files.

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

## Changes compared to gallery 2

The main change only affects model with either size-dependent growth or respiration. While size limits were previously computed based on the absolute size difference between the size classes:
```
// compute size-dependent gamma and rho
if (xi > 0){
    sizelim = exp(xi*(v[i]-v[m]));
} else {
    sizelim = exp(xi*(v[i]-v[1]));
}
```
it is now based on the _relative_ difference:
```
// compute size-dependent gamma and rho
if (xi > 0){
    sizelim = exp(xi*(v[i]-v[m])/(v[m]-v[1]));
} else {
    sizelim = exp(xi*(v[i]-v[1])/(v[m]-v[1]));
}
```
This change makes the prior of `xi` (and `xir` for respiration) less dependent on changes in `v`.

**Note:** Size limits are also pre-computed more efficiently now.


## What all models have in common:
 * All models estimate the initial conditions the same way.
 * All models are currently using the same, updated multinomial to fit the data.
 * All models are using an updated priors compared to the code in [gallery 1](/stancode_gallery1) but the same as in [gallery 2](/stancode_gallery2).

## Model differences:

| version | core model <sup>[\[1\]](#corefootnote) | `delta_max` | using respiration | size-dep respiration | size-dep growth | light-dep division | using net growth <sup>[\[2\]](#netfootnote) | growth/respiration version <sup>[\[3\]](#versionfootnote) |
| ------- | ---------- | ----------  | --- | --- | --- | --- | --- | -------------------------- |
|`m1`     |            | free        |     |     |     |     |     | basic                      |
|`m2`     | ✓          | monotonic   |     |     |     |     |     | basic                      |
|`m3`     |            | monotonic   |     |     | ✓   |     |     | `gammaiv6`                 |
|`m4`     | ✓          | monotonic   | ✓   |     |     |     |     | `respv1`                   |
|`m5`     | ✓          | monotonic   | ✓   |     |     |     | ✓   | `respv2`                   |
|`m6`     |            | monotonic   | ✓   | ✓   | ✓   |     | ✓   | `respiv6`                  |
|`m7`     |            | monotonic   | ✓   | ✓   | ✓   |     |     | `respiv7`                  |
|`m8`     |            | monotonic   | ✓   | ✓   | ✓   | ✓   | ✓   | `respiv6`                  |
|`m9`     |            | monotonic   | ✓   | ✓   | ✓   | ✓   |     | `respiv7`                  |
|`m10`    |            | monotonic   | ✓   |     |     | ✓   | ✓   | `respv2`                   |
|`m11`    |            | monotonic   | ✓   |     | ✓   |     | ✓   | `resp_gammaiv6`            |
|`m12`    | ✓          | monotonic   | ✓   |     | ✓   |     |     | `resp_gammaiv7`            |
|`m13`    |            | monotonic   | ✓   |     | ✓   | ✓   | ✓   | `resp_gammaiv6`            |
|`m14`    | ✓          | monotonic   | ✓   |     | ✓   | ✓   |     | `resp_gammaiv7`            |

<a name="corefootnote">[1]</a> Based on discussion on 2020-09-01.

<a name="netfootnote">[2]</a> Growth and respiration cannot occur at the same time.

<a name="versionfootnote">[3]</a> See [this notebook](/sizedep_formulations.ipynb).

## Initial tests

 * Initial tests are underway.
 
## Proposed naming schemes

[see wiki](https://github.com/fribalet/Bayesian-matrixmodel/wiki/Model-names)




