# Model description

## Notes:
 1. All code has been compiled but not all the models have been tested yet.
 2. 2 models with light-dependent have been added.

## The model versions contained in gallery1:
`m1`:`matrixmodel_multinom_estinilnorm_freedelta_normparam_trackgrowth_xval.stan`
`m2`:`matrixmodel_multinom_estinilnorm_monodelta_normparam_trackgrowth_xval.stan`
`m3`:`matrixmodel_multinom_estinilnorm_monodelta_gammaiv6_normparam_trackgrowth_xval.stan`
`m4`:`matrixmodel_multinom_estinilnorm_monodelta_respv1_normparam_trackgrowth_xval.stan`
`m5`:`matrixmodel_multinom_estinilnorm_monodelta_respv2_normparam_trackgrowth_xval.stan`
`m6`:`matrixmodel_multinom_estinilnorm_monodelta_respiv6_normparam_trackgrowth_xval.stan`
`m7`:`matrixmodel_multinom_estinilnorm_monodelta_respiv7_normparam_trackgrowth_xval.stan`
`m8`:`matrixmodel_multinom_estinilnorm_monodelta-lightsig_respiv6_normparam_trackgrowth_xval.stan`
`m9`:`matrixmodel_multinom_estinilnorm_monodelta-lightsig_respiv7_normparam_trackgrowth_xval.stan`

## What all models have in common:
 * All models estimate the initial conditions the same way.
 * All models are currently using a simple multinomial to fit the data.
 * All models are using an updated `uniform(0.0,10.0)` prior for `gamma_max` and the models with respiration use a `uniform(0.0, 3.0)` for `rho_max`.

## Model differences:

| version | `delta_max` | respiration  | size-dep growth | light-dep division | using net growth <sup>[\[1\]](#netfootnote) | growth/respiration version <sup>[\[2\]](#versionfootnote) |
| ------- | ----------  | --- | --- | --- | --- | -------------------------- |
|`m1`     | free        |     |     |     |     | basic                      |
|`m2`     | monotonic   |     |     |     |     | basic                      |
|`m3`     | monotonic   |     | ✓   |     |     | `gammaiv6`                 |
|`m4`     | monotonic   | ✓   |     |     |     | `respv1`                   |
|`m5`     | monotonic   | ✓   |     |     | ✓   | `respv2`                   |
|`m6`     | monotonic   | ✓   | ✓   |     | ✓   | `respiv6`                  |
|`m7`     | monotonic   | ✓   | ✓   |     |     | `respiv7`                  |
|`m8`     | monotonic   | ✓   | ✓   | ✓   | ✓   | `respiv6`                  |
|`m9`     | monotonic   | ✓   | ✓   | ✓   |     | `respiv7`                  |

<a name="netfootnote">[1]</a> Growth and respiration cannot occur at the same time.

<a name="versionfootnote">[2]</a> See [this notebook](/sizedep_formulations.ipynb).

## Initial tests

 * Initial results of models `m1` to `m7` are presented in [this notebook](/experimental/experimental_zinser_seaflow_20200602_gallery1_test.ipynb)
