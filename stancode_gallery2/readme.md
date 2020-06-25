# Model description

## Notes
 * Code compiles but is otherwise untested, **do not use yet**.

## The model versions contained in gallery1:
`m1`:`matrixmodel_mlmultinom_estinilnorm2_freedelta_normparam_trackgrowth_xval2.stan`
`m2`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_normparam_trackgrowth_xval2.stan`
`m3`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_gammaiv6_normparam_trackgrowth_xval2.stan`
`m4`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_respv1_normparam_trackgrowth_xval2.stan`
`m5`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_respv2_normparam_trackgrowth_xval2.stan`
`m6`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_respiv6_normparam_trackgrowth_xval2.stan`
`m7`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_respiv7_normparam_trackgrowth_xval2.stan`
`m8`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_respiv6_normparam_trackgrowth_xval2.stan`
`m9`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_respiv7_normparam_trackgrowth_xval2.stan`
`m10`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_respv2_normparam_trackgrowth_xval2.stan`

## What all models have in common:
 * All models estimate the initial conditions the same way.
 * All models are currently using the same, updated multinomial to fit the data.
 * All models are using an updated priors compared to the code in [gallery1](/stancode_gallery1).

## Model differences:

| version | `delta_max` | using respiration | size-dep respiration | size-dep growth | light-dep division | using net growth <sup>[\[1\]](#netfootnote) | growth/respiration version <sup>[\[2\]](#versionfootnote) |
| ------- | ----------  | --- | --- | --- | --- | --- | -------------------------- |
|`m1`     | free        |     |     |     |     |     | basic                      |
|`m2`     | monotonic   |     |     |     |     |     | basic                      |
|`m3`     | monotonic   |     |     | ✓   |     |     | `gammaiv6`                 |
|`m4`     | monotonic   | ✓   |     |     |     |     | `respv1`                   |
|`m5`     | monotonic   | ✓   |     |     |     | ✓   | `respv2`                   |
|`m6`     | monotonic   | ✓   | ✓   | ✓   |     | ✓   | `respiv6`                  |
|`m7`     | monotonic   | ✓   | ✓   | ✓   |     |     | `respiv7`                  |
|`m8`     | monotonic   | ✓   | ✓   | ✓   | ✓   | ✓   | `respiv6`                  |
|`m9`     | monotonic   | ✓   | ✓   | ✓   | ✓   |     | `respiv7`                  |
|`m10`    | monotonic   | ✓   |     |     | ✓   | ✓   | `respv2`                   |

<a name="netfootnote">[1]</a> Growth and respiration cannot occur at the same time.

<a name="versionfootnote">[2]</a> See [this notebook](/sizedep_formulations.ipynb).

## Initial tests

 * Tests are underway.
