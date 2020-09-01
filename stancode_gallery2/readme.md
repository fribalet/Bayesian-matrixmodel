# Model description

## Notes
 * Code compiles and converges well in initial tests.

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
`m11`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_resp_gammaiv6_normparam_trackgrowth_xval2.stan`
`m12`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2_resp_gammaiv7_normparam_trackgrowth_xval2.stan`
`m13`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_resp_gammaiv6_normparam_trackgrowth_xval2.stan`
`m14`:`matrixmodel_mlmultinom_estinilnorm2_monodelta2-lightsig_resp_gammaiv7_normparam_trackgrowth_xval2.stan`

## What all models have in common:
 * All models estimate the initial conditions the same way.
 * All models are currently using the same, updated multinomial to fit the data.
 * All models are using an updated priors compared to the code in [gallery1](/stancode_gallery1).

## Model differences:

| version | core model<sup>[\[1\]](#corefootnote) | `delta_max` | using respiration | size-dep respiration | size-dep growth | light-dep division | using net growth <sup>[\[2\]](#netfootnote) | growth/respiration version <sup>[\[3\]](#versionfootnote) |
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

<a name="versionfootnote">[1]</a> See [this notebook](/sizedep_formulations.ipynb).

## Initial tests

 * Initial tests of models `m1` to `m10` shown in [this notebook](/experimental/exp_zs_20200624_gallery2_test.ipynb).
