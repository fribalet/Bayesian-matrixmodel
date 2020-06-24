#!/bin/bash

for f in ../stancode_gallery1/mlmn_version/matrixmodel_mlmultinom_estinilnorm_*delta*_normparam_trackgrowth_xval.stan ; do
    f_new="$(sed -e s/monodelta/monodelta2/ -e s/estinilnorm/estinilnorm2/ -e s/xval/xval2/ <<< "$f")"
    cp -iv "$f" "$(basename $f_new)"
done


