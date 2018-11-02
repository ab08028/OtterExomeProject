####### Rsync to get mapDamage
# do this from laptop
# don't want the MCMC stats .csv -- too big
# location of github:  which may be on remote server or local laptop
#REFPREFIX=Mustela_putorius_furo.MusPutFur1.0.dna.toplevel
REFPREFIX=sea_otter_23May2016_bS9RH.deduped.99
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
ancHeaders=$wd/samples/ancientLibs.Screen.txt
#modHeaders=$wd/samples/modernSamples.txt

mkdir $wd/transfers # files to transfer to laptop
mkdir $wd/transfers/mapDamagePlots
out=$wd/transfers/mapDamagePlots

cat $ancHeaders | while read header
do
# only copying the misincorp figs, length and MCMC posterior stats
cp $wd/paleomix/testMapping/$header/$header.$REFPREFIX.mapDamage/${header}_1a/Fragmisincorporation_plot.pdf $out/${header}.Fragmisincorporation_plot.pdf
cp $wd/paleomix/testMapping/$header/$header.$REFPREFIX.mapDamage/${header}_1a/Length_plot.pdf $out/${header}.Length_plot.pdf
cp $wd/paleomix/testMapping/$header/$header.$REFPREFIX.mapDamage/${header}_1a/Stats_out_MCMC_post_pred.pdf $out/${header}.Stats_out_MCMC_post_pred.pdf
done


#cat $modHeaders | while read header
#do
# only copying the misincorp figs, length and MCMC posterior stats
#cp $wd/paleomix/$header/$header.$REFPREFIX.mapDamage/${header}_1a/Fragmisincorporation_plot.pdf $out/${header}.Fragmisincorporation_plot.pdf
#cp $wd/paleomix/$header/$header.$REFPREFIX.mapDamage/${header}_1a/Length_plot.pdf $out/${header}.Length_plot.pdf
#cp $wd/paleomix/$header/$header.$REFPREFIX.mapDamage/${header}_1a/Stats_out_MCMC_post_pred.pdf $out/${header}.Stats_out_MCMC_post_pred.pdf
#done
# then want to download the transfers to MapDamage with cyberduck or rsync
