models="dadiModel1.OldShallowContraction dadiModel2.RecentExtremeContraction"
for model in $models
do
python make.simulate.${model}.py > macs.Simulation.${model}.sh
done
