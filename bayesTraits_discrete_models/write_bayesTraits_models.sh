#!/bin/bash


for boot in {0..9}
  do

  cd subset_${boot}

  inputs=bayestraits_*.txt
  tree=filter_combine_1000_SUBSET.phy_ROOT.nexus
  scratch_dir=/scratch/anfij/2023_bayestraits_subsets/subset_${boot}
  bayestraits_dir=/scratch/anfij/soft/BayesTraitsV4.0.1-Linux

  for input in $inputs
    do
    for mod in A B C
      do
      nameis=$(echo $input | sed 's/bayestraits_//g' | sed 's/.txt//g')
      command=command_Discrete_${mod}.txt
      echo $nameis
      rm run_bayesTraits_${nameis}_${mod}.sh
      touch run_bayesTraits_${nameis}_${mod}.sh
      echo "#!/bin/bash
#SBATCH -D "${scratch_dir}"
#SBATCH --account=def-clandry
#SBATCH --job-name="${mod}"-"${nameis}"
#SBATCH --output="${mod}"-"${nameis}"-%j.out
#SBATCH -c 1
#SBATCH --time=3-00:00
#SBATCH --mem=10G
DIR="${bayestraits_dir}"
tree="${tree}"
dataset="${nameis}"
model=model_\${dataset}_"${mod}"
mkdir -p \${model}
for run in run1 run2 run3
  do
  cp bayestraits_\${dataset}.txt ./\${model}/bayestraits_\${dataset}.txt
  \$DIR/BayesTraitsV4 \${tree} ./\${model}/bayestraits_\${dataset}.txt < "${command}"
  mv ./\${model}/bayestraits_\${dataset}.txt.Log.txt ./\${model}/bayestraits_\${dataset}_\${run}.txt.Log.txt
  mv ./\${model}/bayestraits_\${dataset}.txt.Schedule.txt ./\${model}/bayestraits_\${dataset}_\${run}.txt.Schedule.txt
  mv ./\${model}/bayestraits_\${dataset}.txt.Stones.txt ./\${model}/bayestraits_\${dataset}_\${run}.txt.Stones.txt
  done" > run_bayesTraits_${nameis}_${mod}.sh
      done
    done

    cd ..

  done
