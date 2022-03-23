

ml git python

cd /sc/arion/projects/CommonMind/hoffman/muscat-comparison

git pull

Rscript setup.R

snakemake -j2

# https://hoffmg01.u.hpc.mssm.edu/sims/

# remove results and plots 
# rm -f results/* plots/*

snakemake --jobs 200 --cluster 'bsub -q premium -R "rusage[mem=64000]" -R span[hosts=1] -W 12:00 -P acc_CommonMind -n 10'





ml git 
git clone https://github.com/GabrielHoffman/qtlPlots.git
R CMD INSTALL qtlPlots