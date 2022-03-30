
# create sims
create_sims.R

cd /sc/arion/projects/CommonMind/hoffman/muscat-comparison

screen -S sims -L

ml git python


git pull


# remove results and plots 
# rm -f results/* plots/*

# remove simulated data
# rm -f data/sim_data/*
# rm -f data/raw_data/sce_nps.rds data/raw_data/ref_nps.rds

# remote snakemake progress
# rm -rf .snakemake/


Rscript setup.R

snakemake -j1 --rerun-incomplete

# https://hoffmg01.u.hpc.mssm.edu/sims/

snakemake --jobs 300 --cluster 'bsub -q premium -R "rusage[mem=36000]" -R span[hosts=1] -W 6:00 -P acc_CommonMind -n 2' 





snakemake -n


 R_LIBS_USER=/hpc/users/hoffmg01/.Rlib/R_410/ /hpc/packages/minerva-centos7/R/4.1.0-cairo/lib64/R/bin/R CMD BATCH --no-restore --no-save         "--args res= wcs=did=nps,x=s            ggp=plots/nps-perf_by_ns.rds fig=plots/nps-perf_by_ns.pdf"              scripts/plot_perf_by_nx.R logs/plot_perf_by_ns-nps.Rout
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
    cluster_jobid: Job <59446237> is submitted to queue <premium>.






    