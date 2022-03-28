
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

snakemake -j1

# https://hoffmg01.u.hpc.mssm.edu/sims/

snakemake --jobs 20 --cluster 'bsub -q premium -R "rusage[mem=24000]" -R span[hosts=1] -W 6:00 -P acc_CommonMind -n 6'






data/raw_data/sce_nps.rds



R_LIBS_USER=/hpc/users/hoffmg01/.Rlib/R_410/ /hpc/packages/minerva-centos7/R/4.1.0-cairo/lib64/R/bin/R CMD BATCH --no-restore --no-save            "--args sim=data/sim_data/nps,db10,3.rds fun=scripts/apply_pb.R wcs=c=x,did=nps,g=x,i=3,j=1,k=x,mid=edgeR-treat.sum.counts,s=x,sid=db10               meth_pars=meta/meth_pars/edgeR-treat.sum.counts.json run_pars=meta/run_pars/nps,db10.json res=results/nps,db10,3,edgeR-treat.sum.counts,1,gx,cx,kx,sx.rds"            scripts/run_meth.R logs/run_meth-nps,db10,3,edgeR-treat.sum.counts,1,gx,cx,kx,sx.Rout



