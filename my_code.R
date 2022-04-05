
# create sims
create_sims.R

cd /sc/arion/projects/CommonMind/hoffman/muscat-comparison

screen -S sims -L

ml git python
git pull


# remove results and plots 
# rm -f results/* plots/* logs/*

# remove simulated data
# rm -f data/sim_data/*
# rm -f data/raw_data/sce_nps.rds data/raw_data/ref_nps.rds

# remote snakemake progress
# rm -rf .snakemake/


Rscript setup.R

mkdir plots

snakemake -n

snakemake -j3 --rerun-incomplete

# https://hoffmg01.u.hpc.mssm.edu/sims/

snakemake --jobs 900 --cluster 'bsub -q premium -R "rusage[mem=36000]" -R span[hosts=1] -W 6:00 -P acc_CommonMind -n 2' 



FILE=/sc/arion/projects/CommonMind/hoffman/muscat-comparison/.snakemake/log/2022-04-03T220534.582458.snakemake.log

grep -n Error $FILE


head -n 1926 $FILE | tail -n  20


grep -n Error   | less


args = list(sim = "data/sim_data/nps,de10,5.rds",
fun = "scripts/apply_scdd.R",
wcs = "c=x,did=nps,g=x,i=5,j=1,k=x,mid=scDD.vstresiduals,s=x,sid=de10",
meth_pars = "meta/meth_pars/scDD.vstresiduals.json",
run_pars = "meta/run_pars/nps,de10.json",
res = "results/nps,de10,5,scDD.vstresiduals,1,gx,cx,kx,sx.rds",
c = "x",
did = "nps",
g = "x",
i = "5",
j = "1",
k = "x",
mid = "scDD.vstresiduals",
s = "x",
sid = "de10")

scripts/run_meth.R


args = list(sim = "data/sim_data/nps,de10_ns,1.rds",
fun= "scripts/apply_pb.R",
wcs= "c=x,did=nps,g=x,i=1,j=1,k=x,mid=limma-voom.sum.counts,s=10,sid=de10_ns",
meth_pars=  "meta/meth_pars/limma-voom.sum.counts.json",
run_pars= "meta/run_pars/nps,de10_ns.json",
res= "results/nps,de10_ns,1,limma-voom.sum.counts,1,gx,cx,kx,s10.rds")

R_LIBS_USER=/hpc/users/hoffmg01/.Rlib/R_410/ /hpc/packages/minerva-centos7/R/4.1.0-cairo/lib64/R/bin/R CMD BATCH --no-restore --no-save              "--args sce=data/raw_data/ref_nps.rds ggp=plots/nps-pb_mean_disp.rds fig=plots/nps-pb_mean_disp.pdf"              scripts/plot_pb_mean_disp.R logs/plot_pb_mean_disp-nps.Rout


