To run Mrbayes, you will need to make .nex file and specify the settings to run Mrbayes.

Example of settings (from Filip and Petra)

```
begin mrbayes;
log start;
prset aamodelpr=fixed(lg) statefreqpr=fixed(empirical);
lset rates=invgamma;
mcmcp ngen=10000000 printfreq=10000 samplefreq=100
mcmcdiagn=yes diagnfreq=500000
nchains=4 savebrlens=yes;
mcmc; sumt;
log stop;
end;
```
