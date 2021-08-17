  ## First step for people using WUR's Anunna:
Check if the cache file for your species already exist in `/lustre/nobackup/SHARED/cache/`. If it doesn't, create it with

```
/usr/bin/perl /cm/shared/apps/SHARED/ensembl-vep/INSTALL.pl --CACHEDIR /lustre/nobackup/SHARED/cache/ --AUTO c -n --SPECIES <species>
```
You might need to run it again with `--ASSEMBLY <assembly name>`
