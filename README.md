# ADSP_SV
The structural variant (SV) calling pipeline for Alzheimer's Disease Sequencing Project (ADSP) R3_17K and R4_36K whole genome sequence data.


## Step 1 : Sample level calling
- [Manta](sample_level_calling/Snakefile-Manta)
- [Smoove](sample_level_calling/Snakefile-Smoove)
- [Strelka](sample_level_calling/Snakefile-Strelka)

Note: Strelka was not performed for R4_36K. The complete flow is on [GCAD_SV_pipeline](https://bitbucket.org/ottov123/sv-pipeline/src/master/).

## Step 2: Sample level VCF merging
- [svimmer](sample_level_calling/svimmer.sh)

## Step 3: Project level VCF merging
- [svimmer](Graphtyper2/svimmer_merge.sh)

## Step 4: Graphtyper2 joint genotyping

## License
The implementation is available for academic and nonprofit use for free [LICENSE.md](LICENSE.md).
