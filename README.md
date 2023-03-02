# ADSP_SV
The structural variant (SV) calling pipeline for Alzheimer's Disease Sequencing Project (ADSP) R3_17K and R4_36K whole genome sequence data.

## Dependencies
- [Manta v1.6.0](https://github.com/Illumina/manta/releases/tag/v1.6.0)
- [Smoove v0.2.5](https://github.com/brentp/smoove/releases/tag/v0.2.5)
- [Strelka v2.9.10](https://github.com/Illumina/strelka/releases/tag/v2.9.10)
- [Svimmer](https://github.com/DecodeGenetics/svimmer)
- [GraphTyper2 v2.7.3](https://github.com/DecodeGenetics/graphtyper/releases/tag/v2.7.3)

## Steps
### Sample level calling
1. Apply SV calling programs on each sample
- [Manta](sample_level_calling/Snakefile-Manta)
- [Smoove](sample_level_calling/Snakefile-Smoove)
- [Strelka](sample_level_calling/Snakefile-Strelka) Note: Strelka was not performed for R4_36K.

The complete flow is on [GCAD_SV_pipeline](https://bitbucket.org/ottov123/sv-pipeline/src/master/).

2. Merge VCFs for each sample
- [svimmer](sample_level_calling/svimmer.sh)

### Project level joint genotyping
3. Merge all samples
- [svimmer](Graphtyper2/svimmer_merge.sh)

4. Graphtyper2 joint genotyping
- [Graphtyper2](Graphtyper2/run_graphtyper2.sh)

Note: Joint genotyping was performed for each 100kb region, but centromere.

## License
The implementation is available for academic and nonprofit use for free [LICENSE.md](LICENSE.md).
