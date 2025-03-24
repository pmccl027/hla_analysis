# HLA-TAPAS pipeline

## NomenCleaner to convert .ped to cleaned .chped HLA type files (to G group resolution)
- `python -m NomenCleaner --hped /path/to/ped --out Ggroup.chped --Ggroup`

### create validation set
- Filter out 3rd degree or > related individuals
    - `plink2 --vcf ../pruned/ALL.common_indep_snps.vcf --king-cutoff 0.01105 --out all_chr.6th-deg`
selected 15 unrelated individuals of each ancestry as validation set to exclude from ref panel

## Prepare snp files
- extract HLA region and filter on quality, missingness
- and include only designated samples (in reference set and with WGS)
    - `bcftools view -S sample_list.txt -r chr6:25000000-35000000 -fPASS -i 'F_MISSING<0.01' -v snps </path/to/SNPs/> > chr6.HLAregion.snps`
## convert to PLINK format
- `plink --vcf BQC19.chr6.HLAregion.snps --set-missing-var-ids @:# --make-bed --out chr6.HLAregion.snps `

## MakeReference to create reference panel
- `sbatch --account=<acct> --time=70:00:00 --ntasks=1 --mem-per-cpu=30G --wrap="python -m MakeReference --variants chr6.HLAregion.snps --chped Ggroup.chped --hg 38 --out MakeReference/Ggroup_out/Ggroup.ref --dict-AA MakeReference/data/hg38/HLA_DICTIONARY_AA.hg38.imgt3320.Ggroup --dict-SNPS MakeReference/data/hg38/HLA_DICTIONARY_SNPS.hg38.imgt3320.Ggroup --phasing --save-intermediates --mem='30g' "`
    - with phasing flag

