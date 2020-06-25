# Merge and filter vole resequencing data

module load bcftools/1.9
 
### 1. Filter individual VCFs
# with custom DP, GQ, FS, MQ, MQ0 filters
bcftools view -e '(INFO/DP < 81) | (INFO/DP > 146) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase1/call_stampy_s05/vcf//Eugene_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/Eugene_UnifiedGenotyper.filtered.vcf.gz --threads 2
bcftools view -e '(INFO/DP < 5) | (INFO/DP > 24) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase2/call_stampy_s05/vcf//MagCHRi22_wh_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/MagCHRi22_wh_UnifiedGenotyper.filtered.vcf.gz --threads 2
bcftools view -e '(INFO/DP < 10) | (INFO/DP > 36) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase2/call_stampy_s05/vcf//MarCHBo17_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/MarCHBo17_UnifiedGenotyper.filtered.vcf.gz --threads 2
bcftools view -e '(INFO/DP < 14) | (INFO/DP > 43) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase2/call_stampy_s05/vcf//MarCZD02_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/MarCZD02_UnifiedGenotyper.filtered.vcf.gz --threads 2
bcftools view -e '(INFO/DP < 5) | (INFO/DP > 23) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase2/call_stampy_s05/vcf//MarFTh497_wh_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/MarFTh497_wh_UnifiedGenotyper.filtered.vcf.gz --threads 2
bcftools view -e '(INFO/DP < 8) | (INFO/DP > 39) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase2/call_stampy_s05/vcf//MbrRClm01_wh_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/MbrRClm01_wh_UnifiedGenotyper.filtered.vcf.gz --threads 2
bcftools view -e '(INFO/DP < 5) | (INFO/DP > 16) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase2/call_stampy_s05/vcf//Mcab10_wh_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/Mcab10_wh_UnifiedGenotyper.filtered.vcf.gz --threads 2
bcftools view -e '(INFO/DP < 9) | (INFO/DP > 37) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase2/call_stampy_s05/vcf//MduP01_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/MduP01_UnifiedGenotyper.filtered.vcf.gz --threads 2
bcftools view -e '(INFO/DP < 5) | (INFO/DP > 28) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase2/call_stampy_s05/vcf//MglCHBk23_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/MglCHBk23_UnifiedGenotyper.filtered.vcf.gz --threads 2
bcftools view -e '(INFO/DP < 5) | (INFO/DP > 21) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase2/call_stampy_s05/vcf//MluP01_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/MluP01_UnifiedGenotyper.filtered.vcf.gz --threads 2
bcftools view -e '(INFO/DP < 5) | (INFO/DP > 28) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase2/call_stampy_s05/vcf//MoePBi05_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/MoePBi05_UnifiedGenotyper.filtered.vcf.gz --threads 2
bcftools view -e '(INFO/DP < 5) | (INFO/DP > 30) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase2/call_stampy_s05/vcf//MpeMi11_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/MpeMi11_UnifiedGenotyper.filtered.vcf.gz --threads 2
bcftools view -e '(INFO/DP < 11) | (INFO/DP > 40) | GQ < 30 | FS > 26 | MQ < 25 | (MQ0 >= 4 & (MQ0/(1.0*INFO/DP) > 0.1))' /storage/data/Voles/phase2/call_stampy_s05/vcf//MroFl38_UnifiedGenotyper.vcf.gz -Oz -o /home/gouy/res/MroFl38_UnifiedGenotyper.filtered.vcf.gz --threads 2

### 2. Index individual VCFs
tabix -p vcf /home/gouy/res/Eugene_UnifiedGenotyper.filtered.vcf.gz
tabix -p vcf /home/gouy/res/MagCHRi22_wh_UnifiedGenotyper.filtered.vcf.gz
tabix -p vcf /home/gouy/res/MarCHBo17_UnifiedGenotyper.filtered.vcf.gz
tabix -p vcf /home/gouy/res/MarCZD02_UnifiedGenotyper.filtered.vcf.gz
tabix -p vcf /home/gouy/res/MarFTh497_wh_UnifiedGenotyper.filtered.vcf.gz
tabix -p vcf /home/gouy/res/MbrRClm01_wh_UnifiedGenotyper.filtered.vcf.gz
tabix -p vcf /home/gouy/res/Mcab10_wh_UnifiedGenotyper.filtered.vcf.gz
tabix -p vcf /home/gouy/res/MduP01_UnifiedGenotyper.filtered.vcf.gz
tabix -p vcf /home/gouy/res/MglCHBk23_UnifiedGenotyper.filtered.vcf.gz
tabix -p vcf /home/gouy/res/MluP01_UnifiedGenotyper.filtered.vcf.gz
tabix -p vcf /home/gouy/res/MoePBi05_UnifiedGenotyper.filtered.vcf.gz
tabix -p vcf /home/gouy/res/MpeMi11_UnifiedGenotyper.filtered.vcf.gz
tabix -p vcf /home/gouy/res/MroFl38_UnifiedGenotyper.filtered.vcf.gz

### 3. Merge filtered individual VCFs and index resulting file
# Remove sites for which REF is N
# -m snps: allow multiallelic SNPs

bcftools merge --threads 2 -m snps -f 'REF!="N"'  \
	/home/gouy/res/Eugene_UnifiedGenotyper.filtered.vcf.gz \
	/home/gouy/res/MagCHRi22_wh_UnifiedGenotyper.filtered.vcf.gz \
	/home/gouy/res/MarCHBo17_UnifiedGenotyper.filtered.vcf.gz \
	/home/gouy/res/MarCZD02_UnifiedGenotyper.filtered.vcf.gz \
	/home/gouy/res/MarFTh497_wh_UnifiedGenotyper.filtered.vcf.gz \
	/home/gouy/res/MbrRClm01_wh_UnifiedGenotyper.filtered.vcf.gz \
	/home/gouy/res/Mcab10_wh_UnifiedGenotyper.filtered.vcf.gz \
	/home/gouy/res/MduP01_UnifiedGenotyper.filtered.vcf.gz \
	/home/gouy/res/MglCHBk23_UnifiedGenotyper.filtered.vcf.gz \
	/home/gouy/res/MluP01_UnifiedGenotyper.filtered.vcf.gz \
	/home/gouy/res/MoePBi05_UnifiedGenotyper.filtered.vcf.gz \
	/home/gouy/res/MpeMi11_UnifiedGenotyper.filtered.vcf.gz \
	/home/gouy/res/MroFl38_UnifiedGenotyper.filtered.vcf.gz -Oz -o /home/gouy/res/13V.vcf.gz

tabix -p vcf /home/gouy/res/13V.vcf.gz

### 4. Remove sites with at least 1 missing genotype, keep diallelic SNPs and index file
bcftools view --threads 2 -a \
	-r ScOZjSD_3547,ScOZjSD_2324,ScOZjSD_3553,ScOZjSD_2488,ScOZjSD_3549,ScOZjSD_3244,ScOZjSD_727,ScOZjSD_3548,ScOZjSD_2843,ScOZjSD_3505,ScOZjSD_2119,ScOZjSD_2657,ScOZjSD_1113,ScOZjSD_541,ScOZjSD_917,ScOZjSD_3552,ScOZjSD_748,ScOZjSD_1489,ScOZjSD_3546,ScOZjSD_3550,ScOZjSD_3551,ScOZjSD_436 \
	--threads 2  -f '.' -c 1:minor -m2 -M2 -g ^miss -Oz -o /home/gouy/res/13V-22sc-Di-Nomiss.vcf.gz /home/gouy/res/13V.vcf.gz
tabix -p vcf /home/gouy/res/13V-22sc-Di-Nomiss.vcf.gz

### 5. Remove all fields except "GT" and index file
bcftools annotate --threads 2 -x INFO,^FORMAT/GT -Oz -o /home/gouy/res/13V-22sc-Di-Nomiss-GT.vcf.gz /home/gouy/res/13V-22sc-Di-Nomiss.vcf.gz
tabix -p vcf /home/gouy/res/13V-22sc-Di-Nomiss-GT.vcf.gz
