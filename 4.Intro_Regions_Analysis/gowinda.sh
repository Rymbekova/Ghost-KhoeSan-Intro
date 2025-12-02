module load Conda HTSlib BEDTools BCFtools R

#### STEP1: PREPARATION OF GENE/GO/GOWINDA INFORMATION (done once for this genome!!)
## in R: preparation of table dataset of all genes, all GO categories, and GTF format for Gowinda
Rscript gowinda_prep.R



#### STEP2: PREPARATION OF RAW DATA
## get files containing only genotypes
# only genotypes
bcftools annotate -x INFO /lisc/scratch/admixlab/aigerim/African/25KS.48RHG.74comp.HCBP_chr"${chrom}"_AF.vcf.gz  | bcftools view -m2 -M2 -v snps | bcftools view -a -S /lisc/scratch/admixlab/mk_data/san/ksp.txt -T /lisc/scratch/admixlab/mk_data/san/allfrag.bed.gz  | bcftools view -m2 -M2 -v snps  | bcftools query -f "%POS [%GT ]\n" | bgzip > /lisc/scratch/admixlab/mk_data/san/chr"${chrom}".gt.gz

## position list
# only positions
for chrom in {1..22}; do
echo $chrom
bcftools query -f "%CHROM\t%POS\n" /lisc/scratch/admixlab/aigerim/African/25KS.48RHG.74comp.HCBP_chr"${chrom}"_AF.vcf.gz \
>> /lisc/scratch/admixlab/mk_data/san/all.pos
done

## get VEP output and turn into nice bed file
awk '$2!~"-" {print $0}' /lisc/scratch/admixlab/aigerim/ensembl-vep/vep_output/chr*.txt | zgrep -v "^#" > /lisc/scratch/admixlab/mk_data/san/tmp.txt
cut -f 2 /lisc/scratch/admixlab/mk_data/san/tmp.txt | sed "s/:/ /g" | awk -v OFS="\t" '{print $0}' | awk -v OFS="\t" '{print $1,$2-1, $2}' | paste - /lisc/scratch/admixlab/mk_data/san/tmp.txt | bgzip -c > /lisc/scratch/admixlab/mk_data/san/anno.bed.gz

# bed file with target coordinates
allfrag=/lisc/scratch/admixlab/mk_data/san/allfrag.bed.gz
# final annotation file
anno=/lisc/scratch/admixlab/mk_data/san/anno.bed.gz

## intersect to get annotations only within target regions
intersectBed -wb -a $allfrag -b $anno | bgzip -c > /lisc/scratch/admixlab/mk_data/san/annosub.txt.gz

## get missense mutations or a simplified file of all mutations
zcat /lisc/scratch/admixlab/mk_data/san/annosub.txt.gz | grep "missense" > /lisc/scratch/admixlab/mk_data/san/allmissense.txt
zcat /lisc/scratch/admixlab/mk_data/san/annosub.txt.gz | cut -f1-3,12,15 > /lisc/scratch/admixlab/mk_data/san/allsimp.txt



#### STEP3: PREPARATION OF GOWINDA DATASET
## again in R: preparation of table dataset of all genes, all GO categories, and GTF format for Gowinda
Rscript gowinda_getgenes.R getdata



#### STEP4: RUNNING GOWINDA ON THE DATASET
module load BEDTools HTSlib BCFtools Java

dadir=/lisc/scratch/admixlab/mk_data/san/
tedir=/lisc/scratch/admixlab/mktest/

# either genwithmut or genwithfun
snp=genwithmut.txt
java -Xmx4g -jar $tedir/Gowinda-1.12.jar --snp-file $dadir/all.pos \
--candidate-snp-file $dadir/$snp --gene-set-file $tedir/hu_annotation.txt \
--annotation-file $tedir/gencode.v49.exons.cor.gtf --simulations 1000000 \
--min-significance 1 --gene-definition gene --threads 8 --output-file $dadir/results_go_gene.txt --mode gene --min-genes 3

snp=genwithfun.txt
java -Xmx4g -jar $tedir/Gowinda-1.12.jar --snp-file $dadir/all.pos \
--candidate-snp-file $dadir/$snp --gene-set-file $tedir/hu_annotation.txt \
--annotation-file $tedir/gencode.v49.exons.cor.gtf --simulations 1000000 \
--min-significance 1 --gene-definition gene --threads 8 --output-file $dadir/results_go_fun.txt --mode gene --min-genes 3



#### STEP5: GET FINAL LIST
## again in R: get gene names & significant categories (0.05 corrected p-value)
Rscript gowinda_getgenes.R getsig
