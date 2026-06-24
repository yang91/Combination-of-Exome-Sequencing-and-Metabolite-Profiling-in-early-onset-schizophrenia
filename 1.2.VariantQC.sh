## Before running, please make sure that you've already added all need tools excuceable path to your environment, including Picard ToolKit, GATK, BWA, Bedtools.
## Please store or use the soft links to guarantee all the seqdata, reference annotations and output results would be stored in the same project directory.

#!/usr/bin/bash
set -euo pipefail

###################################################################################
# Set environment parameters
# input/output path
RESULT_DIR="${HOME}/ExomeSeq-Result/variant"
REF_INDEX_DIR="${HOME}/Ref_and_Index/"

# reference and annotations
REF_FA="${REF_INDEX_DIR}/Homo_sapiens_assembly38.fa"
DBSNP="${REF_INDEX_DIR}/GATK-hg38bundle/dbsnp_144.hg38.withchr.vcf"

# input vcf files
NOR_GT_GVCF="${RESULT_DIR}/GATK/NOR.GenotypeGVCFs.vcf.gz"
SCZ_GT_GVCF="${RESULT_DIR}/GATK/SCZ.GenotypeGVCFs.vcf.gz"
NOR_VQSR_VCF="${RESULT_DIR}/GATK/NOR.SNPs_and_indels.recalibrated.vcf.gz"
SCZ_VQSR_VCF="${RESULT_DIR}/GATK/SCZ.SNPs_and_indels.recalibrated.vcf.gz"

# low complexity bed file
LOW_COMPLEXITY_BED="${REF_INDEX_DIR}/low_complexity_targets.bed"

# log file settings
LOGDIR="./Logs/"
mkdir -p "$LOGDIR"

###################################################################################

# 日志函数
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOGDIR}/pipeline.log"
}

# 检查文件是否存在
check_file() {
    if [[ ! -f "$1" ]]; then
        log "ERROR: File not found: $1"
        exit 1
    fi
}

###################################################################################
# Step1: merge SCZ and NOR vcf files #
STEP=1
log "=== Step $STEP: Merging samples ==="

bcftools merge -0 -m none -o "${RESULT_DIR}/GATK/merge.vcf.gz" -O z "${NOR_GT_GVCF}" "${SCZ_GT_GVCF}"
bcftools merge -0 -m none -o "${RESULT_DIR}/GATK/merge.VQSR.vcf.gz" -O z "${NOR_VQSR_VCF}" "${SCZ_VQSR_VCF}"
log "Step $STEP completed: merge.vcf.gz,  merge.VQSR.vcf.gz"

# Step2: split mulit-allele variants #
STEP=2
log "=== Step $STEP: LeftAlignAndTrimVariants (merge.vcf.gz) ==="

gatk LeftAlignAndTrimVariants -R "$REF_FA" -V ${RESULT_DIR}/GATK/merge.vcf.gz" -O "${RESULT_DIR}/GATK/merge.split-multi-allelics.vcf.gz" \
   --max-indel-length 500 --split-multi-allelics --dont-trim-alleles --keep-original-ac
log "Step $STEP completed"

# Step3: collect variant calling metrics #
STEP=3
log "=== Step $STEP: CollectVariantCallingMetrics (merge) ==="

java -jar picard.jar CollectVariantCallingMetrics I=merge.split-multi-allelics.vcf.gz \
    DBSNP="$DBSNP" O="${RESULT_DIR}/picard/merge.split-multi-allelics_metrics"
log "Step $STEP completed"

# ==================== 步骤4: VQSR后的拆分（注意：需要merge.VQSR.vcf.gz存在） ====================
# Step4: split mulit-allele variants (merge.VQSR.vcf.gz) #
STEP=4
log "=== Step $STEP: LeftAlignAndTrimVariants (VQSR) ==="
gatk LeftAlignAndTrimVariants -R "$REF_FA" -V "${RESULT_DIR}/GATK/merge.VQSR.vcf.gz" -O "${RESULT_DIR}/GATK/merge.VQSR.split-multi-allelics.vcf.gz" \
   --max-indel-length 500 --split-multi-allelics --dont-trim-alleles --keep-original-ac
log "Step $STEP completed"

# Step5: CollectVariantCallingMetric (merge.VQSR.vcf.gz) #
STEP=5
log "=== Step $STEP: CollectVariantCallingMetrics (VQSR) ==="
java -jar picard.jar CollectVariantCallingMetrics I=merge.VQSR.split-multi-allelics.vcf.gz \
    DBSNP="$DBSNP" O="${RESULT_DIR}/picard/merge.VQSR.split-multi-allelics_metrics"
log "Step $STEP completed"

# Step6: Remove low complexity region #
STEP=6
log "=== Step $STEP: Remove low complexity regions ==="
check_file "$LOW_COMPLEXITY_BED"
gatk VariantFiltration -R "$REF_FA" -XL "$LOW_COMPLEXITY_BED" \
    -V "${RESULT_DIR}/GATK/merge.VQSR.split-multi-allelics.vcf.gz" -O "${RESULT_DIR}/GATK/merge.VQSR.low_complexity.vcf.gz"
log "Step $STEP completed"

# Step7: Remove low complexity region #
STEP=7
log "=== Step $STEP: CollectVariantCallingMetrics (low complexity filtered) ==="
java -jar picard.jar CollectVariantCallingMetrics I="${RESULT_DIR}/GATK/merge.VQSR.low_complexity.vcf.gz" \
    DBSNP="$DBSNP" O="${RESULT_DIR}/picard/merge.VQSR.low_complexity_metrics"
log "Step $STEP completed"

# Step8: Hail QC analysis by python#
STEP=8
log "=== Step $STEP: Hail QC Analysis ==="

if [[ ! -f "merge.VQSR.low_complexity.vcf.bgz" ]]; then
    log "Converting vcf.gz to bgz format for Hail..."
    zcat "${RESULT_DIR}/GATK/merge.VQSR.low_complexity.vcf.gz" | bgzip -c > "${RESULT_DIR}/GATK/merge.VQSR.low_complexity.vcf.bgz"
fi

export RESULT_DIR="${RESULT_DIR}"
cat > run_hail_qc.py << "PYEOF"
import hail as hl
import sys

hl.init(log='hail_qc.log')

vcf_path = "${RESULT_DIR}/GATK/merge.VQSR.low_complexity.vcf.bgz"
mt_path = "${RESULT_DIR}/picard/merge.VQSR.low_complexity_metrics'

hl.import_vcf(vcf_path, reference_genome='GRCh38').write(mt_path, overwrite=True)
mt = hl.read_matrix_table(mt_path)

print('variants QC1: Samples: %d  Variants: %d' % (mt.count_cols(), mt.count_rows()))

#Filtering samples with call rate < 0.8
mt = hl.sample_qc(mt)
mt = mt.filter_cols(mt.sample_qc.call_rate >= 0.8)
print('sample QC : Samples: %d  Variants: %d' % (mt.count_cols(), mt.count_rows()))

# Filtering variants with mean depth <10 or >1000
mt = hl.variant_qc(mt)
mt = mt.filter_rows((mt.variant_qc.dp_stats.mean >= 10) & (mt.variant_qc.dp_stats.mean <= 1000))

# Filtering variants with imbalance allele genotype
ab = mt.AD[1] / hl.sum(mt.AD)
filter_condition_ab = (
    (mt.GT.is_hom_ref() & (ab <= 0.1)) |
    (mt.GT.is_het() & (ab >= 0.25) & (ab <= 0.75)) |
    (mt.GT.is_hom_var() & (ab >= 0.9))
)
mt = mt.filter_entries(filter_condition_ab)
print('Genotype QC: Samples: %d  Variants: %d' % (mt.count_cols(), mt.count_rows()))

# Filtering common variants to remain rare variants (MAF < 0.05 and mean GQ >= 25)
mt = hl.variant_qc(mt)
mt = mt.filter_rows((mt.variant_qc.AF[1] < 0.05) & (mt.variant_qc.gq_stats.mean >= 25))
print('rare variants: Samples: %d  Variants: %d' % (mt.count_cols(), mt.count_rows()))

# Export the final after-filtering vcf file
hl.export_vcf(mt, '${RESULT_DIR}/final.QC.maf_0.05.vcf.bgz')
print('Final VCF exported to: ${RESULT_DIR}/final.QC.maf_0.05.vcf.bgz')
PYEOF

python run_hail_qc.py
log "Step $STEP completed: Hail QC finished"
log "=== Pipeline completed successfully! ==="
