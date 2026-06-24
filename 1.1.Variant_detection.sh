## Please run this pipeline for patient samples and normal sampls seperately
## Before running, please make sure that you've already added all need tools excuceable path to your environment, including Picard ToolKit, GATK, BWA, Bedtools.
## Please store or use the soft links to guarantee all the seqdata, reference annotations and output results would be stored in the same project directory.

#!/usr/bin/bash
# bwa + bedtools + picard + GATK

# tag
GROUP="SCZ" # modify this to NOR if you are running normal sample analysis

# input/output path
PROJECT_DIR="${HOME}/EOSCZ"
SEQDATA_DIR="${HOME}/EOSCZ/ExomeSeq-Data"
RESULT_DIR="${HOME}/EOSCZ/ExomeSeq-Result"
REF_INDEX_DIR="${HOME}/EOSCZ/Ref_and_Index/"

# Reference/index path
BWA_INDEX="${REF_INDEX_DIR}/BWA_INDEX"
BED_FILE="${REF_INDEX_DIR}/T086V4_MT.merged.success.liftover.to.hg38.bed"
REF_FA="${REF_INDEX_DIR}/GATK-hg38bundle/Homo_sapiens_assembly38.fa"
DBSNP="${REF_INDEX_DIR}/GATK-hg38bundle/dbsnp_144.hg38.withchr.vcf"
KNOWN_INDEL="${REF_INDEX_DIR}/GATK-hg38bundle/Homo_sapiens_assembly38.known_indels.vcf"
GOLDEN_INDEL="${REF_INDEX_DIR}/GATK-hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf"
HAPMAP="${REF_INDEX_DIR}/GATK-hg38bundle/hapmap_3.3.hg38.vcf.gz"
OMNI="${REF_INDEX_DIR}/GATK-hg38bundle/1000G_omni2.5.hg38.vcf.gz"
OKG="${REF_INDEX_DIR}/GATK-hg38bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
AXIOM="${REF_INDEX_DIR}/GATK-hg38bundle/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"

# parameters
THREADS=6

SAMPLE_LIST="${HOME}/EOSCZ/${GROUP}.sample_list.txt" 
# sample_list.txt contains all patient samples or all normal samples name used for exome-seq analysis
# please not list both patient sample and normal sample in the sample_list.txt at the same time
mapfile -t SAMPLES < "${SAMPLE_LIST}"
TOTAL=${#SAMPLES[@]}
COUNT=0

for SAMPLE in "${SAMPLES[@]}"; do
    COUNT=$((COUNT + 1))
    echo "============================================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Start to dealing with ${SAMPLE} (${COUNT}/${TOTAL}) ..."
    echo "============================================================"

    # Define the input fastq files
    R1="${SEQDATA_DIR}/${SAMPLE}.R1.fastq.gz"
    R2="${SEQDATA_DIR}/${SAMPLE}.R2.fastq.gz"

    # Check if the input fastq file exists
    if [[ ! -f "${R1}" ]]; then
        echo "[ERROR] Cannot find R1 file: ${R1}" >&2
        continue
    fi
    if [[ ! -f "${R2}" ]]; then
        echo "[ERROR] Cannot find R2 file: ${R2}" >&2
        continue
    fi

    # Define and build the output directory if not exists
    BWA_DIR="${RESULT_DIR}/bwa"
    PICARD_DIR="${RESULT_DIR}/picard/${SAMPLE}"
    BEDTOOLS_DIR="${RESULT_DIR}/bedtools"
    GATK_DIR="${RESULT_DIR}/GATK/${SAMPLE}"

    mkdir -p "${BWA_DIR}" "${PICARD_DIR}" "${BEDTOOLS_DIR}" "${GATK_DIR}"

    # DEFINE the intermediate file
    SAM="${BWA_DIR}/${SAMPLE}.aln-pe.sam"
    UNALI_BAM="${PICARD_DIR}/unaligned_read_pairs.bam"
    MERGED_BAM="${PICARD_DIR}/MergeBamAlignment.bam"
    MARKDUP_BAM="${PICARD_DIR}/Mark_duplicates.bam"
    MARKDUP_METRICS="${PICARD_DIR}/marked_dup_metrics.txt"
    GENOMECOV="${BEDTOOLS_DIR}/${SAMPLE}.genomecov"
    RECAL_TABLE="${GATK_DIR}/recal_data.table"
    BQSR_BAM="${GATK_DIR}/AfterBQSR.bam"
    GVCF="${GATK_DIR}/HaplotypeCaller.g.vcf.gz"

    # Log file
    LOG="${PROJECT_DIR}/logs/${SAMPLE}.log"
    mkdir -p "${PROJECT_DIR}/logs"

    #-------------------------------------------------------------------------------
    # Step 1: BWA MEM alignment
    #-------------------------------------------------------------------------------
    echo "[${SAMPLE}] Step 1/6: BWA MEM alignment ..."
    if [[ ! -f "${SAM}" ]]; then
        bwa mem -t ${THREADS} -M \
            "${BWA_INDEX}" "${R1}" "${R2}" > "${SAM}" 2>> "${LOG}"
        
        if [[ $? -ne 0 ]]; then
            echo "[ERROR_1] Failed to run BWA alignment: ${SAMPLE}" >&2
            continue
        fi
    else
        echo "[${SAMPLE}] SAM already exists，skip the BWA alignment."
    fi

    #-------------------------------------------------------------------------------
    # Step 2: Picard MarkDuplicates
    #-------------------------------------------------------------------------------
    echo "[${SAMPLE}] Step 2/6: Picard MarkDuplicates ..."
    if [[ ! -f "${UNALI_BAM}" ]]; then
        java -jar picard.jar FastqToSam FASTQ="${R1}" FASTQ2="${R2}" \
            OUTPUT="{$UNALI_BAM}" READ_GROUP_NAME=A00682 SAMPLE_NAME=${SAMPLE} \
            LIBRARY_NAME=Illumina PLATFORM_UNIT=HJVWVDSXX PLATFORM=Illumina \
             >> "${LOG}" 2>&1
        
        if [[ $? -ne 0 ]]; then
            echo "[ERROR_2] Failed to run Picard FastqToSam: ${SAMPLE}" >&2
            continue
        fi
    else
        echo "[${SAMPLE}] FastqToSam already exists，skip"
    fi

    if [[ ! -f "${MERGED_BAM}" ]]; then
        java -jar picard.jar MergeBamAlignment ALIGNED="${SAM}" UNMAPPED="${UNALI_BAM}" \
        O="${MERGED_BAM}" R="${REF_FA}" >> "${LOG}" 2>&1
        
        if [[ $? -ne 0 ]]; then
            echo "[ERROR_3] Failed to run Picard MergeBamAlignment: ${SAMPLE}" >&2
            continue
        fi
    else
        echo "[${SAMPLE}] FastqToSam already exists，skip"
    fi

    if [[ ! -f "${MARKDUP_BAM}" ]]; then
        java -jar picard.jar MarkDuplicates I="${MERGED_BAM}" O="${MARKDUP_BAM}" \
            M="${MARKDUP_METRICS}" TAGGING_POLICY=All ASSUME_SORT_ORDER=coordinate \
            CREATE_INDEX=TRUE >> "${LOG}" 2>&1
        
        if [[ $? -ne 0 ]]; then
            echo "[ERROR_4] Failed to run Picard MarkDuplicates : ${SAMPLE}" >&2
            continue
        fi
    else
        echo "[${SAMPLE}] MarkDuplicates already exists，skip"
    fi

    #-------------------------------------------------------------------------------
    # Step 3: bedtools genomecov to estimate the genome coverage
    #-------------------------------------------------------------------------------
    echo "[${SAMPLE}] Step 3/6: bedtools genomecov ..."
    if [[ ! -f "${GENOMECOV}" ]]; then
        bedtools genomecov -ibam "${MARKDUP_BAM}" -bga -g "${REF_FA}" > "${GENOMECOV}" \
            2>> "${LOG}"
        
        if [[ $? -ne 0 ]]; then
            echo "[ERROR_5] Faild to run bedtools genomecov: ${SAMPLE}" >&2
            continue
        fi
    else
        echo "[${SAMPLE}] genomecov already exists, skip"
    fi

    #-------------------------------------------------------------------------------
    # Step 4: GATK for single sample
    #-------------------------------------------------------------------------------
    echo "[${SAMPLE}] Step 4/6: GATK for single sample ..."
    if [[ ! -f "${RECAL_TABLE}" ]]; then
        gatk BaseRecalibrator -R "${REF_FA}" -I "${MARKDUP_BAM}" -O "${RECAL_TABLE}" \
            --known-sites "${DBSNP}" --known-sites "${KNOWN_INDEL}" --known-sites "${GOLDEN_INDEL}" \
            -L "${BED_FILE}" --interval-padding 200 >> "${LOG}" 2>&1
        
        if [[ $? -ne 0 ]]; then
            echo "[ERROR_6.1] Failed to run GATK BaseRecalibrator: ${SAMPLE}" >&2
            continue
        fi
    else
        echo "[${SAMPLE}] GATK Recalibration result already exists, skip"
    fi

    if [[ ! -f "${BQSR_BAM}" ]]; then
        gatk ApplyBQSR -R "${REF_FA}" -I "${MARKDUP_BAM}" --bqsr-recal-file "${RECAL_TABLE}" \
            -O "${BQSR_BAM}" >> "${LOG}" 2>&1
        
        if [[ $? -ne 0 ]]; then
            echo "[ERROR_6.2] Failed to run GATK ApplyBQSR: ${SAMPLE}" >&2
            continue
        fi
    else
        echo "[${SAMPLE}] GATK BQSR result already exists, skip"
    fi

    if [[ ! -f "${GVCF}" ]]; then
        gatk --java-options "-Xmx4g" HaplotypeCaller -R "${REF_FA}" \
            -I "${BQSR_BAM}" -O "${GVCF}" -ERC GVCF -G StandardAnnotation \
            -G StandardHCAnnotation -G AS_StandardAnnotation \
            -L "${BED_FILE}" --interval-padding 200  >> "${LOG}" 2>&1
        
        if [[ $? -ne 0 ]]; then
            echo "[ERROR_6.3] Failed to run GATK HaplotypeCaller: ${SAMPLE}" >&2
            continue
        fi
    else
        echo "[${SAMPLE}] GATK HaplotypeCaller result already exists, skip"
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${SAMPLE} 处理完成！"
    echo ""

done

echo "============================================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished (${TOTAL} single file process.)"
echo "============================================================"


# Step 5: GATK CombineGVCF + GenotypeGVCF + GatherVCF + VQSR
echo "============================================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Start to collaspe GVCFs of all samples and run GATK further analysis ..."
echo "============================================================"

COMBINE_DIR="${RESULT_DIR}/GATK/"
COMBINED_GVCF="${COMBINE_DIR}/${GROUP}.combined.g.vcf.gz"
GT_GVCF="${COMBINE_DIR}/${GROUP}.GenotypeGVCFs.vcf.gz"

mkdir -p "${COMBINE_DIR}"

if [[ -f "${COMBINED_GVCF}" ]]; then
    echo "[Notice] Combined GVCF already exists: ${COMBINED_GVCF}"
    echo "If you should re-combine all samples, please delete this file at first."
else
    VARIANT_ARGS=()
    MISSING_GVCF=0

    for SAMPLE in "${SAMPLES[@]}"; do
        GVCF="${RESULT_DIR}/GATK/${SAMPLE}/HaplotypeCaller.g.vcf.gz"
        
        if [[ -f "${GVCF}" ]]; then
            VARIANT_ARGS+=("--variant" "${GVCF}")
        else
            echo "[Warning] Cannot find GVCF 文件，this sample will be skipped: ${GVCF}" >&2
            MISSING_GVCF=$((MISSING_GVCF + 1))
        fi
    done

    if [[ ${#VARIANT_ARGS[@]} -eq 0 ]]; then
        echo "[ERROR_7] No avaliable GVCF file was found, CombineGVCFs will not be proceed!" >&2
        exit 1
    fi

    if [[ ${MISSING_GVCF} -gt 0 ]]; then
        echo "[Warning] ${MISSING_GVCF} GVCF were not found，remained $(( ${#SAMPLES[@]} - ${MISSING_GVCF} )) will be proceed in CombineGVCF ..."
    fi

    echo "[Information] CombineGVCF ${#VARIANT_ARGS[@]} files ..."

    gatk --java-options "-Xmx4g -Xms4g" CombineGVCFs -R "${REF_FA}" \
        "${VARIANT_ARGS[@]}" \
        -O "${COMBINED_GVCF}" >> "${COMBINE_DIR}/CombineGVCFs.log" 2>&1

    if [[ $? -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] GATK CombineGVCFs finished!"
    else
        echo "[ERROR_7] Failed to run GATK CombineGVCFs! Detailed logs in ${COMBINE_DIR}/CombineGVCFs.log" >&2
        exit 1
    fi
fi

gatk GenotypeGVCFs --java-options "-Xmx4g" -R "${REF_FA}" \
  -V "${COMBINED_GVCF}" -O "${GT_GVCF}" >> "${COMBINE_DIR}/GenotypeGVCFs.log" 2>&1

gatk VariantRecalibrator --java-options "-Xmx4g -Xms4g" \
  -V "${GT_GVCF}" \
  -O "${COMBINE_DIR}/${GROUP}.SNPs.recal" \
  --tranches-file "${COMBINE_DIR}/${GROUP}.SNPs.tranches" \
  --trust-all-polymorphic -tranche 100 -tranche 99.95 -tranche 99.9 \
  -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 \
  -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
  -mode SNP --max-gaussians 6 \
  --resource:hapmap,known=false,training=true,truth=true,prior=15 "${HAPMAP}" \
  --resource:omni,known=false,training=true,truth=true,prior=12 "${OMNI}" \
  --resource:1000G,known=false,training=true,truth=false,prior=10 "${OKG}" \
  --resource:dbsnp,known=true,training=false,truth=false,prior=7 "${DBSNP}" \
  >> "${COMBINE_DIR}/VQSR.log" 2>&1

gatk VariantRecalibrator --java-options "-Xmx24g -Xms24g" \
  -V "${GT_GVCF}" \
  -O "${COMBINE_DIR}/${GROUP}.indels.recal" \
  --tranches-file "${COMBINE_DIR}/${GROUP}.INDELs.tranches" \
  --trust-all-polymorphic -tranche 100 -tranche 99.95 -tranche 99.9 \
  -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 \
  -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
  -mode INDEL --max-gaussians 4 \
  --resource:mills,known=false,training=true,truth=true,prior=12 "${GOLDEN_INDEL}" \
  --resource:axiomPoly,known=false,training=true,truth=false,prior=10 "${AXIOM}" \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2 "${DBSNP}" \
  >> "${COMBINE_DIR}/VQSR.log" 2>&1

gatk ApplyVQSR --java-options "-Xmx5g -Xms5g" \
  -O "${COMBINE_DIR}/${GROUP}.indels.recalibrated.vcf" \
  -V "${GT_GVCF}" \
  --recal-file "${COMBINE_DIR}/${GROUP}.indels.recal" \
  --tranches-file "${COMBINE_DIR}/${GROUP}.INDELs.tranches" \
  --truth-sensitivity-filter-level 99.9 --create-output-variant-index true -mode INDEL \
  >> "${COMBINE_DIR}/VQSR.log" 2>&1

gatk ApplyVQSR --java-options "-Xmx5g -Xms5g" \
  -O "${COMBINE_DIR}/${GROUP}.SNPs_and_indels.recalibrated.vcf" \
  -V "${COMBINE_DIR}/${GROUP}.indels.recalibrated.vcf" \
  --recal-file "${COMBINE_DIR}/${GROUP}.indels.recal" \
  --tranches-file "${COMBINE_DIR}/${GROUP}.SNPs.tranches" \
  --truth-sensitivity-filter-level 99.9 --create-output-variant-index true \
  >> "${COMBINE_DIR}/VQSR.log" 2>&1

echo ""
echo "============================================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] All processes to detect variants from ${GROUP} have been finished!"
