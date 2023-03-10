#!/bin/sh
gatk='singularity run --cleanenv -B /nas:/nas -B /ruby:/ruby -B /the_dysk:/the_dysk /nas/usr/local/opt/gatk4.simg gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx16G"'
ref='/ruby/WES_reference/hg38/Homo_sapiens_assembly38'
echo ""
echo "bed to interval 1 and 2"
echo ""
${gatk} BedToIntervalList --INPUT KAPA_HyperExome_hg38_capture_targets.bed --OUTPUT KAPA_HyperExome_hg38_bait.interval_list --SEQUENCE_DICTIONARY ${ref}.dict
${gatk} BedToIntervalList --INPUT KAPA_HyperExome_hg38_primary_targets.bed --OUTPUT KAPA_HyperExome_hg38_target.interval_list --SEQUENCE_DICTIONARY ${ref}.dict
echo ""
echo "cat"
echo ""
cat KAPA_HyperExome_hg38_capture_targets.bed KAPA_HyperExome_hg38_primary_targets.bed | sort -k1,1 -k2,2n |  /nas/software/bedtools2/bin/bedtools merge -i - > KAPA_HyperExome_hg38_capture_primary_targets_union.bed
echo ""
echo "bed to interval 3"
echo ""
${gatk} BedToIntervalList --INPUT KAPA_HyperExome_hg38_capture_primary_targets_union.bed --OUTPUT KAPA_HyperExome_hg38_union.interval_list --SEQUENCE_DICTIONARY ${ref}.dict
