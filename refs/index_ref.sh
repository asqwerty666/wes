#!/bin/sh
ref=$1
shift
/nas/usr/local/bin/bwa index -a bwtsw ${ref}
/nas/software/samtools/bin/samtools faidx ${ref}
singularity run --cleanenv -B /nas:/nas -B /ruby:/ruby -B /the_dysk:/the_dysk /nas/usr/local/opt/gatk4.simg gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx16G"  CreateSequenceDictionary --REFERENCE ${ref}
