nextflow.enable.dsl=2

// baseline parameters
params.basedir = "/project/guedira_seq_map/nico/SunRILs_population_description"
params.population_table = "${params.basedir}/data/cross_info.csv"
params.output_dir = "${params.basedir}/data/processed_vcf"

// discovery and production parameters
params.fastq_dir  = "/90daydata/guedira_seq_map/nico/fastq"
params.work_dir     = "/90daydata/guedira_seq_map/nico/SunRILs_population_calling"
params.study   = "SunRILs"
params.keyfile      = "${params.basedir}/data/HILLIARD_keyfile.csv"
//params.keyfile      = "${params.basedir}/data/SunRILs_keyfile_merged.csv"
params.ref     = "/90daydata/guedira_seq_map/RefCS_2.1/iwgsc_refseqv2.1_assembly.fa"
params.enzymes      = "PstI-MspI"
params.taglength    = "85"
params.ram          = "300g"
params.ncores       = "40"
params.tassel_run   = " /90daydata/guedira_seq_map/nico/UX1992_reseq_calling/HPC-GBS-Pipeline/tassel-5-standalone/run_pipeline.pl"
params.tassel_pipeline = "${params.basedir}/scripts/external_dependencies/tassel_disc_plus_prod_bwa.sh"

//params.vcf_files = "${params.basedir}/data/raw_vcf"

// Input: fastq directory
Channel.fromPath("${params.work_dir}/fastq/*.{fq,fastq,fastq.gz,fq.gz}").ifEmpty { error "No FASTQ files found" }.set { fastq_files }



// bcftools filtering parameters
params.depth = '3'
params.quality = '20'
params.MAF = '0.05'
params.missing = '0.4'

// individual filtering parameters
//params.hz = '0.2'
//params.missing_line = '0.5'


// imputation parameters
params.monotonic = "${params.basedir}/linkage_map/monotonic"


workflow {
    new File(params.output_dir).mkdirs()

    // Load input VCFs
    Channel.fromPath("${params.fastq_dir}")
        .set { fastq_dir }

    // Broadcast population table
    Channel.value(file(params.population_table))
        .set { pop_table_ch }

    // Step 1: discovery and production
	raw_vcf = tassel_discovery_and_production(fastq_dir)

    // Step 2: Combine filtered VCF with the population table
    vcf_with_pop_table = raw_vcf.combine(pop_table_ch)

    // Step 3: Split filtered VCFs by subpopulations
    subpop_vcfs = split_by_subpop(vcf_with_pop_table)
        .flatten()
        .map { file -> 
            def cross_id = file.getBaseName().tokenize("_")[-2]
            tuple(cross_id, file)
        }

    // Step 4: Group subset VCFs by Cross_ID
    grouped_subpop_vcfs = subpop_vcfs
        .map { cross_id, file -> tuple(cross_id, file) }
        .groupTuple()

	// Step 5: filter and impute vcfs
    grouped_subpop_vcfs
        .map { cross_id, vcf ->
            println "Sending ${vcf} for imputation of ${cross_id}"
            tuple(cross_id, vcf)  // ensure tuple structure
        }
        | bcftools_filter
		| beagle_impute
}

process tassel_discovery_and_production {
    tag "${params.study}"

    input:
    path fastq_dir

    output:
    path "${params.study}_production.vcf.gz"

    script:
    """
    mkdir -p tassel_tmp && cd tassel_tmp

    # Absolute paths
    keyfile=\$(realpath ${params.keyfile})
    fastq=\$(realpath ../${fastq_dir})
    ref=\$(realpath ${params.ref})

    # Ensure slash
    [[ \${fastq: -1} != "/" ]] && fastq="\${fastq}/"

    # Discovery
    \${params.tassel_run} -Xms25g -Xmx${params.ram} -fork1 \
        -GBSSeqToTagDBPlugin -e ${params.enzymes} -i \${fastq} \
        -db ${params.study}.db -k \${keyfile} -kmerLength ${params.taglength} \
        -mnQS ${params.minq} -c 5 -mxKmerNum 50000000 -deleteOldData true \
        -endPlugin -runfork1

    \${params.tassel_run} -Xms25g -Xmx${params.ram} -fork1 \
        -TagExportToFastqPlugin -db ${params.study}.db \
        -o ${params.study}_MasterGBStags.fa.gz -endPlugin -runfork1

    bwa mem -t ${task.cpus} \${ref} ${params.study}_MasterGBStags.fa.gz > ${params.study}.sam

    \${params.tassel_run} -Xms25g -Xmx${params.ram} -fork1 \
        -SAMToGBSdbPlugin -i ${params.study}.sam -db ${params.study}.db \
        -aLen 0 -aProp 0.0 -endPlugin -runfork1

    # Discovery SNP calling
    DISC_ARGS="-fork1 -DiscoverySNPCallerPluginV2 -db ${params.study}.db -mnMAF 0.01 -mnLCov 0.1 -deleteOldData true"

    \${params.tassel_run} -Xms25g -Xmx${params.ram} \$DISC_ARGS -endPlugin -runfork1

    # Production SNP calling
    \${params.tassel_run} -Xms25g -Xmx${params.ram} -fork1 \
        -ProductionSNPCallerPluginV2 -db ${params.study}.db -e ${params.enzymes} \
        -i \${fastq} -k \${keyfile} -kmerLength ${params.taglength} \
        -o ${params.study}_production.vcf -endPlugin -runfork1

    # Compress and index
    bgzip -f ${params.study}_production.vcf
    bcftools index -f ${params.study}_production.vcf.gz

    mv ${params.study}_production.vcf.gz ../
    """
}

process split_by_subpop {
    errorStrategy 'ignore'

    input:
    tuple path(vcf), path(pop_table)

    output:
    path("*_subset.vcf.gz")
	
    script:
    sample = vcf.getBaseName().replaceAll(/\.vcf(\.gz)?$/, "")
    """
    # Extract sample names from the VCF file header
    bcftools query -l ${vcf} > all_samples.txt
    grep '^UX[0-9]\\{4\\}-' all_samples.txt | sed 's/-.*//' | sort -u > cross_ids_present.txt
    awk -F',' 'NR==FNR {keep[\$1]; next} \$1 in keep' cross_ids_present.txt ${pop_table} > subset_pop_table.txt
	

    ##set internal field separator to comma
    IFS=","
    # For each population, subset out the population samples and parents
    while read Cross_ID Parent_1 Parent_2; do
        grep "^\${Cross_ID}-" all_samples.txt > temp_list.txt
        echo "\${Parent_1}" >> temp_list.txt
        echo "\${Parent_2}" >> temp_list.txt
        
        grep -Fxf temp_list.txt all_samples.txt > "\${Cross_ID}_list.txt"

        bcftools view -S \${Cross_ID}_list.txt ${vcf} -Oz -o ${sample}_\${Cross_ID}_subset.vcf.gz
    done < subset_pop_table.txt

    ls -l *.vcf.gz || echo "No matching VCFs created."
    """
}

process bcftools_filter {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    //path vcf
    tuple val(cross_id), path(vcf)

    output:
    //path "${sample}_filt.vcf.gz"
	tuple val(cross_id), path("${sample}_filt.vcf.gz")
	
    script:
    sample = vcf.getBaseName().replaceAll(/\.vcf(\.gz)?$/, "")
    """
    bcftools view -i 'FORMAT/DP > ${params.depth} && 'FORMAT/GQ' >= ${params.quality} && MAF > ${params.MAF} && F_MISSING < ${params.missing}' ${vcf} -Oz -o DP_MAF_MISS_filter.vcf.gz
    bcftools view -m2 -M2 -v snps DP_MAF_MISS_filter.vcf.gz -Oz -o biallelic.vcf.gz
    bcftools index -c biallelic.vcf.gz
    bcftools view -t "^UNKNOWN" biallelic.vcf.gz -Oz -o "${sample}_filt.vcf.gz"
    
    ##create table of marker numbers
    echo -e "VCF File\tTotal Markers" > filtering_marker_table.txt
    for vcf in *.vcf.gz; do
        total_markers=\$(bcftools stats "\$vcf" | grep "^SN" | grep "number of records" | awk '{print \$6}')
        vcf_basename=\$(basename "\$vcf")
        echo -e "\${vcf_basename}\t\${total_markers}" >> filtering_marker_table.txt
    done
    """
}

process beagle_impute {
    publishDir "${params.output_dir}", mode: 'copy'
    memory '100 GB'
    time '24h'

    input:
    tuple val(cross_id), path(vcf)

    output: 
    path "${cross_id}_imp.vcf.gz"

    script:
    """
    beagle gt=${vcf} \
            out=${cross_id}_imp \
            map=${params.monotonic}/${cross_id}_GBS_monotonic.map \
            nthreads=40 \
            window=350
    """
}