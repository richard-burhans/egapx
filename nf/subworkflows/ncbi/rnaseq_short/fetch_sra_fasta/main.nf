#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'
params.inpdir = ""
workflow fetch_sra_fasta  {
    take:
        sra_run_file  // Channel: sra run file
        parameters     // Map : extra parameter and parameter update
    main:
        a = read_sra_file(sra_run_file)
        b =  a.map { it -> it.split()}
        c =  run_fetch_sra_fasta(b.flatten())

    emit:
        //fasta_pair_list = run_fetch_sra_fasta.out.fasta_pair_list
        fasta_pair_list = c
}


process read_sra_file {
    input:
        path sra_run_file
    output:
        env  exitvar
    script:

    """
    exitvar=()
    while read -r line; do  [[ \$line = \\#* ]] && continue; exitvar+=(\"\$line\"); done < ${sra_run_file}
    """

    stub:
    """
    exitvar="SRA000001 SRA000002"
    """
}

process run_fetch_sra_fasta {
    input:
        val sra
    output:
        tuple val (sra),  path ('output/*_{1,2}.fasta')  , emit: 'fasta_pair_list'
    script:
    """
    mkdir -p output
    curl -fL --retry 5 -C - -o ${sra}.sra \$(srapath ${sra})
    fasterq-dump --skip-technical --threads 6 --split-files --seq-defline ">gnl|SRA|\\\$ac.\\\$si.\\\$ri" --fasta --outdir output  ./${sra}.sra
    ls output/${sra}_*.fasta
    rm -f ${sra}.sra
    """
    stub:
    """
    mkdir -p output
    touch output/${sra}_1.fasta
    touch output/${sra}_2.fasta
    """
}
