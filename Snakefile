configfile: "config.yaml"

chroms = [str(x) for x in range(1,24)] + ['X']
chroms = ['chr' + i for i in chroms]

rule all:
    input:
        vcf=f"output/{config['project_name']}.{config['new_build']}.all_chrs.vcf.gz",
        index=f"output/{config['project_name']}.{config['new_build']}.all_chrs.vcf.gz.csi"

rule rename_chr:
    input: config['target_vcf_path']

    output: f"output/{config['project_name']}.{{old_build}}.rename_chrs.vcf.gz"

    shell:
        """
        bcftools annotate --rename-chrs <(for chr in {{1..23}} X; do echo -e $chr"\t"chr"${{chr}}"; done) {input}  | bcftools view -Oz > {output}
        """

rule index_rename:
    input: f"output/{config['project_name']}.{{old_build}}.rename_chrs.vcf.gz"

    output: f"output/{config['project_name']}.{{old_build}}.rename_chrs.vcf.gz.csi"

    shell:
        """
        bcftools index {input}
        """

rule split_chr:
    input:
        vcf=f"output/{config['project_name']}.{{old_build}}.rename_chrs.vcf.gz",
        index=f"output/{config['project_name']}.{{old_build}}.rename_chrs.vcf.gz.csi"

    output: f"output/{config['project_name']}.{{old_build}}.rename_chrs.{{chrom}}.vcf.gz"

    shell:
        """
        bcftools view -r {wildcards.chrom} {input.vcf} -Oz > {output}
        """

rule get_ref:
    output: ref_path=temp(f"resources/{config['new_build']}.fasta")

    params:
        zipped=f"resources/{config['new_build']}.fa.gz"

    shell:
        f"""
        wget https://hgdownload.cse.ucsc.edu/goldenpath/{config['new_build']}/bigZips/{config['new_build']}.fa.gz -P resources
        zcat {{params.zipped}} > {{output.ref_path}}
        """

rule get_chain:
    output:
        chain_path=f"resources/{config['old_build']}To{config['new_build'].capitalize()}.over.chain.gz"

    shell:
        f"""
        wget https://hgdownload.cse.ucsc.edu/goldenpath/{config['old_build']}/liftOver/{config['old_build']}To{config['new_build'].capitalize()}.over.chain.gz -P resources
        """

rule seq_dict_path:
    input:
        ref_path=f"resources/{config['new_build']}.fasta"

    output:
        dict_path=f"resources/{config['new_build']}.fasta.dict"

    shell:
        """
        java -jar -Xmx32g bin/picard.jar \
        CreateSequenceDictionary \
        R={input.ref_path} \
        O={output.dict_path} 
        """
#
rule liftover:
    input:
        target_vcf_path=f"output/{config['project_name']}.{config['old_build']}.rename_chrs.{{chrom}}.vcf.gz",
        ref_path=f"resources/{config['new_build']}.fasta",
        dict_path=f"resources/{config['new_build']}.fasta.dict",
        chain_path=f"resources/{config['old_build']}To{config['new_build'].capitalize()}.over.chain.gz"

    output:
        rejected_vcf_path=temp(f"output/{config['project_name']}.rejected.{{chrom}}.vcf"),
        lifted_vcf_path=temp(f"output/{config['project_name']}.{config['new_build']}.{{chrom}}.vcf"),
        lifted_vcf_index_path=temp(f"output/{config['project_name']}.{config['new_build']}.{{chrom}}.vcf.idx")

    shell:
        """
         java -jar -Xmx64g bin/picard.jar LiftoverVcf \
         I={input.target_vcf_path} \
         O={output.lifted_vcf_path} \
         CHAIN={input.chain_path} \
         REJECT={output.rejected_vcf_path} \
         MAX_RECORDS_IN_RAM=10000 \
         R={input.ref_path}
        """

rule bgzip:
    input: f"output/{config['project_name']}.{config['new_build']}.{{chrom}}.vcf"

    output: f"output/{config['project_name']}.{config['new_build']}.{{chrom}}.vcf.gz"

    shell:
        """
        bcftools sort {input} | bcftools view -Oz > {output}
        """

rule index_sorted:
    input: f"output/{config['project_name']}.{config['new_build']}.{{chrom}}.vcf.gz"

    output: f"output/{config['project_name']}.{config['new_build']}.{{chrom}}.vcf.gz.csi"

    shell:
        """
        bcftools index {input} 
        """

rule concat:
    input:
        vcf=expand(f"output/{config['project_name']}.{config['new_build']}.{{chrom}}.vcf.gz", chrom=chroms),
        index=expand(f"output/{config['project_name']}.{config['new_build']}.{{chrom}}.vcf.gz.csi",chrom=chroms)

    output: f"output/{config['project_name']}.{config['new_build']}.all_chrs.vcf.gz"

    shell:
        """
        bcftools concat -a {input.vcf} | bcftools sort  | bcftools view -Oz > {output}
        """

rule index_concat:
    input: f"output/{config['project_name']}.{config['new_build']}.all_chrs.vcf.gz"

    output: f"output/{config['project_name']}.{config['new_build']}.all_chrs.vcf.gz.csi"

    shell:
        """
        bcftools index {input}
        """
