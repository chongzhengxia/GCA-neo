import argparse
import os
import sys


def receive_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('R1_path', type=str)
    parser.add_argument('R2_path', type=str)
    parser.add_argument('out_dir', type=str)
    parser.add_argument('sample_name', type=str)
    parser.add_argument('--star_thread', type=int, default=4)
    parser.add_argument('--fastp_thread', type=int, default=4)
    parser.add_argument('--star_genomeDir', type=str, default=os.path.dirname(__file__) + '/ref/STAR_index')
    parser.add_argument('--ucsc_genomeDir', type=str, default=os.path.dirname(__file__) + '/ref/ucsc_hg38/hg38.fa')
    args = parser.parse_args()
    return args


def check_args(args):
    if (not os.path.exists(args.R1_path)) or (not os.path.exists(args.R2_path)):
        sys.exit('there is no fastq data')


def mkdir_process(args):
    out_dir = os.path.abspath(args.out_dir)
    command = ['/01_fastpout', '/02_starout', '/03_sambambaout', '/04_splitncigarreadsout', '/05_BQSRout', '/06_haplotypecaller']
    for i in command:
        os.makedirs(os.path.abspath(out_dir) + i)


def fastp_process(args):
    out_dir = os.path.abspath(args.out_dir)
    command = 'fastp -w {fastp_thread} -i {R1_path} -I {R2_path} -o {out_dir}/01_fastpout/{output_R1} -O {out_dir}/01_fastpout/{output_R2} 2>>{out_dir}/01_fastpout/fastp_err.log '
    exit_code = os.system(
        command.format(fastp_thread=args.fastp_thread, R1_path=args.R1_path, R2_path=args.R2_path, out_dir=out_dir,
                       output_R1=os.path.basename(args.R1_path), output_R2=os.path.basename(args.R2_path)))
    if exit_code != 0:
        sys.exit('an error has occurred, log in {out_dir}/01_fastpout/fastp_err.log'.format(out_dir=out_dir))


def star_two_pass(args):
    out_dir = os.path.abspath(args.out_dir)
    command = """STAR --outFileNamePrefix {out_dir}/02_starout/{sample_name}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 1 \
    --outSAMstrandField intronMotif \
    --genomeDir {genomeDir} \
    --runThreadN {star_thread} \
    --readFilesCommand zcat  \
    --readFilesIn  {out_dir}/01_fastpout/{output_R1}  {out_dir}/01_fastpout/{output_R2} \
    --twopassMode Basic 2>>{out_dir}/02_starout/star_err.log"""
    exit_code = os.system(
        command.format(star_thread=args.star_thread, genomeDir=args.star_genomeDir,
                       out_dir=out_dir, sample_name=args.sample_name,
                       output_R1=os.path.basename(args.R1_path), output_R2=os.path.basename(args.R2_path)))
    if exit_code != 0:
        sys.exit('an error has occurred, log in {out_dir}/02_bwaout/star_err.log'.format(out_dir=out_dir))


def markdup(args):
    out_dir = os.path.abspath(args.out_dir)
    command = """sambamba markdup --overflow-list-size 600000 \
    --tmpdir='./' \
    -r {out_dir}/02_starout/{sample_name}_Aligned.sortedByCoord.out.bam {out_dir}/03_sambambaout/{sample_name}_rmPCR.bam \
    2>>{out_dir}/03_sambambaout/sambamba_err.log"""
    exit_code = os.system(command.format(out_dir=out_dir, sample_name=args.sample_name))
    if exit_code != 0:
        sys.exit('an error has occurred, log in {out_dir}/03_sambambaout/sambamba_err.log'.format(out_dir=out_dir))

def splitncigarreads(args):
    out_dir = os.path.abspath(args.out_dir)
    command = """gatk SplitNCigarReads -R {genomeDir} \
    -I {out_dir}/03_sambambaout/{sample_name}_rmPCR.bam \
    -O {out_dir}/04_splitncigarreadsout/{sample_name}_rmPCR_split.bam 2>>{out_dir}/04_splitncigarreadsout/splitncigarreads_err.log"""
    exit_code = os.system(command.format(genomeDir=args.ucsc_genomeDir, out_dir=out_dir,sample_name=args.sample_name))
    if exit_code != 0:
        sys.exit('an error has occurred, log in {out_dir}/04_splitncigarreadsout/splitncigarreads_err.log'.format(out_dir=out_dir))

def BQSR(args):
    out_dir = os.path.abspath(args.out_dir)
    command1 = """gatk AddOrReplaceReadGroups   \
    -I  {out_dir}/04_splitncigarreadsout/{sample_name}_rmPCR_split.bam  \
    -O  {out_dir}/05_BQSRout/{sample_name}_rmPCR_split_add.bam \
    -LB {sample_name} -PL ILLUMINA \
    -PU {sample_name} -SM {sample_name} 2>>{out_dir}/05_BQSRout/AddOrReplaceReadGroups_err.log"""
    exit_code1 = os.system(command1.format(out_dir=out_dir, sample_name=args.sample_name))
    if exit_code1 != 0:
        sys.exit('an error has occurred, log in {out_dir}/05_BQSRout/AddOrReplaceReadGroups_err.log'.format(out_dir=out_dir))
    command2 = """gatk BaseRecalibrator \
        -I {out_dir}/05_BQSRout/{sample_name}_rmPCR_split_add.bam \
        -R {genomeDir} \
        -O {out_dir}/05_BQSRout/{sample_name}_recal.table \
        --known-sites {dbsnp_146_hg38_vcf_gz} \
        --known-sites {Mills_and_1000G_gold_standard_indels_hg38_vcf_gz} 2>>{out_dir}/05_BQSRout/BaseRecalibrator_err.log"""
    exit_code2 = os.system(command2.format(out_dir=out_dir, sample_name=args.sample_name, genomeDir=args.ucsc_genomeDir,
                                           dbsnp_146_hg38_vcf_gz=os.path.dirname(__file__) + '/ref/gatk/dbsnp_146.hg38.vcf.gz',
                                           Mills_and_1000G_gold_standard_indels_hg38_vcf_gz=os.path.dirname(__file__) + '/ref/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'))
    if exit_code2 != 0:
        sys.exit('an error has occurred, log in {out_dir}/05_BQSRout/BaseRecalibrator_err.log'.format(out_dir=out_dir))

    command3 = """gatk ApplyBQSR \
        -I {out_dir}/05_BQSRout/{sample_name}_rmPCR_split_add.bam \
        -R {genomeDir} \
        -bqsr {out_dir}/05_BQSRout/{sample_name}_recal.table \
        --output {out_dir}/05_BQSRout/{sample_name}_recal.bam 2>>{out_dir}/05_BQSRout/ApplyBQSR_err.log"""
    exit_code3 = os.system(command3.format(out_dir=out_dir, sample_name=args.sample_name, genomeDir=args.ucsc_genomeDir))
    if exit_code3 != 0:
        sys.exit('an error has occurred, log in {out_dir}/05_BQSRout/ApplyBQSR_err.log'.format(out_dir=out_dir))


def haplotypecaller(args):
    out_dir = os.path.abspath(args.out_dir)
    command = """gatk  HaplotypeCaller  \
        -ERC GVCF --dbsnp {dbsnp_146_hg38_vcf_gz}  -R {genomeDir} \
        -I {out_dir}/05_BQSRout/{sample_name}_recal.bam  \
        -O  {out_dir}/06_haplotypecaller/{sample_name}.gvcf 2>>{out_dir}/06_haplotypecaller/haplotypecaller_err.log"""
    exit_code = os.system(command.format(out_dir=out_dir, sample_name=args.sample_name, genomeDir=args.ucsc_genomeDir,
                                         dbsnp_146_hg38_vcf_gz=os.path.dirname(__file__) + '/ref/gatk/dbsnp_146.hg38.vcf.gz'))
    if exit_code != 0:
        sys.exit('an error has occurred, log in {out_dir}/06_haplotypecaller/haplotypecaller_err.log'.format(out_dir=out_dir))




def main():
    args = receive_args()
    check_args(args)
    mkdir_process(args)
    fastp_process(args)
    star_two_pass(args)
    markdup(args)
    splitncigarreads(args)
    BQSR(args)
    haplotypecaller(args)


if __name__ == '__main__':
    main()
