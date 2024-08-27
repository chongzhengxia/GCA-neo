import os
import argparse
import sys

def receive_args():
    parser = argparse.ArgumentParser(description="Accepts fasta raw data absolute path, number of software threads, "
                                                 "sample_ID, library_info, lane_ID, platform_info")

    parser.add_argument('-R1_normal_path', '--R1_normal_path', type=str, nargs='+', required=True,
                        help="the absolute path of normal data R1")

    parser.add_argument('-R2_normal_path', '--R2_normal_path', type=str, nargs='+', required=True,
                        help="the absolute path of normal data R2")

    parser.add_argument('-R1_tumor_path', '--R1_tumor_path', type=str, nargs='+', required=True,
                        help="the absolute path of tumor data R1")

    parser.add_argument('-R2_tumor_path', '--R2_tumor_path', type=str, nargs='+', required=True,
                        help="the absolute path of tumor data R2")

    parser.add_argument('-sample_ID', '--sample_ID', type=str, nargs='+', required=True,
                        help="the sample ID ")

    parser.add_argument('-library_info', '--library_info', type=str, nargs='+', required=True,
                        help="the library information")

    parser.add_argument('-lane_ID', '--lane_ID', type=str, nargs='+', required=True,
                        help="the lane_ID")

    parser.add_argument('-platform_info', '--platform_info', type=str, nargs='+', required=True,
                        help="the platform information")

    parser.add_argument('-out_dir', '--out_dir', type=str, required=True,
                        help="the absolute path of output")

    parser.add_argument('-fastp_thread', '--fastp_thread', type=int, default=4,
                        help='fastp thread number, default:4, max:16')

    parser.add_argument('-bwa_thread', '--bwa_thread', type=int, default=4,
                        help='bwa thread number, default:4')

    parser.add_argument('-samtools_thread', '--samtools_thread', type=int, default=4,
                        help='samtools thread number, default:4')

    parser.add_argument('-mutect2_thread', '--mutect2_thread', type=int, default=4,
                        help='mutect2 thread number, default:4')

    # parser.add_argument('-bwa_ref', '--bwa_ref', type=str, default="/ref/bwa_index/hg38.fa",
    #                     help='the absolute path of bwa reference, default:/PATH/TO/THIS/FILE/ref/bwa_index/hg38.fa')
    # bwa索引在GCA-neo/ref/bwa_index目录下
    args = parser.parse_args()
    return args


def edit_name(R1, R2, TN_flag, sample_ID, library_info, platform_info, lane_ID):
    """

    [{'normal_sample1_xxxx_illuminal1_R1': 'normal_sample1_xxxx_illumina_l1_R1.gz', 'normal_sample1_xxxx_illuminal1_R2': 'normal_sample1_xxxx_illumina_l1_R2.gz'},
     {'normal_sample1_xxxx_illuminal2_R1': 'normal_sample1_xxxx_illumina_l2_R1.gz', 'normal_sample1_xxxx_illuminal2_R2': 'normal_sample1_xxxx_illumina_l2_R2.gz'},
     {'normal_sample1_xxxx_illuminal3_R1': 'normal_sample1_xxxx_illumina_l3_R1.gz', 'normal_sample1_xxxx_illuminal3_R2': 'normal_sample1_xxxx_illumina_l3_R2.gz'}]
    """

    if len(R1) == 1:
        name_path_list = []
        name_R1 = TN_flag + '_' + sample_ID[0] + '_' + library_info[0] + '_' + platform_info[0] + '_' + lane_ID[
            0] + '_R1'
        name_R2 = TN_flag + '_' + sample_ID[0] + '_' + library_info[0] + '_' + platform_info[0] + '_' + lane_ID[
            0] + '_R2'
        name_path_dict = {name_R1: R1[0], name_R2: R2[0]}
        name_path_list.append(name_path_dict)
        return name_path_list

    elif len(R1) > 1:
        name_path_list = []
        paired_data = list(zip(R1, R2))
        for index, (segmental_R1, segmental_R2) in enumerate(paired_data):
            name_R1 = TN_flag + '_' + sample_ID[index] + '_' + library_info[index] + '_' + platform_info[index] + '_' + \
                      lane_ID[index] + '_R1'
            name_R2 = TN_flag + '_' + sample_ID[index] + '_' + library_info[index] + '_' + platform_info[index] + '_' + \
                      lane_ID[index] + '_R2'
            name_path_dict = {name_R1: segmental_R1, name_R2: segmental_R2}
            name_path_list.append(name_path_dict)
        return name_path_list


def mkdir(out_dir, name_path_list):
    dir_demo1 = ['01_cleandata', '02_bwaout', '03_samtoolsout', '04_picardout', '05_markdupout', '06_BQSRout',
                 '07_mutect2out', '08_annovarout']

    # dir_demo2 = ['01_cleandata', '02_bwaout', '03_samtoolsout', '04_makedupout', '05_BQSRout', '06_mutect2out',
    #              '07_annovarout']

    out_dir = out_dir.rstrip('/')

    # if len(name_path_list) == 1:
    #     for i in dir_demo2:
    #         os.makedirs(out_dir + '/' + i)
    # elif len(name_path_list) > 1:
    #     for i in dir_demo1:
    #         os.makedirs(out_dir + '/' + i)
    for i in dir_demo1:
        os.makedirs(out_dir + '/' + i)


def quality_control(name_path_list, out_dir, fastp_thread):
    out_dir = out_dir.rstrip('/')
    command_demo1 = 'fastp -w {fastp_thread} -i ${input_R1} -I ${input_R2} -o {out_dir}/01_cleandata/{output_R1} -O {out_dir}/01_cleandata/{output_R2} 1>>{out_dir}/01_cleandata/rig_log.txt 2>>{out_dir}/01_cleandata/err_log.txt '
    for d in name_path_list:
        ((name_R1, path_R1), (name_R2, path_R2)) = d.items()
        exit_code = os.system(
            command_demo1.format(fastp_thread=fastp_thread, input_R1=path_R1, input_R2=path_R2, out_dir=out_dir,
                                 output_R1=name_R1, output_R2=name_R2))
        if exit_code != 0:
            sys.exit("An error occurred because fastp, please check log, log in {out_dir}/01_cleandata/".format(
                out_dir=out_dir))


def bwa_align(name_path_list, out_dir, bwa_thread):
    out_dir = out_dir.rstrip('/')
    command_demo1 = r'bwa mem -t {bwa_thread} -M -R "@RG\tID:${lane_ID}\tPL:${platform_info}\tLB:${library_info}\tSM:${sample_ID}" {bwa_ref} {input_R1} {input_R2} > {out_dir}/02_bwaout/{output_name_del_R1}.sam 1>>{out_dir/02_bwaout/rig_log.txt 2>>{out_dir}/02_bwaout/err_log.txt}'
    for d in name_path_list:
        ((name_R1, path_R1), (name_R2, path_R2)) = d.items()
        RG_info = name_R1.rstrip('_R1').split("_")
        [TN_flag, sample_ID, library_info, platform_info, lane_ID] = RG_info
        input_R1 = out_dir + '/01_cleandata/' + name_R1
        input_R2 = out_dir + '/01_cleandata/' + name_R2
        exit_code = os.system(command_demo1.format(bwa_thread=bwa_thread, lane_ID=lane_ID, platform_info=platform_info,
                                                   library_info=library_info, sample_ID=sample_ID,
                                                   bwa_ref=os.path.dirname(__file__) + '/ref/bwa_index/hg38.fa',
                                                   input_R1=input_R1, input_R2=input_R2, out_dir=out_dir,
                                                   output_name_del_R1=name_R1.rstrip('_R1')))

        if exit_code != 0:
            sys.exit(
                'An error occurred because bwa, please check log, log in {out_dir}/02_bwaout?'.format(out_dir=out_dir))


def samtools_process(name_path_list, out_dir, samtools_thread):
    out_dir = out_dir.rstrip('/')
    # view
    command_demo1 = 'samtools view -@ {samtools_thread} -b {input_sam} > {out_dir}/03_samtoolsout/{output_name_del_R1}.bam 1>>{out_dir}/03_samtoolsout/rig_view_log.txt 2>>{out_dir}/03_samtoolsout/err_view_log.txt'
    for d in name_path_list:
        ((name_R1, path_R1), (name_R2, path_R2)) = d.items()
        input_sam = out_dir + '/02_bwaout/' + name_R1.rstrip('_R1') + '.sam'
        exit_code_view = os.system(
            command_demo1.format(samtools_thread=samtools_thread, input_sam=input_sam, out_dir=out_dir,
                                 output_name_del_R1=name_R1.rstrip('_R1')))
        if exit_code_view != 0:
            sys.exit(
                'An error occurred because samtools view, please check log, log in {out_dir}/03_samtoolsout'.format(
                    out_dir=out_dir))
    # sort
    command_demo2 = 'samtools sort -@ {samtools_thread} -O bam -o {out_dir}/03_samtoolsout/{output_name_del_R1}.sort.bam {input_bam} 1>>{out_dir}/03_samtoolsout/rig_sort_log.txt 2>>{out_dir}/03_samtoolsout/err_sort_log.txt'
    for d in name_path_list:
        ((name_R1, path_R1), (name_R2, path_R2)) = d.items()
        input_bam = out_dir + '/03_samtoolsout/' + name_R1.rstrip('_R1') + '.bam'
        exit_code_sort = os.system(
            command_demo2.format(samtools_thread=samtools_thread, input_bam=input_bam, out_dir=out_dir,
                                 output_name=name_R1.rstrip('_R1')))
        if exit_code_sort != 0:
            sys.exit(
                'An error occurred because samtools sort, please check log, log in {out_dir}/03_samtoolsout'.format(
                    out_dir=out_dir))


def picard_process(name_path_list, out_dir):
    out_dir = out_dir.rstrip('/')
    input_bam_list = []
    command = 'picard MergeSamFiles {picard_command} O={out_dir}/04_picardout/{output_name_NT_flag_samplename}.all.sort.bam 1>{out_dir}/04_picardout/rig_log.txt 2>{out_dir}/04_picardout/err_log.txt}'
    for d in name_path_list:
        ((name_R1, path_R1), (name_R2, path_R2)) = d.items()
        input_bam_list.append(name_R1.rstrip('_R1') + '.sort.bam')
    input_bam_path_list = ['I=' + out_dir + '/03_samtoolsout/' + i for i in input_bam_list]
    picard_command = ' '.join(input_bam_path_list)
    exit_code = os.system(command.format(picard_command=picard_command, out_dir=out_dir,
                                         output_name_NT_flag_samplename='_'.join(name_R1.split('_')[0:2])))
    merged_bam_path_name = (
    out_dir + '/04_picardout/' + '_'.join(name_R1.split('_')[0:2]) + '.all.sort.bam', '_'.join(name_R1.split('_')[0:2]))
    if exit_code != 0:
        sys.exit(
            'An error occurred because picard, please check log, log in {out_dir}/04_picardout'.format(out_dir=out_dir))
    return merged_bam_path_name


def markduplicates_process(merged_bam_path_name, out_dir):
    """

    :param merged_bam_path_name:
    :param out_dir:
    :return: (标记去重后的全路径 + 包含后缀的文件名, 不包含后缀的文件名)
    """
    out_dir = out_dir.rstrip('/')
    command = 'gatk MarkDuplicates -I {input_bam} -O {out_dir}/05_markdupout/{output_bam_name}.markdup.bam -M {out_dir}/05_markdupout/{output_bam_name}.markdup.txt 1>>{out_dir}/05_markdupout/rig_log.txt 2>>{out_dir}/05_markdupout/err_log.txt'
    exit_code = os.system(
        command.format(input_bam=merged_bam_path_name[0], out_dir=out_dir, output_bam_name=merged_bam_path_name[1]))
    if exit_code != 0:
        sys.exit('An error occurred because markduplicates, please check log, log in {out_dir}/05_markdupout'.format(
            out_dir=out_dir))
    markdup_bam_path_name = (out_dir + '/05_markdupout/' + merged_bam_path_name[1] + '.markdup.bam',
                             merged_bam_path_name[1])  # 路径加名字拼接，作为下一个流程的输入文件
    return markdup_bam_path_name


def BQSR_process(markdup_bam_path_name, out_dir):
    """

    :param markdup_bam_path_name:
    :param out_dir:
    :return: (碱基质控后的全路径 + 包含后缀的文件名, 不包含后缀的文件名)
    """
    out_dir = out_dir.rstrip('/')
    command1 = r' gatk  BaseRecalibrator \
        -R {UCSC_hg38_ref} \
        -I {input_bam} \
        -O {out_dir}/06_BQSRout/{output_name}.recal_data.table \
        --known-sites {dbsnp_146_hg38_vcf_gz} \
        --known-sites {Mills_and_1000G_gold_standard_indels_hg38_vcf_gz} 1>>{out_dir}/06_BQSRout/rig_baserecalibrator_log.txt 2>>{out_dir}/06_BQSRout/err_baserecalibrator_log.txt'
    exit_code1 = os.system(command1.format(UCSC_hg38_ref=os.path.dirname(__file__) + '/ref/ucsc_hg38/hg38.fa',
                                           input_bam=markdup_bam_path_name[0], output_name=markdup_bam_path_name[1],
                                           dbsnp_146_hg38_vcf_gz=os.path.dirname(
                                               __file__) + '/ref/gatk/dbsnp_146.hg38.vcf.gz',
                                           Mills_and_1000G_gold_standard_indels_hg38_vcf_gz=os.path.dirname(
                                               __file__) + '/ref/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
                                           out_dir=out_dir))
    if exit_code1 != 0:
        sys.exit('An error occurred because BaseRecalibrator, please check log, log in {out_dir}/06_BQSRout'.format(
            out_dir=out_dir))
    command2 = r'gatk  ApplyBQSR \
        -R {UCSC_hg38_ref} \
        --bqsr-recal-file {out_dir}/06_BQSRout/{output_name}.recal_data.table \
        -I {input_bam} \
        -O {out_dir}/06_BQSRout/{output_name}.BQSR.bam 1>>{out_dir}/06_BQSRout/rig_ApplyBQSR_log.txt 2>>{out_dir}/06_BQSRout/err_ApplyBQSR_log.txt'
    exit_code2 = os.system(command2.format(UCSC_hg38_ref=os.path.dirname(__file__) + '/ref/ucsc_hg38/hg38.fa',
                                           input_bam=markdup_bam_path_name[0], output_name=markdup_bam_path_name[1],
                                           out_dir=out_dir))
    if exit_code2 != 0:
        sys.exit('An error occurred because ApplyBQSR, please check log, log in {out_dir}/06_BQSRout'.format(
            out_dir=out_dir))
    BQSR_bam_path_name = (out_dir + '/06_BQSRout/' + markdup_bam_path_name[1] + '.BQSR.bam', markdup_bam_path_name[1])
    return BQSR_bam_path_name


def mutect2_process(normal_BQSR_bam_path_name, tumor_BQSR_bam_path_name, out_dir, mutect2_thread):
    out_dir = out_dir.rstrip('/')
    command1 = r'''
    gatk Mutect2 --native-pair-hmm-threads {mutect2_thread} \
        -R {UCSC_hg38_ref} \
        -I {normal_input_bam} -normal {normal_name} \
        -I {tumor_input_bam} -tumor {tumor_name} \
        -O {out_dir}/07_mutect2out/{output_name}_mutect2.vcf 1>>{out_dir}/07_mutect2out/rig_Mutect2_log.txt 2>>{out_dir}/07_mutect2out/err_Mutect2_log.txt
        '''
    exit_code1 = os.system(command1.format(mutect2_thread=mutect2_thread,
                                           UCSC_hg38_ref=os.path.dirname(__file__) + '/ref/ucsc_hg38/hg38.fa',
                                           normal_input_bam=normal_BQSR_bam_path_name[0],
                                           normal_name=normal_BQSR_bam_path_name[1],
                                           tumor_input_bam=tumor_BQSR_bam_path_name[0],
                                           tumor_name=tumor_BQSR_bam_path_name[1],
                                           out_dir=out_dir, output_name=normal_BQSR_bam_path_name[0].split("_")[1]
                                           ))
    if exit_code1 != 0:
        sys.exit('An error occurred because Mutect2, please check log, log in {out_dir}/06_mutect2out'.format(
            out_dir=out_dir))
    command2 = r'''
    gatk FilterMutectCalls \
        -R {UCSC_hg38_ref} \
        -V {input_vcf} \
        -O {out_dir}/07_mutect2out/{output_name}_filter.vcf 
        cat {input_filter_vcf} | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }' > {out_dir}/07_mutect2out/{output_name}_final.vcf 1>>{out_dir}/07_mutect2out/rig_FilterMutectCalls_log.txt 2>>{out_dir}/07_mutect2out/err_FilterMutectCalls_log.txt
        '''
    input_vcf = out_dir + '/07_mutect2out/' + normal_BQSR_bam_path_name[0].split("_")[1] + '_mutect2.vcf'
    input_filter_vcf = out_dir + '/07_mutect2out/' + normal_BQSR_bam_path_name[0].split("_")[1] + '_filter.vcf'
    exit_code2 = os.system(command2.format(UCSC_hg38_ref=os.path.dirname(__file__) + '/ref/ucsc_hg38/hg38.fa',
                                           input_vcf=input_vcf, out_dir=out_dir,
                                           output_name=normal_BQSR_bam_path_name[0].split("_")[1],
                                           input_filter_vcf=input_filter_vcf))
    if exit_code2 != 0:
        sys.exit('An error occurred because FilterMutectCalls, please check log, log in {out_dir}/06_mutect2out'.format(
            out_dir=out_dir))
    mutect2_vcf_path_name = (out_dir + '/07_mutect2out/' + normal_BQSR_bam_path_name[0].split("_")[1] + '_final.vcf', normal_BQSR_bam_path_name[0].split("_")[1])
    return mutect2_vcf_path_name

def annovar_process(mutect2_vcf_path_name, out_dir, annovar_thread):
    out_dir = out_dir.rstrip('/')
    command = r"""
    perl {path_to_annovar}/table_annovar.pl \
        {input_vcf} \
        {path_to_annovar_db} \
        -buildver hg38 \
        -outfile {out_dir}/08_annovarout/{output_name} \
        -remove \
        -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a \
        -operation g,r,f,f,f \
        -nastring . \
        -vcfinput \
        -thread {annovar_thread} 1>>{out_dir}/07_annovarout/rig_log.txt 2>>{out_dir}/07_annovarout/err_log.txt
    """
    exit_code = os.system(command.format(path_to_annovar=os.path.dirname(__file__) + '/software/annovar',
                             input_vcf=mutect2_vcf_path_name[0],
                             path_to_annovar_db=os.path.abspath(__file__) + '/ref/annovar/humandb/',
                             out_dir=out_dir, output_name=mutect2_vcf_path_name[1], annovar_thread=annovar_thread))
    if exit_code != 0:
        sys.exit('An error occurred because annovar, please check log, log in {out_dir}/07_annovarout'.format(
            out_dir=out_dir))


def main():
    args = receive_args()
    normal_name_path_list = edit_name(args.R1_normal_path, args.R2_normal_path, 'normal', args.sample_ID,
                                      args.library_info, args.platform_info, args.lane_ID)
    tumor_name_path_list = edit_name(args.R1_tumor_path, args.R2_tumor_path, 'tumor', args.sample_ID, args.library_info,
                                     args.platform_info, args.lane_ID)
    mkdir(args.out_dir, normal_name_path_list)
    quality_control(normal_name_path_list, args.out_dir, args.fastp_thread)
    quality_control(tumor_name_path_list, args.out_dir, args.fastp_thread)
    bwa_align(normal_name_path_list, args.out_dir, args.bwa_thread)
    bwa_align(tumor_name_path_list, args.out_dir, args.bwa_thread)
    samtools_process(normal_name_path_list, args.out_dir, args.samtools_thread)
    samtools_process(tumor_name_path_list, args.out_dir, args.samtools_thread)
    normal_merged_bam_path_name = picard_process(normal_name_path_list, args.out_dir)
    tumor_merged_bam_path_name = picard_process(tumor_name_path_list, args.out_dir)
    normal_markdup_bam_path_name = markduplicates_process(normal_merged_bam_path_name, args.out_dir)
    tumor_markdup_bam_path_name = markduplicates_process(tumor_merged_bam_path_name, args.out_dir)
    normal_BQSR_bam_path_name = BQSR_process(normal_markdup_bam_path_name, args.out_dir)
    tumor_BQSR_bam_path_name = BQSR_process(tumor_markdup_bam_path_name, args.out_dir)
    mutect2_vcf_path_name = mutect2_process(normal_BQSR_bam_path_name, tumor_BQSR_bam_path_name, args.out_dir, args.mutect2_thread)
    annovar_process(mutect2_vcf_path_name, args.out_dir, args.annovar_thread)


if __name__ == '__main__':
    main()
