import csv
import subprocess
import sys
import os
import shutil


def get_sample_name(file):
    if file[-6:] == '.fq.gz' or file[-6:] == '.fastq':
        return file[0: -6]
    elif file[-3:] == '.fq':
        return file[0: -3]
    elif file[-9:] == '.fastq.gz':
        return file[0: -9]


def run_fastp(input1, input2, output, sample1, sample2, threads, sample):
    trim_dir = os.path.join(output, "1.trimmed")
    os.makedirs(trim_dir, exist_ok=True)

    # 设置fastp参数
    params = "-l 90 -q 20 -u 30 -y --trim_poly_g"

    cmd = f"""
    fastp --in1 {input1} --in2 {input2} \
          --out1 {trim_dir}/{sample1}.fq.gz \
          --out2 {trim_dir}/{sample2}.fq.gz \
          {params} \
          --thread {threads} \
          --html {trim_dir}/{sample}report.html \
          --json {trim_dir}/report.json
    """
    s = subprocess.call(cmd, shell=True)
    if s != 0:
        print("ERRO run_fastp!")


def run_bowtie2(output, threads, sample1, sample2, index_path, sample):
    bowtie2_dir = os.path.join(output, "2.bowtie2")
    sample_dir = os.path.join(bowtie2_dir, sample)
    trim_dir = os.path.join(output, "1.trimmed")

    if os.path.exists(sample_dir):
        shutil.rmtree(sample_dir)  # 修复：使用shutil.rmtree删除非空目录
    os.makedirs(sample_dir, exist_ok=True)

    try:
        shutil.copy(os.path.join(trim_dir, f"{sample1}.fq.gz"),
                    os.path.join(sample_dir, f"{sample1}.fq.gz"))
        shutil.copy(os.path.join(trim_dir, f"{sample2}.fq.gz"),
                    os.path.join(sample_dir, f"{sample2}.fq.gz"))
    except Exception as e:
        sys.exit(f"File copy failed: {e}")

    for index in index_path:
        input1 = os.path.join(sample_dir, f"{sample1}.fq.gz")
        input2 = os.path.join(sample_dir, f"{sample2}.fq.gz")
        tmp_prefix = os.path.join(sample_dir, "tmp")
        sam_output = os.path.join(sample_dir, "tmp.sam")
        cmd = [
            "bowtie2", "-p", str(threads),
            "-x", index,
            "-1", input1, "-2", input2,
            "--un-conc", tmp_prefix,
            "-S", sam_output
        ]
        try:
            subprocess.run(cmd, check=True)
            shutil.move(f"{tmp_prefix}.1", os.path.join(sample_dir, f"{sample1}.fastq"))
            shutil.move(f"{tmp_prefix}.2", os.path.join(sample_dir, f"{sample2}.fastq"))
        except Exception as e:
            sys.exit(f"bowtie2 failed: {e}")

    # 使用安全路径处理
    subprocess.run(["rm", "-rf", f"{output}/1.trimmed/{sample}*"], check=False)
    subprocess.run(["rm", "-rf", f"{output}/2.bowtie2/{sample}/tmp.sam"], check=False)


def run_megahit(output, threads, sample1, sample2, sample):
    assembly_dir = os.path.join(output, "3.assembly", sample)
    if os.path.exists(assembly_dir):
        shutil.rmtree(assembly_dir)

    # 设置k-mer参数
    k_list = "21,29,39,59,79,99,119"

    cmd = f"""
    megahit -1 {output}/2.bowtie2/{sample}/{sample1}.fastq \
            -2 {output}/2.bowtie2/{sample}/{sample2}.fastq \
            -o {assembly_dir} \
            --k-list {k_list} \
            --num-cpu-threads {threads} \
            --min-contig-len 1000
    """
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"MEGAHIT failed: {e}")


def run_checkv(output, sample, db, threads):
    checkv_dir = os.path.join(output, "4.checkv", sample)
    input_fasta = os.path.join(output, "3.assembly", sample, "final.contigs.fa")
    if os.path.exists(checkv_dir):
        shutil.rmtree(checkv_dir)
    os.makedirs(checkv_dir, exist_ok=True)

    # Step 1: 运行CheckV初步筛选 SRR10983056
    cmd = f"""
    checkv end_to_end {input_fasta} {checkv_dir} \
        -d {os.path.join(db, 'checkvdb/checkv-db-v1.0')} \
        -t {threads}
    """
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"CheckV failed: {e}")

    # Step 2: 解析CheckV结果进行过滤
    quality_summary_path = os.path.join(checkv_dir, 'quality_summary.tsv')
    if not os.path.exists(quality_summary_path):
        sys.exit(f"CheckV quality summary file not found: {quality_summary_path}")

    contigs_to_remove = []
    contigs_to_save = []
    with open(quality_summary_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                host_genes = int(row['host_genes'])
                viral_genes = int(row['viral_genes'])
                gene_count = int(row['gene_count'])
            except KeyError as e:
                sys.exit(f"Missing required column in CheckV quality summary: {e}")
            except ValueError as e:
                sys.exit(f"Invalid value in CheckV quality summary: {e}")
            # 应用过滤条件：原核基因超过总基因的50%且超过病毒基因10倍
            if (host_genes > 0.5 * gene_count) and (host_genes > 10 * viral_genes):
                contigs_to_remove.append(row['contig_id'])
            if (viral_genes > host_genes) and (viral_genes > 10):
                contigs_to_save.append(row['contig_id'])

    # 定义读取fasta文件的辅助函数
    def read_fasta(file_path):
        contigs = {}
        with open(file_path, 'r') as f:
            header = ''
            seq = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if header:
                        contigs[header] = ''.join(seq)
                        seq = []
                    header = line[1:].split()[0]  # 提取contig ID（第一个字段）
                else:
                    seq.append(line)
            if header:  # 添加最后一个contig
                contigs[header] = ''.join(seq)
        return contigs

    # 读取原始fasta并过滤
    contigs = read_fasta(input_fasta)
    filtered_contigs = {k: v for k, v in contigs.items() if k not in contigs_to_remove}
    save_contigs = {k: v for k, v in contigs.items() if k in contigs_to_save}

    # 写入过滤后的fasta文件
    output_fasta = os.path.join(checkv_dir, 'filtered_contigs.fa')
    with open(output_fasta, 'w') as f:
        for header, seq in filtered_contigs.items():
            f.write(f'>{header}\n{seq}\n')
    print(f"Filtered contigs saved to: {output_fasta}")

    output_fasta = os.path.join(checkv_dir, 'checkv.fa')
    with open(output_fasta, 'w') as f:
        for header, seq in save_contigs.items():
            f.write(f'>{header}\n{seq}\n')

    print(f"checkv contigs saved to: {output_fasta}")


def run_dvf(output, sample, db, threads):
    dvf_dir = os.path.join(output, "5.dvf", sample)
    input_fasta = os.path.join(output, "4.checkv", sample, "filtered_contigs.fa")
    if os.path.exists(dvf_dir):
        shutil.rmtree(dvf_dir)
    os.makedirs(dvf_dir, exist_ok=True)
    cmd = f"""
    /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/rj_24212030018/miniconda3/envs/viroprofiler-dvf/bin/python \
    /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/rj_24212030018/miniconda3/envs/viroprofiler-dvf/bin/dvf.py \
    -i {input_fasta} -o {dvf_dir} -c {threads}   \
    -m /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/rj_24212030018/miniconda3/envs/viroprofiler-dvf/share/deepvirfinder/models
"""

    s1 = subprocess.call(cmd, shell=True)
    if s1 != 0:
        sys.exit(f"DVF failed: {s1}")

    # # qvalue过滤结果
    # grep = f"""
    # /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/rj_24212030018/miniconda3/envs/viroprofiler-dvf/bin/Rscript \
    # /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/rj_24212030018/testflow/nf-core-nextvirus/bin/calc_qvalue.r \
    # {dvf_dir}/*_dvfpred.txt 0.9 {dvf_dir}/dvf_virus.tsv
    # """
    # s2 = subprocess.call(grep, shell=True)
    # if s2 != 0:
    #     print(f"grep failed: {s2}")
    # cmd = f"sed 1d {dvf_dir}/dvf_virus.tsv | cut -f1 > {dvf_dir}/virus_dvf.list"
    # s3 = subprocess.call(cmd, shell=True)
    # if s3 != 0:
    #     print(f"sed failed: {s3}")
    cmd = f"awk 'NR>1 && $3 > 0.9 && $4 < 0.01 {{print $1}}' {dvf_dir}/*_dvfpred.txt > {dvf_dir}/virus_dvf.list"
    s3 = subprocess.call(cmd, shell=True)
    if s3 != 0:
        print(f"sed failed: {s3}")

    cmd = f"seqkit grep -f {dvf_dir}/virus_dvf.list {input_fasta} > {dvf_dir}/dvf.fasta"
    s4 = subprocess.call(cmd, shell=True)
    if s4 != 0:
        print(f"seqkit failed: {s4}")


def run_vibrant(output, sample, db, threads):
    vibrant_dir = os.path.join(output, "6.vibrant", sample)
    input_fasta = os.path.join(output, "4.checkv", sample, "filtered_contigs.fa")
    if os.path.exists(vibrant_dir):
        shutil.rmtree(vibrant_dir)
    os.makedirs(vibrant_dir, exist_ok=True)

    cmd = f"""
        source activate /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/rj_24212030018/miniconda3/envs/vibrant &&
        VIBRANT_run.py -i {input_fasta} -folder {vibrant_dir}  -t {threads} \
        -d {db}/vibrant/databases  -m {db}/vibrant/files
    """

    s1 = subprocess.call(cmd, shell=True)
    if s1 != 0:
        sys.exit(f"DVF failed: {s1}")


def run_combine(output, sample, db, threads):
    # 定义输入文件路径
    checkv_file = os.path.join(output, "4.checkv", sample, "checkv.fa")
    dvf_file = os.path.join(output, "5.dvf", sample, "dvf.fasta")
    vibrant_file = os.path.join(output, "6.vibrant", sample,
                                "VIBRANT_filtered_contigs",
                                "VIBRANT_phages_filtered_contigs",
                                "filtered_contigs.phages_combined.fna")

    # 创建输出目录
    combine_dir = os.path.join(output, "7.combine", sample)
    combined_file = os.path.join(combine_dir, "virus_combine.fa")

    # 清空并重建目录
    if os.path.exists(combine_dir):
        shutil.rmtree(combine_dir)
    os.makedirs(combine_dir, exist_ok=True)

    # 合并文件
    with open(combined_file, 'wb') as f_out:
        # 依次合并三个文件
        for input_file in [checkv_file, dvf_file, vibrant_file]:
            if os.path.exists(input_file):
                with open(input_file, 'rb') as f_in:
                    shutil.copyfileobj(f_in, f_out)
            else:
                print(f"Warning: Missing input file {input_file}")

    print(f"Merged files saved to: {combined_file}")


def run_filter(output, sample, db, threads):
    filter_dir = os.path.join(output, "8.filter", sample)
    input_fasta = os.path.join(output, "7.combine", sample, "virus_combine.fa")
    if os.path.exists(filter_dir):
        shutil.rmtree(filter_dir)
    os.makedirs(filter_dir, exist_ok=True)

    # Step 1: 运行CheckV初步筛选
    cmd = f"""
    checkv end_to_end {input_fasta} {filter_dir} \
        -d {os.path.join(db, 'checkvdb/checkv-db-v1.0')} \
        -t {threads}
    """
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"CheckV failed: {e}")

    # Step 2: 解析CheckV结果进行过滤
    quality_summary_path = os.path.join(filter_dir, 'quality_summary.tsv')
    if not os.path.exists(quality_summary_path):
        sys.exit(f"CheckV quality summary file not found: {quality_summary_path}")

    contigs_to_remove = []
    with open(quality_summary_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                completeness_str = row['completeness'].strip()
                # 处理 NA 值（替换为 0）
                if completeness_str == 'NA':
                    completeness_str = '0'
                completeness = float(completeness_str)
            except KeyError as e:
                sys.exit(f"Missing required column in CheckV quality summary: {e}")
            except ValueError as e:
                sys.exit(f"Invalid value in CheckV quality summary: {e}")
            # 应用过滤条件：completeness score < 50%
            if completeness < 50:
                contigs_to_remove.append(row['contig_id'])

            # 定义读取fasta文件的辅助函数

        def read_fasta(file_path):
            contigs = {}
            with open(file_path, 'r') as f:
                header = ''
                seq = []
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if header:
                            contigs[header] = ''.join(seq)
                            seq = []
                        header = line[1:].split()[0]  # 提取contig ID（第一个字段）
                    else:
                        seq.append(line)
                if header:  # 添加最后一个contig
                    contigs[header] = ''.join(seq)
            return contigs

    # 读取原始fasta并过滤
    contigs = read_fasta(input_fasta)
    filtered_contigs = {k: v for k, v in contigs.items() if k not in contigs_to_remove}

    # 写入过滤后的fasta文件
    output_fasta = os.path.join(filter_dir, 'filtered_contigs.fa')
    with open(output_fasta, 'w') as f:
        for header, seq in filtered_contigs.items():
            f.write(f'>{header}\n{seq}\n')
    print(f"Filtered contigs saved to: {output_fasta}")


def run_busco_filter(output, sample, db, threads):
    busco_dir = os.path.join(output, "9.busco_filter", sample)
    input_fasta = os.path.join(output, "8.filter", sample, "filtered_contigs.fa")

    # 清理并创建目录
    if os.path.exists(busco_dir):
        shutil.rmtree(busco_dir)
    os.makedirs(busco_dir, exist_ok=True)

    # Step 1: 使用BUSCO
    cmd = f"""
    source activate /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/rj_24212030018/miniconda3/envs/busco &&
     busco -f -i {input_fasta} -c {threads} -o {busco_dir} -m geno -l {db}/bacteria_odb12 --offline
     """
    s1 = subprocess.call(cmd, shell=True)
    if s1 != 0:
        sys.exit(f"busco filter failed: {s1}")
    # Step 2: 解析基因预测结果统计总基因数
    predicted_file = os.path.join(busco_dir, r"prodigal_output/predicted_genes/predicted.fna")
    # Count genes per contig
    predicted_counts = {}
    with open(predicted_file, 'r') as pf:
        for line in pf:
            if line.startswith('>'):
                header = line[1:].split()[0]
                contig = '_'.join(header.split('_')[:-1])  # contig = part before first underscore
                predicted_counts[contig] = predicted_counts.get(contig, 0) + 1
    if not predicted_counts:
        sys.exit("Error: No predicted genes found in predicted.fna.")
    print(f"Total contigs with predicted genes: {len(predicted_counts)}")
    total_genes = sum(predicted_counts.values())
    print(f"Total predicted genes: {total_genes}")
    full_table = os.path.join(busco_dir, r"run_bacteria_odb12/full_table.tsv")
    busco_counts = {}
    with open(full_table, 'r') as ft:
        reader = csv.reader(ft, delimiter='\t')
        headers = next(reader, None)
        for row in reader:
            if len(row) < 3:
                continue
            status = row[1].strip()
            if status in ("Complete", "Fragmented"):
                seq_field = row[2]
                # If sequence field includes "file:contig:start-end", extract contig
                if ':' in seq_field:
                    parts = seq_field.split(':')
                    contig_name = parts[-2] if len(parts) >= 2 else parts[0]
                else:
                    contig_name = seq_field
                contig_name = contig_name.split()[0]
                busco_counts[contig_name] = busco_counts.get(contig_name, 0) + 1

    total_busco_hits = sum(busco_counts.values()) if busco_counts else 0
    print(f"Total BUSCO genes (Complete/Frag): {total_busco_hits}")
    print(f"Contigs with BUSCO hits: {len(busco_counts)}")
    contigs_to_remove = []
    for contig, gene_count in predicted_counts.items():
        if gene_count == 0:
            continue
        busco_genes = busco_counts.get(contig, 0)
        ratio = busco_genes / gene_count
        print(f"Contig {contig}: {busco_genes}/{gene_count} BUSCO genes (ratio {ratio:.2%})")
        if ratio > 0.05:
            contigs_to_remove.append(contig)
    if contigs_to_remove:
        print(f"Removing {len(contigs_to_remove)} contigs with BUSCO ratio > 5%: {contigs_to_remove}")
    else:
        print("No contigs exceed BUSCO ratio threshold (5%).")
    input_fasta = os.path.join(output, "8.filter", sample, "filtered_contigs.fa")
    output_fasta = os.path.join(busco_dir, "filtered_contigs.fa")
    # Read input FASTA sequences
    contig_seqs = {}
    with open(input_fasta, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    contig_seqs[header] = ''.join(seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            contig_seqs[header] = ''.join(seq_lines)

    # Write filtered sequences
    with open(output_fasta, 'w') as out:
        for hdr, seq in contig_seqs.items():
            if hdr not in contigs_to_remove:
                out.write(f">{hdr}\n{seq}\n")

    print(f"Filtered contigs saved to: {output_fasta}")
