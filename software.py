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
