# pipeline.py
import argparse
import os
import shutil
from software import *


def parameter_input():
    parser = argparse.ArgumentParser(description='Viral Genome Assembly Pipeline')
    parser.add_argument('input1', help='Path to read1 fastq(.gz)')
    parser.add_argument('input2', help='Path to read2 fastq(.gz)')
    parser.add_argument('--host',
                        help='Name(s) of Bowtie2 index(es) (default: None, skip host removal)',
                        default=None)
    parser.add_argument('-o', '--output', help='Output directory path',
                        default=os.path.join(os.getcwd(), 'result'))

    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads')
    parser.add_argument('-k', '--keep_log', action='store_true',
                        help='Resume interrupted run (do not modify output dir)',
                        default=False)
    parser.add_argument('--db', help='Path to database directory',
                        default="/cpfs01/projects-HDD/cfff-47998b01bebd_HDD/rj_24212030018/db")
    return parser.parse_args()


def main():
    args = parameter_input()
    threads = args.threads
    output = args.output
    db = args.db
    sample1 = get_sample_name(args.input1.split('/')[-1])
    sample2 = get_sample_name(args.input2.split('/')[-1])
    sample = sample1[0: -2]

    # 创建输出目录
    os.makedirs(os.path.join(output, "logs"), exist_ok=True)
    log_path = os.path.join(output, "logs", f"{sample}log.txt")
    if not args.keep_log:
        with open(log_path, "w") as f:
            f.write("0\n")

    with open(log_path, "r") as f:
        log = int(f.readline().strip())

    # 步骤1: fastp质控
    if log < 1:
        run_fastp(args.input1, args.input2, output, sample1, sample2, threads, sample)
        log = 1
        with open(log_path, "w") as f:
            f.write(f"{log}\n")

    # 步骤2: 去除人类污染
    if log < 2:
        if args.host:
            host_list = args.host.split(',')
            index_path = [os.path.join(db, "bowtie2_index", na, na) for na in host_list]
            run_bowtie2(output, threads, sample1, sample2, index_path, sample)
        else:
            print("Skipping host removal: No host index provided.")
        log = 2
        with open(log_path, "w") as f:
            f.write(f"{log}\n")

    # 步骤3: MEGAHIT组装
    if log < 3:
        run_megahit(args.output, threads, sample1, sample2, sample)
        log = 3
        with open(log_path, "w") as f:
            f.write(f"{log}\n")

    if log < 4:
        run_checkv(output, sample, db, threads)
        log = 4
        with open(log_path, "w") as f:
            f.write(f"{log}\n")

    if log < 5:
        run_dvf(output, sample, db, threads)
        log = 5
        with open(log_path, "w") as f:
            f.write(f"{log}\n")

    if log < 6:
        run_vibrant(output, sample, db, threads)
        log = 6
        with open(log_path, "w") as f:
            f.write(f"{log}\n")

    print("Pipeline completed successfully")


if __name__ == "__main__":
    main()
