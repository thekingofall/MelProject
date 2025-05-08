#!/usr/bin/env python3
import os
import gzip
import argparse
import subprocess
import shutil # Added import for shutil.which
from Bio import SeqIO
# from Bio.Seq import Seq # No longer needed as SeqRecord is not manually constructed
# from Bio.SeqRecord import SeqRecord # No longer needed
from tqdm import tqdm

# --- Default Configuration ---
DEFAULT_SEQ_C = "ATGTCGGAACTGTTGCTTGTCCGACT"
DEFAULT_SEQ_D = "AGTCGGACAAGCAACAGTTCCGACAT"
DEFAULT_OUTPUT_DIR = "Run1_fastq"
DEFAULT_NUM_CPUS = 10
# --- End Default Configuration ---

# --- Important Note on Conda Environment ---
# This script is designed to be run AFTER you have activated the 'hicpro3' conda environment.
# For example:
# conda activate hicpro3
# python your_script_name.py -i <input_dir> --num_cpus <N>
#
# The script will then attempt to run 'run_hicpro' at the end.
# ParaFly should also be available in this environment.
# ---

def find_fastq_pairs(input_dir):
    """
    Finds paired-end FASTQ files.
    """
    pairs = []
    all_gz_files = sorted([f for f in os.listdir(input_dir) if f.endswith(".fq.gz")])
    potential_r1_files = [f for f in all_gz_files if f.endswith("_1.fq.gz")]

    for r1_name in potential_r1_files:
        r2_name_candidate = r1_name[:-len("_1.fq.gz")] + "_2.fq.gz"
        if r2_name_candidate in all_gz_files:
            r1_path = os.path.abspath(os.path.join(input_dir, r1_name))
            r2_path = os.path.abspath(os.path.join(input_dir, r2_name_candidate))
            pairs.append((r1_path, r2_path))
        else:
            print(f"警告：未找到R1文件 '{r1_name}' 对应的R2文件 '{r2_name_candidate}'。将跳过此R1文件。") # Warning: R2 file not found for R1 file. Skipping.

    return pairs

def core_process_single_fastq_pair(r1_path, r2_path, output_dir, seq_c_target, seq_d_target):
    """
    Core logic to process a single pair of FASTQ files.
    Keeps pairs if SEQ_C or SEQ_D is found in either read. Sequences are NOT removed.
    """
    base_r1_name = os.path.basename(r1_path)
    base_r2_name = os.path.basename(r2_path)
    output_r1_path = os.path.join(output_dir, base_r1_name)
    output_r2_path = os.path.join(output_dir, base_r2_name)

    worker_pid = os.getpid()
    # print(f"[Worker PID: {worker_pid}] 开始处理文件对: {base_r1_name}, {base_r2_name}") # Starting processing pair
    # print(f"[Worker PID: {worker_pid}] 输出到: {output_r1_path}, {output_r2_path}") # Outputting to

    reads_processed_count = 0
    reads_kept_count = 0
    reads_discarded_empty_count = 0 # Counts pairs discarded because one or both reads were empty

    try:
        with gzip.open(r1_path, "rt") as r1_in_handle, \
             gzip.open(r2_path, "rt") as r2_in_handle, \
             gzip.open(output_r1_path, "wt") as r1_out_handle, \
             gzip.open(output_r2_path, "wt") as r2_out_handle:

            r1_parser = SeqIO.parse(r1_in_handle, "fastq")
            r2_parser = SeqIO.parse(r2_in_handle, "fastq")

            num_lines_r1 = 0
            # Count lines for progress bar more efficiently if possible, or estimate
            # For simplicity, keeping the original line counting method for now.
            with gzip.open(r1_path, "rt") as temp_r1_check:
                for _ in temp_r1_check:
                    num_lines_r1 += 1
            total_reads_in_pair = num_lines_r1 // 4

            paired_reads_iterator = zip(r1_parser, r2_parser)
            progress_bar_desc = f"[Worker PID: {worker_pid}] 过滤中 {base_r1_name[:15]}... & {base_r2_name[:15]}..." # Filtering

            for read1, read2 in tqdm(paired_reads_iterator, total=total_reads_in_pair, desc=progress_bar_desc, unit=" 对", leave=False, position=worker_pid % 10 if worker_pid else 0): # "pairs"
                reads_processed_count += 1

                seq1_str_orig = str(read1.seq)
                seq2_str_orig = str(read2.seq)

                s1_has_target = seq_c_target in seq1_str_orig or seq_d_target in seq1_str_orig
                s2_has_target = seq_c_target in seq2_str_orig or seq_d_target in seq2_str_orig

                keep_pair = s1_has_target or s2_has_target

                if keep_pair:
                    # Sequences are NOT removed. Original reads are written if the pair is kept.
                    if len(read1.seq) > 0 and len(read2.seq) > 0:
                        SeqIO.write(read1, r1_out_handle, "fastq")
                        SeqIO.write(read2, r2_out_handle, "fastq")
                        reads_kept_count += 1
                    else:
                        # This case means one of the original reads in the pair was empty.
                        reads_discarded_empty_count +=1

        print(f"[Worker PID: {worker_pid}] 完成处理 {base_r1_name} & {base_r2_name}. "
              f"总共处理: {reads_processed_count} 对, 保留: {reads_kept_count} 对, 因原始序列为空丢弃: {reads_discarded_empty_count} 对.")
        return True

    except Exception as e:
        print(f"[Worker PID: {worker_pid}] 处理 {base_r1_name} 和 {base_r2_name} 时发生错误: {e}") # Error processing
        return False

def run_parafly(command_file_path, num_cpus, log_dir):
    """
    Runs ParaFly with the generated command file.
    """
    failed_cmds_log = os.path.join(log_dir, "parafly_failed_cmds.log")
    parafly_command = [
        "ParaFly",
        "-c", command_file_path,
        "-CPU", str(num_cpus),
        "-verbose",
        "-failed_cmds", failed_cmds_log
    ]
    print(f"\n正在执行 ParaFly: {' '.join(parafly_command)}") # Executing ParaFly
    try:
        # Using subprocess.run for ParaFly itself for better control and output capture
        result = subprocess.run(parafly_command, capture_output=True, text=True, check=False)

        print("ParaFly STDOUT:")
        print(result.stdout)
        if result.stderr:
            print("ParaFly STDERR:")
            print(result.stderr)

        if result.returncode == 0:
            print("ParaFly 主进程成功完成。") # ParaFly main process completed successfully.
            if os.path.exists(failed_cmds_log) and os.path.getsize(failed_cmds_log) > 0:
                print(f"警告: ParaFly 报告部分命令执行失败。请检查日志: {failed_cmds_log}") # Warning: ParaFly reported some commands failed. Check log.
                return False
            return True
        else:
            print(f"ParaFly 主进程失败，退出码: {result.returncode}.") # ParaFly main process failed, exit code.
            if os.path.exists(failed_cmds_log) and os.path.getsize(failed_cmds_log) > 0:
                 print(f"ParaFly 报告了命令执行失败。请检查日志: {failed_cmds_log}") # ParaFly reported command failures. Check log.
            return False
    except FileNotFoundError:
        print("错误: ParaFly 命令未找到。请确保它已安装并且在您的PATH中 (通常在hicpro3 conda环境中)。") # Error: ParaFly command not found.
        return False
    except Exception as e:
        print(f"运行 ParaFly 时发生未知错误: {e}") # Unknown error running ParaFly.
        return False

def main():
    parser = argparse.ArgumentParser(
        description="使用ParaFly并行过滤成对的FASTQ文件。如果read包含特定序列则保留该read对 (但不移除这些序列)。随后使用 os.system 执行run_hicpro命令。",
        epilog="重要提示: 请先激活 'hicpro3' conda 环境再运行此脚本。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-i", "--input_dir",
        help="包含gzipped FASTQ文件的输入目录 (例如: sample_L1_1.fq.gz, sample_L1_2.fq.gz)。 (编排器模式下必需)"
    )
    parser.add_argument(
        "-o", "--output_dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"用于保存过滤后的gzipped FASTQ文件以及ParaFly日志的输出目录。\n(默认: {DEFAULT_OUTPUT_DIR})"
    )
    parser.add_argument(
        "--seq_c", default=DEFAULT_SEQ_C,
        help=f"用于筛选的第一个目标序列。\n(默认: {DEFAULT_SEQ_C})"
    )
    parser.add_argument(
        "--seq_d", default=DEFAULT_SEQ_D,
        help=f"用于筛选的第二个目标序列。\n(默认: {DEFAULT_SEQ_D})"
    )
    parser.add_argument(
        "--num_cpus", type=int,
        default=DEFAULT_NUM_CPUS,
        help=f"供ParaFly使用的CPU核心数量。\n(默认: {DEFAULT_NUM_CPUS})"
    )
    # Suppress help for worker arguments
    parser.add_argument("--run_worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--worker_r1", help=argparse.SUPPRESS)
    parser.add_argument("--worker_r2", help=argparse.SUPPRESS)

    args = parser.parse_args()

    if args.run_worker:
        if not (args.worker_r1 and args.worker_r2 and args.output_dir and args.seq_c and args.seq_d):
            # This error message is for the worker, so it's okay if it's brief
            print("[Worker Error] Missing necessary parameters (worker_r1, worker_r2, output_dir, seq_c, seq_d).")
            exit(1)

        if not os.path.exists(args.output_dir):
            try:
                os.makedirs(args.output_dir, exist_ok=True)
            except OSError as e:
                print(f"[Worker Error] Cannot create output directory '{args.output_dir}': {e}")
                exit(1)

        success = core_process_single_fastq_pair(
            args.worker_r1, args.worker_r2,
            args.output_dir, args.seq_c, args.seq_d
        )
        exit(0 if success else 1)

    # Orchestrator mode from here
    print("--- 启动FASTQ过滤与HiC-Pro执行脚本 (Orchestrator模式) ---")
    print(f"重要提示: 请确保您已激活 'hicpro3' conda 环境，并且 ParaFly 和 run_hicpro 可用。")

    if not args.input_dir:
        parser.error("-i/--input_dir is required in orchestrator mode.")

    if not os.path.isdir(args.input_dir):
        print(f"错误: 输入目录 '{args.input_dir}' 未找到或不是一个目录。")
        return

    absolute_output_dir = os.path.abspath(args.output_dir)
    if not os.path.exists(absolute_output_dir):
        try:
            os.makedirs(absolute_output_dir)
            print(f"已创建输出目录: {absolute_output_dir}")
        except OSError as e:
            print(f"错误: 无法创建输出目录 '{absolute_output_dir}': {e}")
            return
    elif not os.path.isdir(absolute_output_dir):
        print(f"错误: 输出路径 '{absolute_output_dir}' 已存在但不是一个目录。")
        return

    fastq_file_pairs = find_fastq_pairs(args.input_dir)
    if not fastq_file_pairs:
        print(f"在 '{args.input_dir}' 中未找到以 '_1.fq.gz' 和 '_2.fq.gz' 结尾的FASTQ文件对。")
        return

    print(f"\n找到 {len(fastq_file_pairs)} 对FASTQ文件准备通过ParaFly进行处理。")
    print(f"将根据以下序列进行筛选 (序列本身不会被移除): SEQ_C='{args.seq_c}', SEQ_D='{args.seq_d}'")
    print(f"处理后的文件将保存到: {absolute_output_dir}")
    print(f"ParaFly将使用 {args.num_cpus} 个CPU核心。")

    command_file_path = os.path.join(absolute_output_dir, "parafly_commands.txt")
    current_script_path = os.path.abspath(__file__) # Get absolute path of the current script

    with open(command_file_path, "w") as cmd_file:
        for r1_path, r2_path in fastq_file_pairs:
            # Ensure paths in the command file are properly quoted if they contain spaces
            command = (
                f"python \"{current_script_path}\" --run_worker "
                f"--worker_r1 \"{r1_path}\" --worker_r2 \"{r2_path}\" "
                f"-o \"{absolute_output_dir}\" "
                f"--seq_c \"{args.seq_c}\" --seq_d \"{args.seq_d}\""
            )
            cmd_file.write(command + "\n")

    print(f"\nParaFly 命令文件已生成: {command_file_path}")

    parafly_successful = run_parafly(command_file_path, args.num_cpus, absolute_output_dir)

    if not parafly_successful:
        print("\nParaFly 处理步骤失败或报告了错误。将跳过 run_hicpro 命令。")
        return

    print("\n所有FASTQ文件通过ParaFly处理完成。")

    print("\n准备使用 os.system 执行 run_hicpro 命令...")
    # IMPORTANT: Customize this command as needed for your HiC-Pro setup.
    # This is a placeholder. You will likely need to provide a config file and output dir for HiC-Pro.
    hicpro_command_to_run_list = ["run_hicpro"]
    # Example of a more complete HiC-Pro command:
    # Ensure these paths are correct and accessible.
    # hicpro_config_file = "/path/to/your/hicpro_config.txt" # Make this an argument or a fixed known path
    # hicpro_results_dir = os.path.join(os.path.dirname(absolute_output_dir), "hicpro_analysis_results") # Example output for HiC-Pro
    # if not os.path.exists(hicpro_results_dir):
    #     os.makedirs(hicpro_results_dir, exist_ok=True)
    #
    # hicpro_command_to_run_list = [
    #     "run_hicpro",
    #     "-c", hicpro_config_file,
    #     "-i", absolute_output_dir, # Input for HiC-Pro is the filtered FASTQ dir from this script
    #     "-o", hicpro_results_dir
    # ]

    command_name = hicpro_command_to_run_list[0]
    absolute_command_path = shutil.which(command_name)

    if absolute_command_path is None:
        print(f"错误: run_hicpro 的命令 '{command_name}' 未在系统PATH中找到。")
        print(f"请确保 '{command_name}' 已安装，并且正确的conda环境已激活，或者提供命令的完整路径。")
        print(f"当前系统PATH: {os.environ.get('PATH')}")
        return # Exit if command not found

    # Construct the command string for os.system
    # The first element is the command itself (now with absolute path), followed by its arguments
    full_hicpro_command_list = [absolute_command_path] + hicpro_command_to_run_list[1:]
    hicpro_command_string = ' '.join(f'"{arg}"' if ' ' in arg else arg for arg in full_hicpro_command_list) # Quote arguments with spaces

    print(f"将要执行的 run_hicpro 命令: {hicpro_command_string}")

    try:
        # Execute the command using os.system
        exit_status = os.system(hicpro_command_string)
        
        # os.system returns the exit status of the command.
        # On Unix, a 0 exit status typically means success.
        if exit_status == 0:
            print(f"run_hicpro 命令 '{hicpro_command_string}' 已成功启动/完成 (os.system 返回码: {exit_status})。")
            print("请检查 run_hicpro 的具体输出来确认其实际完成情况和结果。")
        else:
            # The return value is shell-dependent. For bash, 127 means command not found (though shutil.which should prevent this).
            # Other non-zero values indicate errors from the command itself.
            print(f"run_hicpro 命令 '{hicpro_command_string}' 执行失败或返回了非零退出码: {exit_status}。")
            print("请检查 run_hicpro 的输出以及相关的日志文件。")

    except Exception as e: # Catch any other unexpected errors during os.system call itself
        print(f"尝试执行 run_hicpro 命令时发生意外错误: {e}")


    print("\n--- 脚本执行完毕 ---") # Script execution finished.

if __name__ == "__main__":
    main()
