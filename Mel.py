#!/usr/bin/env python3
import os
import gzip
import argparse
import subprocess
import shutil # Added import for shutil.which
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

def remove_substring_and_adjust_qual(record, substring_to_remove):
    """
    Removes all non-overlapping occurrences of substring_to_remove from a SeqRecord's sequence
    and adjusts its PHRED quality scores accordingly.
    """
    original_seq_str = str(record.seq)
    if not substring_to_remove:
        return record

    original_qualities = record.letter_annotations.get("phred_quality")
    
    if original_qualities is None or len(original_qualities) != len(original_seq_str):
        new_seq_str_only = original_seq_str.replace(substring_to_remove, "")
        return SeqRecord(Seq(new_seq_str_only),
                         id=record.id, name=record.name, description=record.description,
                         letter_annotations={}) 

    new_seq_parts = []
    new_qual_parts = []
    current_pos = 0
    modified = False

    while current_pos < len(original_seq_str):
        found_at = original_seq_str.find(substring_to_remove, current_pos)
        if found_at != -1:
            new_seq_parts.append(original_seq_str[current_pos:found_at])
            new_qual_parts.extend(original_qualities[current_pos:found_at])
            current_pos = found_at + len(substring_to_remove)
            modified = True
        else:
            new_seq_parts.append(original_seq_str[current_pos:])
            new_qual_parts.extend(original_qualities[current_pos:])
            break
    
    if not modified:
        return record

    final_seq_str = "".join(new_seq_parts)
    final_qualities = new_qual_parts

    if len(final_seq_str) != len(final_qualities):
        worker_pid = os.getpid()
        print(f"[Worker PID: {worker_pid}] Warning: Length mismatch after modifying read {record.id}. Seq len: {len(final_seq_str)}, Qual len: {len(final_qualities)}. Qualities will be cleared.")
        return SeqRecord(Seq(final_seq_str),
                         id=record.id, name=record.name, description=record.description,
                         letter_annotations={})

    return SeqRecord(Seq(final_seq_str),
                     id=record.id,
                     name=record.name,
                     description=record.description,
                     letter_annotations={"phred_quality": final_qualities})


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
    Keeps pairs if SEQ_C or SEQ_D is found in either read, and removes these sequences.
    """
    base_r1_name = os.path.basename(r1_path)
    base_r2_name = os.path.basename(r2_path)
    output_r1_path = os.path.join(output_dir, base_r1_name)
    output_r2_path = os.path.join(output_dir, base_r2_name)

    worker_pid = os.getpid()
    print(f"[Worker PID: {worker_pid}] 开始处理文件对: {base_r1_name}, {base_r2_name}") # Starting processing pair
    print(f"[Worker PID: {worker_pid}] 输出到: {output_r1_path}, {output_r2_path}") # Outputting to

    reads_processed_count = 0
    reads_kept_count = 0
    reads_became_empty_count = 0

    try:
        with gzip.open(r1_path, "rt") as r1_in_handle, \
             gzip.open(r2_path, "rt") as r2_in_handle, \
             gzip.open(output_r1_path, "wt") as r1_out_handle, \
             gzip.open(output_r2_path, "wt") as r2_out_handle:

            r1_parser = SeqIO.parse(r1_in_handle, "fastq")
            r2_parser = SeqIO.parse(r2_in_handle, "fastq")

            num_lines_r1 = 0
            with gzip.open(r1_path, "rt") as temp_r1_check:
                for _ in temp_r1_check:
                    num_lines_r1 += 1
            total_reads_in_pair = num_lines_r1 // 4
            
            paired_reads_iterator = zip(r1_parser, r2_parser)
            progress_bar_desc = f"[Worker PID: {worker_pid}] 过滤中 {base_r1_name[:10]}... & {base_r2_name[:10]}..." # Filtering
            
            for read1, read2 in tqdm(paired_reads_iterator, total=total_reads_in_pair, desc=progress_bar_desc, unit=" 对", leave=False, position=worker_pid % 10): # "pairs"
                reads_processed_count += 1
                
                seq1_str_orig = str(read1.seq)
                seq2_str_orig = str(read2.seq)

                s1_has_target = seq_c_target in seq1_str_orig or seq_d_target in seq1_str_orig
                s2_has_target = seq_c_target in seq2_str_orig or seq_d_target in seq2_str_orig

                keep_pair = s1_has_target or s2_has_target

                if keep_pair:
                    current_read1 = read1
                    current_read2 = read2

                    current_read1 = remove_substring_and_adjust_qual(current_read1, seq_c_target)
                    current_read1 = remove_substring_and_adjust_qual(current_read1, seq_d_target)
                    
                    current_read2 = remove_substring_and_adjust_qual(current_read2, seq_c_target)
                    current_read2 = remove_substring_and_adjust_qual(current_read2, seq_d_target)
                    
                    if len(current_read1.seq) > 0 and len(current_read2.seq) > 0:
                        SeqIO.write(current_read1, r1_out_handle, "fastq")
                        SeqIO.write(current_read2, r2_out_handle, "fastq")
                        reads_kept_count += 1
                    else:
                        reads_became_empty_count +=1
            
        print(f"[Worker PID: {worker_pid}] 完成处理 {base_r1_name} & {base_r2_name}. " # Finished processing
              f"总共处理: {reads_processed_count} 对, 保留(修改后): {reads_kept_count} 对, 修改后变空丢弃: {reads_became_empty_count} 对.") # Total processed, Kept (modified), Discarded (empty after modification)
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

def run_final_command(command_list_original, command_name="命令"): # command_name="Command"
    """
    Runs a final external shell command, trying to find the absolute path of the command.
    """
    if not command_list_original:
        print(f"错误：提供给 {command_name} 的命令列表为空。") # Error: Command list provided to command_name is empty.
        return False

    command_to_find = command_list_original[0]
    absolute_command_path = shutil.which(command_to_find)

    if absolute_command_path is None:
        print(f"错误: {command_name} 的命令 '{command_to_find}' 未在系统PATH中找到。") # Error: command_name's command not found in system PATH.
        print(f"请确保 '{command_to_find}' 已安装，并且正确的conda环境已激活，或者提供命令的完整路径。") # Please ensure it's installed and correct conda env is active, or provide full path.
        print(f"当前系统PATH: {os.environ.get('PATH')}") # Current system PATH
        return False
    
    # Use the absolute path for the command
    actual_command_list = [absolute_command_path] + command_list_original[1:]
    
    print(f"\n正在执行 {command_name} (使用路径 {absolute_command_path}): {' '.join(actual_command_list)}") # Executing command_name (using path ...)

    try:
        # It's good practice to pass the environment explicitly if issues persist,
        # but for now, rely on the environment inherited by the Python script.
        # env=os.environ.copy() # To pass the current environment
        result = subprocess.run(actual_command_list, check=True, capture_output=True, text=True) # Add env=env if needed
        
        print(f"{command_name} STDOUT:")
        print(result.stdout)
        if result.stderr: # Print stderr even on success, as some tools use it for info
            print(f"{command_name} STDERR (如果有):") # if any
            print(result.stderr)
        print(f"{command_name} '{' '.join(actual_command_list)}' 执行成功。") # executed successfully.
        return True
    except subprocess.CalledProcessError as e:
        print(f"{command_name} '{' '.join(actual_command_list)}' 执行失败。返回码: {e.returncode}") # failed. Return code.
        print(f"{command_name} STDOUT (如有):") # STDOUT (if any)
        if e.stdout:
            print(e.stdout)
        print(f"{command_name} STDERR (如有):") # STDERR (if any)
        if e.stderr:
            print(e.stderr)
        return False
    # FileNotFoundError should be less likely now that we use shutil.which, but keep for robustness
    except FileNotFoundError: 
        print(f"错误: {command_name} 的命令 '{actual_command_list[0]}' 在尝试执行时未找到（即使shutil.which找到了它，这很奇怪）。") # Error: command_name's command not found when trying to execute (strange if shutil.which found it).
        return False
    except PermissionError as e: # Catching PermissionError explicitly
        print(f"{command_name} '{' '.join(actual_command_list)}' 执行时遇到权限错误: {e}") # Permission error when executing.
        print(f"请检查用户 '{os.geteuid()}' 对 '{actual_command_list[0]}' 及其所需文件的执行/访问权限。") # Please check permissions for user on the command and its required files.
        return False
    except Exception as e:
        print(f"执行 {command_name} 时发生未知错误: {e}") # Unknown error executing command_name.
        return False

def main():
    parser = argparse.ArgumentParser(
        description="使用ParaFly并行过滤成对的FASTQ文件。如果read包含特定序列则保留该read对，并在写入前从中移除这些特定序列。随后执行run_hicpro命令。" # Filter paired FASTQ files in parallel using ParaFly...
                    "重要提示: 请先激活 'hicpro3' conda 环境再运行此脚本。", # Important: Activate 'hicpro3' conda environment first.
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-i", "--input_dir", 
        help="包含gzipped FASTQ文件的输入目录 (例如: sample_L1_1.fq.gz, sample_L1_2.fq.gz)。 (编排器模式下必需)" # Input directory... (Required in orchestrator mode)
    )
    parser.add_argument(
        "-o", "--output_dir", 
        default=DEFAULT_OUTPUT_DIR,
        help=f"用于保存过滤后的gzipped FASTQ文件以及ParaFly日志的输出目录。\n(默认: {DEFAULT_OUTPUT_DIR})" # Output directory... (Default: ...)
    )
    parser.add_argument(
        "--seq_c", default=DEFAULT_SEQ_C,
        help=f"用于过滤和移除的第一个目标序列。\n(默认: {DEFAULT_SEQ_C})" # First target sequence... (Default: ...)
    )
    parser.add_argument(
        "--seq_d", default=DEFAULT_SEQ_D,
        help=f"用于过滤和移除的第二个目标序列。\n(默认: {DEFAULT_SEQ_D})" # Second target sequence... (Default: ...)
    )
    parser.add_argument(
        "--num_cpus", type=int, 
        default=DEFAULT_NUM_CPUS, 
        help=f"供ParaFly使用的CPU核心数量。\n(默认: {DEFAULT_NUM_CPUS})" # Number of CPUs for ParaFly... (Default: ...)
    )
    parser.add_argument("--run_worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--worker_r1", help=argparse.SUPPRESS)
    parser.add_argument("--worker_r2", help=argparse.SUPPRESS)
    
    args = parser.parse_args()

    if args.run_worker:
        if not (args.worker_r1 and args.worker_r2 and args.output_dir and args.seq_c and args.seq_d):
            print("[Worker Error] 缺少必要的参数 (worker_r1, worker_r2, output_dir, seq_c, seq_d)。") # Missing necessary parameters
            exit(1)
        
        if not os.path.exists(args.output_dir):
            try:
                os.makedirs(args.output_dir, exist_ok=True)
            except OSError as e:
                print(f"[Worker Error] 无法创建输出目录 '{args.output_dir}': {e}") # Cannot create output directory
                exit(1)

        success = core_process_single_fastq_pair(
            args.worker_r1, args.worker_r2, 
            args.output_dir, args.seq_c, args.seq_d
        )
        exit(0 if success else 1)

    print("--- 启动FASTQ过滤与HiC-Pro执行脚本 (Orchestrator模式) ---") # Starting script (Orchestrator mode)
    print(f"重要提示: 请确保您已激活 'hicpro3' conda 环境，并且 ParaFly 和 run_hicpro 可用。") # Important: Ensure hicpro3 conda env is active and ParaFly & run_hicpro are available.

    if not args.input_dir:
        parser.error("-i/--input_dir is required in orchestrator mode.")
    # num_cpus has a default, so it will always be present for the orchestrator.

    if not os.path.isdir(args.input_dir):
        print(f"错误: 输入目录 '{args.input_dir}' 未找到或不是一个目录。") # Error: Input directory not found or not a directory.
        return
        
    absolute_output_dir = os.path.abspath(args.output_dir)
    if not os.path.exists(absolute_output_dir):
        try:
            os.makedirs(absolute_output_dir)
            print(f"已创建输出目录: {absolute_output_dir}") # Created output directory
        except OSError as e:
            print(f"错误: 无法创建输出目录 '{absolute_output_dir}': {e}") # Error: Cannot create output directory
            return
    elif not os.path.isdir(absolute_output_dir):
        print(f"错误: 输出路径 '{absolute_output_dir}' 已存在但不是一个目录。") # Error: Output path exists but is not a directory.
        return

    fastq_file_pairs = find_fastq_pairs(args.input_dir)
    if not fastq_file_pairs:
        print(f"在 '{args.input_dir}' 中未找到以 '_1.fq.gz' 和 '_2.fq.gz' 结尾的FASTQ文件对。") # No FASTQ pairs found.
        return

    print(f"\n找到 {len(fastq_file_pairs)} 对FASTQ文件准备通过ParaFly进行处理。") # Found pairs to process with ParaFly.
    print(f"将筛选并移除以下序列: SEQ_C='{args.seq_c}', SEQ_D='{args.seq_d}'") # Will filter and remove sequences...
    print(f"处理后的文件将保存到: {absolute_output_dir}") # Processed files will be saved to...
    print(f"ParaFly将使用 {args.num_cpus} 个CPU核心。") # ParaFly will use CPUs.

    command_file_path = os.path.join(absolute_output_dir, "parafly_commands.txt")
    current_script_path = os.path.abspath(__file__)

    with open(command_file_path, "w") as cmd_file:
        for r1_path, r2_path in fastq_file_pairs:
            command = (
                f"python {current_script_path} --run_worker "
                f"--worker_r1 \"{r1_path}\" --worker_r2 \"{r2_path}\" "
                f"-o \"{absolute_output_dir}\" "
                f"--seq_c \"{args.seq_c}\" --seq_d \"{args.seq_d}\""
            )
            cmd_file.write(command + "\n")
    
    print(f"\nParaFly 命令文件已生成: {command_file_path}") # ParaFly command file generated.

    parafly_successful = run_parafly(command_file_path, args.num_cpus, absolute_output_dir)

    if not parafly_successful:
        print("\nParaFly 处理步骤失败或报告了错误。将跳过 run_hicpro 命令。") # ParaFly step failed or reported errors. Skipping run_hicpro.
        return
    
    print("\n所有FASTQ文件通过ParaFly处理完成。") # All FASTQ files processed via ParaFly.

    print("\n准备执行 run_hicpro 命令...") # Preparing to execute run_hicpro.
    # IMPORTANT: Customize this command as needed for your HiC-Pro setup.
    # This is a placeholder. You will likely need to provide a config file and output dir for HiC-Pro.
    hicpro_command_to_run = ["run_hicpro"] 
    # Example of a more complete HiC-Pro command:
    # Ensure these paths are correct and accessible.
    # hicpro_config_file = "/path/to/your/hicpro_config.txt" # Make this an argument or a fixed known path
    # hicpro_results_dir = os.path.join(os.path.dirname(absolute_output_dir), "hicpro_analysis_results") # Example output for HiC-Pro
    # if not os.path.exists(hicpro_results_dir):
    #     os.makedirs(hicpro_results_dir, exist_ok=True)
    #
    # hicpro_command_to_run = [
    #     "run_hicpro", 
    #     "-c", hicpro_config_file,
    #     "-i", absolute_output_dir, # Input for HiC-Pro is the filtered FASTQ dir from this script
    #     "-o", hicpro_results_dir 
    # ]

    if run_final_command(hicpro_command_to_run, command_name="run_hicpro"):
        print("run_hicpro 命令已成功启动 (请检查其具体输出以确认完成情况)。") # run_hicpro command started successfully (check its output for completion).
    else:
        print("run_hicpro 命令执行失败或未找到。请检查您的环境和命令。") # run_hicpro command failed or not found. Check environment and command.
    
    print("\n--- 脚本执行完毕 ---") # Script execution finished.

if __name__ == "__main__":
    main()
