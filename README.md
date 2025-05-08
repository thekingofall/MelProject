# MelProject
### version1
```
python runmel.py -i raw --seq_c ATGTCGGAGTTCTTAGGCGCTGGTTGCGTTGAGAACTCCGACT --seq_d AGTCGGAGTTCTTCAACGCAACCAGCGCCTAGAACTCCGACAT
```

### version2

```
python runmel2.py -i raw --seq_c ATGTCGGAGTTCTTAGGCGCTGGTTGCGTTGAGAACTCCGACT --seq_d AGTCGGAGTTCTTCAACGCAACCAGCGCCTAGAACTCCGACAT

```

--- 启动FASTQ过滤与HiC-Pro执行脚本 (Orchestrator模式) ---
重要提示: 请确保您已激活 'hicpro3' conda 环境，并且 ParaFly 和 run_hicpro 可用。
已创建输出目录: /data3/maolp/All_ZengXi_data5/20250502_fq/fastq/L7/Run1_fastq

找到 2 对FASTQ文件准备通过ParaFly进行处理。
将根据以下序列进行筛选 (序列本身不会被移除): SEQ_C='ATGTCGGAGTTCTTAGGCGCTGGTTGCGTTGAGAACTCCGACT', SEQ_D='AGTCGGAGTTCTTCAACGCAACCAGCGCCTAGAACTCCGACAT'
处理后的文件将保存到: /data3/maolp/All_ZengXi_data5/20250502_fq/fastq/L7/Run1_fastq
ParaFly将使用 10 个CPU核心。

ParaFly 命令文件已生成: /data3/maolp/All_ZengXi_data5/20250502_fq/fastq/L7/Run1_fastq/parafly_commands.txt

正在执行 ParaFly: ParaFly -c /data3/maolp/All_ZengXi_data5/20250502_fq/fastq/L7/Run1_fastq/parafly_commands.txt -CPU 10 -verbose -failed_cmds /data3/maolp/All_ZengXi_data5/20250502_fq/fastq/L7/Run1_fastq/parafly_failed_cmds.lo
