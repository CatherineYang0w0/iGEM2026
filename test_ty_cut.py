import re

import csv

def write_peptides_csv(peptides, filename):
    """
    将肽段写入CSV文件，包含length和sequence两列
    """
    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["length", "sequence"])
        for pep in peptides:
            writer.writerow([len(pep), pep])

def write_peptides_fasta(peptides, filename):
    """
    将肽段写入FASTA文件，描述行包含length信息
    """
    with open(filename, "w") as f:
        for i, pep in enumerate(peptides, 1):
            f.write(f">peptide_{i} length={len(pep)}\n{pep}\n")

def virtual_trypsin_digest(protein_seq, min_len=8, max_len=12):
    # 胰蛋白酶切点：K 或 R 的 C 端（即 K 或 R 之后切断）
    # 正则表达式：(?<=[KR])(?=[^P]) 表示匹配 K 或 R 后面不跟 P 的位置
    cleavage_sites = [m.start() + 1 for m in re.finditer(r'(?<=[KR])(?=[^P])', protein_seq)]

    # 按切点分割蛋白序列，收集所有满足长度要求的肽段
    peptides = []
    start = 0
    for site in cleavage_sites:
        peptide = protein_seq[start:site]
        if min_len <= len(peptide) <= max_len:
            peptides.append(peptide)
        start = site
    # 处理最后一段
    if start < len(protein_seq):
        peptide = protein_seq[start:]
        if min_len <= len(peptide) <= max_len:
            peptides.append(peptide)

    return peptides

# 新增：读取FASTA文件，拼接所有序列行
def read_fasta_sequence(filepath):
    """
    读取FASTA文件，拼接所有非注释行，返回蛋白质序列字符串
    """
    sequence = ''
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            sequence += line.strip()
    return sequence

if __name__ == "__main__":
    # 指定FASTA文件名
    fasta_file = "sequences.FASTA"
    # 读取蛋白质序列
    protein_seq = read_fasta_sequence(fasta_file)
    # 进行trypsin酶切，获得所有肽段
    peptides = virtual_trypsin_digest(protein_seq)

    # 选择输出格式：csv 或 fasta
    output_format = input("请选择输出格式（csv 或 fasta）：").strip().lower()
    if output_format == "csv":
        write_peptides_csv(peptides, "digestedpeptides.csv")
        print("已输出到 digestedpeptides.csv")
    elif output_format == "fasta":
        write_peptides_fasta(peptides, "digestedpeptides.fasta")
        print("已输出到 digestedpeptides.fasta")
    else:
        print("不支持的格式，请输入 csv 或 fasta")