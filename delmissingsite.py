#!/usr/bin/python3

import numpy as np
from pandas.core.frame import DataFrame
from multiprocessing import Pool
from Bio import SeqIO
import os, time, random
import argparse
import sys

alignment_len = 50000
normal_site_list = ["A","a","T","t","C","c","G","g"]
'''
解析参数
'''
parser = argparse.ArgumentParser(description="Options for cp_alignment.py",
                                 add_help=True)
additional = parser.add_argument_group("additional arguments")    
additional.add_argument('-p', '--proportion', action="store", type=float, 
                        default=0.2, metavar='\b', help='''The proportion of 
                        missing data allowed in each site, default = 0.2''')
additional.add_argument('-n', '--num_cpu', action="store", type=int, default=40,
                        metavar='\b', help='''The maximum number of CPUs that 
                        this script can be used (Two usable CPUs means that the
                        script will analyze two matrices at the same time), 
                        default = 12''')
args = parser.parse_args()
proportion = args.proportion
th         = args.num_cpu

if proportion == 0: #该参数会作为被除数使用，故不能为零
    proportion = 0.0001



def get_file_list(ext):
    '''
    函数get_file_list: 
        获取当前文件夹中符合目标扩展名的文件
    输入: 
        当前工作文件夹中的所有文件名名称(脚本自动获取)
    输出: 
        file_name：扩展名符合的文件的名称列表
    '''          
    file_name = []
    ext_len = len(ext)       
    for each in os.listdir(os.getcwd()):      
        if ext == each[-ext_len:]:
            file_name.append(each)
    return file_name

def get_th_list(file_name):
    '''
    函数 get_th_list：
        将所有文件平均的分给各个进程
    输入 
        file_name: 需要进行分析的所有文件的名称列表
    输出 
        th_list：列表，列表中的各个元素为各个进程所应该处理的文件
    '''
    th_list = []
    for each_num in range(th):
        th_list.append([])
    print("本程序共使用 " + str(th) + " 个进程")
    n = 0
    for each_file in file_name:       
        th_list[n].append(each_file)
        if n == th - 1:
            n = 0
        else:
            n = n + 1
    return th_list

def if80to1(fasta_name):
    '''
    函数if80to1: 
        如果fasta序列是80列换行的话则修改成不换行
        同时，删除所有序列都是“-”“?”或“N”的样本
    输入：
        待处理的fasta文件名称
    输出：
        如果原始fasta中的序列是80行换行的则将所有序列集中至一行，并写出一个新的
        后缀为“.fa”格式的文件
        如果序列的长度超过“alignment_len”中设定的长度，则返回一个真值，准备开始分割这个文件
    '''
    tmp_len_list = []
    valid_records = []  # 存储有效样本的记录
    for each_record in SeqIO.parse(fasta_name, "fasta"):
        seq_str = str(each_record.seq)
        # 检查序列是否完全由“-”、“?”或“N”组成
        if not all(c in "-?N" for c in seq_str):
            valid_records.append(each_record)
            tmp_len_list.append(len(seq_str))
    
    # 如果没有有效样本，删除原文件并返回 False
    if not valid_records:
        print(f"Warning: {fasta_name} has no valid samples after filtering.")
        os.remove(fasta_name)
        return False
    
    # 计算最大序列长度
    tmp_len_list.sort()
    max_len = tmp_len_list[-1]
    
    # 写入新的 FASTA 文件（后缀为 .fa）
    new_fasta_name = fasta_name[:-6] + ".fa"  # 假设输入文件后缀是 .fasta
    with open(new_fasta_name, "w") as write_file:
        for each_record in valid_records:
            seq_str = str(each_record.seq)
            # 如果序列长度小于最大长度，用“-”补齐
            if len(seq_str) < max_len:
                gap = "-" * (max_len - len(seq_str))
                each_record.seq = each_record.seq + gap
            write_file.write(">" + str(each_record.id) + "\n")
            write_file.write(str(each_record.seq) + "\n")
    
    
    # 检查是否需要分割文件
    if max_len >= alignment_len:
        return True
    else:
        return False                         
             
def split_fasta(fasta_name):
    '''
    函数split_fasta: 
        将排序文件每隔alignment_len bp分割一次
    输入：
        待处理的fasta文件名称（后缀为 .fa）
    输出：
        后缀为“.fa”，并且文件名中带有“.split.”字样的文件
    '''
    with open(fasta_name, "r") as read_file:
        sequences = read_file.readlines()
        length = len(sequences[1]) - 1
        left = 0
        right = alignment_len
        while True:
            new_name = (fasta_name[:-3] + ".split." + str(left) + ".fa")
            if right <= length:  
                with open(new_name, "w") as write_file:
                    for each_line in sequences:
                        if each_line[0] == ">":
                            write_file.write(each_line)
                        else:
                            write_file.write(each_line[left:right] + "\n")
                left = left + alignment_len
                right = right + alignment_len
            else:
                with open(new_name, "w") as write_file:
                    for each_line in sequences:
                        if each_line[0] == ">":
                            write_file.write(each_line)
                        else:
                            write_file.write(each_line[left:])
                break

def get_seq_name_list(fasta_name):
    '''
    函数get_seq_name_list: 
        获得alignment文件中各个序列的序列名称
    输入：
        待处理的fasta文件名称
    输出：
        seq_name_list包含有所有序列名称的一个列表文件
    '''   
    seq_name_list = []
    with open(fasta_name) as read_file:
        for each_line in read_file:
            if each_line[0] == ">":
                seq_name_list.append(each_line)
    return seq_name_list

def calculate(fasta_name, proportion):
    '''
    函数calculate: 
        计算alignment文件各位点上gap所占的比例，如果大于设定的阈值，则删除该位点
    输入：
        待处理的fasta文件名称
    输出：
        result_list: 包含所有被删除过gap的序列的列表
    '''      
    seq_array_temp_list = []
    with open(fasta_name) as read_file:
        for each_line in read_file:
            if each_line[0] != ">":
                if each_line[-1] == "\n":
                    seq_array_temp_list.append(list(each_line[:-1]))
                else:
                    seq_array_temp_list.append(list(each_line))
    seq_array = DataFrame(seq_array_temp_list)
    row = seq_array.shape[0]
    column = seq_array.shape[1]
    temp_list = []
    for each_num in range(0,column):
        gap_num = 0
        for site in list(seq_array[each_num]):
            if site not in normal_site_list:
                gap_num = gap_num + 1
        if gap_num/row >= proportion:   
            pass
        else:
            temp_list.append(list(seq_array[each_num]))
    temp_seq_array = DataFrame(temp_list)
    result_list = []
    for each_num in range(0,row):
        try:
            result_list.append("".join(list(temp_seq_array[each_num].values)))
        except:
            result_list = []
    return(result_list)

def white2file(fasta_name, seq_name_list, seq_list):
    '''
    函数white2file: 将计算结果写入输出文件中
    输入：
        fasta_name: 待处理的fasta文件名称； 
        seq_name_list：包含有所有序列名称的列表文件；
        seq_list：包含所有被删除过gap的序列的列表
    输出：
        写出一个后缀为.fas的结果文件
    '''   
    with open(fasta_name + "s", "a") as write_file:
        seq_num = len(seq_name_list)
        for num in range(0,seq_num):
            judge1 = seq_list[num].count("-") == len(seq_list[num])
            judge2 = seq_list[num].count("?") == len(seq_list[num])
            if judge1 or judge2:
                pass
            else:
                write_file.write(seq_name_list[num])
                write_file.write(seq_list[num] + "\n")

def main_get_homo(file_name):
    '''
    函数main_get_homo: 
        控制前面函数运行的函数,删掉alignment中missing site
    输入：
        “.fa”格式文件，该文件由preprocessing函数生成
    输出：
        “.fas”结果文件
    '''  
    for fasta_name in file_name:
        seq_name_list = get_seq_name_list(fasta_name)
        seq_list = calculate(fasta_name, proportion)
        if len(seq_list) != 0:  #删完gap之后删掉只剩下"-"的物种
            white2file(fasta_name, seq_name_list, seq_list)
            os.remove(fasta_name)
        else:
            os.remove(fasta_name)

def preprocessing(file_name):
    '''
    函数preprocessing: 
        1：如果是80行换行的fasta，则转化为不换行
        2：如果矩阵的长度超过"alignment_len"bp，则将矩阵分割成最长为"alignment_len"bp的小矩阵
        3：删除所有序列都是“-”“?”或“N”的样本
    输入：
        自动输入文件夹中“.fasta”格式文件
    输出：
        “.fa”格式的中间文件，最终生成“.fas”格式的结果文件
    '''  
    for fasta_name in file_name:
        if if80to1(fasta_name):
            new_fasta_name = fasta_name[:-6] + ".fa"
            split_fasta(new_fasta_name)
            os.remove(new_fasta_name)


def fasta2dict(fasta_file):
    #将fasta文件转化成一个python列表
    res_dict = {}
    with open(fasta_file, "r") as read_file:
        for each_line in read_file:
            if each_line[0] == ">":
                seq_name = each_line.split(" ")[0][1:].replace("\n", "").replace("/","_").replace("\\","_")
                res_dict[seq_name] = ""
            else:
                res_dict[seq_name] = res_dict[seq_name] + each_line.replace("\n", "")
    return res_dict, len(res_dict[seq_name])

def concat(gene_name_list):
    '''
    函数concat: 
        脚本最初会将过长的矩阵分割成较小的部分进行并行处理，这里需要将这些被分割的
    文件重新合并成一个文件
    输入：
        基因名称的列表，不含有后缀
    输出：
        “.fas”格式的结果文件
    '''  
    for each_gene in gene_name_list:
        fasta_file_list = []
        for each_file_name in os.listdir(os.getcwd()):
            if each_gene + ".split." in each_file_name:
                fasta_file_list.append(each_file_name)
        
        concat_list = []
        name_list = []
        for each_file in fasta_file_list:
            fasta_dict = fasta2dict(each_file)
            concat_list.append(fasta_dict)
            for each_name in fasta_dict[0]:
                if each_name not in name_list:
                    name_list.append(each_name)
        for each_name in name_list:
            with open(each_gene + each_name + ".temp", "a") as write_file:
                write_file.write(">" + each_name + "\n")
                for each_dict in concat_list:
                    if each_name in each_dict[0]:
                        write_file.write(each_dict[0][each_name])
                    else:
                        write_file.write("?"*each_dict[1])
                write_file.write("\n")
        os.system("cat " + each_gene + "*.temp > " + each_gene + ".fas")
        os.system("rm " + each_gene + "*.temp")
        os.system("rm " + each_gene + ".split*")
        


if __name__=="__main__":
    '''
    预处理文件：
    1：如果是80行换行的fasta，则转化为不换行的
    2：如果矩阵的长度超过2000bp，则将矩阵分割成最长为2000bp的小矩阵
    '''
    file_name = get_file_list(".fasta")
    th_list = get_th_list(file_name)
    p = Pool(th)
    for i in range(th):
        p.apply_async(preprocessing,(th_list[i],))
    p.close()  # 关闭进程池，关闭后po不再接收新的请求
    p.join()  # 等待po中所有子进程执行完成，再执行下面的代码

    '''
    并行删除alignment格式文件中的missing site
    '''
    file_name = get_file_list(".fa")
    th_list = get_th_list(file_name)
    p = Pool(th)
    for i in range(th):
        p.apply_async(main_get_homo,(th_list[i],))
    p.close()  # 关闭进程池，关闭后po不再接收新的请求
    p.join()  # 等待po中所有子进程执行完成，再执行下面的代码

    '''
    合并最开始被分成小片段的矩阵
    '''
    gene_name_list = []
    for each_file_name in os.listdir(os.getcwd()): 
        if ".split." in each_file_name:
            gene_name = each_file_name.split(".split")[0]
            if gene_name not in gene_name_list:
                gene_name_list.append(gene_name)
    th_list = get_th_list(gene_name_list)
    p = Pool(th)
    for i in range(th):
        p.apply_async(concat,(th_list[i],))
    p.close()  # 关闭进程池，关闭后po不再接收新的请求
    p.join()  # 等待po中所有子进程执行完成，再执行下面的代码

