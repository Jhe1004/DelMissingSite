import numpy as np
from pandas.core.frame import DataFrame
from multiprocessing import Pool
import os, time, random
from Bio import SeqIO
import argparse
import sys

def get_file_list():
    '''
    函数get_file_list: 获取当前文件夹中符合目标扩展名的文件
    输入: 无，将本脚本放置在目标文件夹中即可
    输出: file_name：所有文件的名称列表
    '''          
    file_name = []       
    for each in os.listdir(os.getcwd()):        
        if ".fasta" in each:
            file_name.append(each)
    return file_name


def if80to1(fasta_name):
    '''
    函数if80to1: 如果序列是80列换行的话则修改成不换行
    输入：待处理的fasta文件名称
    输出：如果原始fasta中的序列是80行换行的则将所有序列集中至一行。并写出一个新的后缀为fa格式的文件
    '''
    tmp_len_list = []
    for each_record in SeqIO.parse(fasta_name, "fasta"):
        tmp_len_list.append(len(str(each_record.seq)))
    tmp_len_list.sort()
    with open(fasta_name[:-3], "a") as write_file:
        for each_record in SeqIO.parse(fasta_name, "fasta"):
            if len(str(each_record.seq)) != tmp_len_list[-1]:
                each_record.seq = each_record.seq + "-"*(tmp_len_list[-1] - len(str(each_record.seq)))
                write_file.write(">" + str(each_record.id) + "\n")
                write_file.write(str(each_record.seq) + "\n")    
            else:
                write_file.write(">" + str(each_record.id) + "\n")
                write_file.write(str(each_record.seq) + "\n")                            
             

def get_seq_name_list(fasta_name):
    '''
    函数get_seq_name_list: 获得alignment文件中各个序列的序列名称
    输入：待处理的fasta文件名称
    输出：seq_name_list包含有所有序列名称的一个列表文件
    '''   
    seq_name_list = []
    with open(fasta_name[:-3]) as read_file:
        for each_line in read_file:
            if each_line[0] == ">":
                seq_name_list.append(each_line)
    return seq_name_list


def calculate(fasta_name, proportion):
    '''
    函数calculate: 计算alignment文件各位点上gap所占的比例，如果大于设定好的阈值，则删除该位点
    输入：待处理的fasta文件名称
    输出：result_list: 包含所有被删除过gap的序列的列表
    '''      
    seq_array_temp_list = []
    with open(fasta_name[:-3]) as read_file:
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
        gap_num = list(seq_array[each_num]).count("-") + list(seq_array[each_num]).count("?")
        if gap_num/row >= proportion:            #gap占比阈值，如果某个位点中gap所占的比例大于该阈值，则删除该位点
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
    输入：fasta_name: 待处理的fasta文件名称； seq_name_list：包含有所有序列名称的列表文件；seq_list：包含所有被删除过gap的序列的列表
    输出：写出一个后缀为.fas的结果文件
    '''     
    with open(fasta_name[:-2], "a") as write_file:
        seq_num = len(seq_name_list)
        for num in range(0,seq_num):
            if seq_list[num].count("-") == len(seq_list[num]) or seq_list[num].count("?") == len(seq_list[num]):
                pass
            else:
                write_file.write(seq_name_list[num])
                write_file.write(seq_list[num] + "\n")


def main_get_homo(file_name):
    '''
    函数main_get_homo: 控制前面函数运行的函数
    输入：待处理的文件名称
    输出：写出的.fas结果文件
    '''  


    for fasta_name in file_name:
        if80to1(fasta_name)
        seq_name_list = get_seq_name_list(fasta_name)
        seq_list = calculate(fasta_name, proportion)
        if len(seq_list) != 0:  #删完gap之后有些物种的序列可能会只剩下"-"，这些物种需要删掉
            white2file(fasta_name, seq_name_list, seq_list)
            os.remove(fasta_name[:-3])
        else:
            os.remove(fasta_name[:-3])


'''
解析参数
'''
parser = argparse.ArgumentParser(description="Options for cp_alignment.py",
                                            add_help=True)
additional = parser.add_argument_group("additional arguments")    
additional.add_argument('-p', '--proportion', action="store", type=float, default=0.1,
                        metavar='\b', help="The proportion of missing data allowed in each site, default = 0.2, If you do not\
                            allow any missing data in the matrix, set this parameter to 0.00001, but do not set it to 0")
additional.add_argument('-n', '--num_cpu', action="store", type=int, default=12,
                        metavar='\b', help="The maximum number of CPUs that this script can be used (Two usable CPUs means that the\
                             script will analyze two matrices at the same time), default = 6")
args = parser.parse_args()
proportion = args.proportion
num_cpu    = args.num_cpu



th = num_cpu

if proportion == 0:
    proportion = 0.0001

def get_th_list(file_name):
    '''
    函数 get_th_list: 将所有文件平均的分给各个进程
    输入 file_name: 所有文件的名称列表
    输出 th_list： 一个列表，列表中的各个元素为各个进程所应该处理的文件
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

file_name = get_file_list()

th_list = get_th_list(file_name)
print(th_list)
p = Pool(th)
for i in range(th):
    p.apply_async(main_get_homo,(th_list[i],))

print("----start----")
p.close()  # 关闭进程池，关闭后po不再接收新的请求
p.join()  # 等待po中所有子进程执行完成，再执行下面的代码,可以设置超时时间join(timeout=)
print("-----end-----") 

