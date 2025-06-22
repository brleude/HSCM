import pandas as pd
import numpy as np
from itertools import combinations
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)


def process_network_data(file_path):
    """
    读取网络节点数据文件，打印节点数量，并返回DataFrame
    
    参数:
        file_path (str): 数据文件路径（需包含双层方括号列表数据）
    
    返回:
        pandas.DataFrame: 包含节点和中心性得分的DataFrame
    """
    # 读取文件内容
    with open(file_path, 'r') as f:
        content = f.read()
    
    # 提取内层列表数据
    start_idx = content.find('[["')
    end_idx = content.rfind(']]')
    
    if start_idx == -1 or end_idx == -1:
        raise ValueError("文件中未找到有效的双层列表数据格式")
    
    raw_data = eval(content[start_idx:end_idx+2])
    
    # 转换为DataFrame并去重
    df = pd.DataFrame(raw_data, columns=['node', 'centrality_score'])
    df = df.drop_duplicates(subset=['node'])  # 按节点去重
    
    # 打印节点数量
    print(f"{file_path}中的总节点数为：{len(df)}")
    
    return df

def rank_correlation_matrix(df_dict):
    """
    计算各组织节点得分排名的相关性矩阵（基于共同节点的Spearman秩相关）
    
    参数:
        df_dict (dict): 组织名称到DataFrame的映射，每个DataFrame含'node'和'centrality_score'列
    
    返回:
        pandas.DataFrame: 相关性矩阵
    """
    #  对每个组织的节点按得分排名
    rank_dict = {}
    for tissue, df in df_dict.items():
        # 去重并排序（确保每个节点唯一）
        df_unique = df.drop_duplicates(subset=['node']).sort_values(
            'centrality_score', ascending=False, ignore_index=True
        )
        # 生成排名（从1开始，得分相同则取平均排名）
        df_unique['rank'] = df_unique['centrality_score'].rank(
            method='average', ascending=False
        )
        rank_dict[tissue] = df_unique.set_index('node')['rank'].to_dict()
    
    # 生成所有组织对组合
    tissues = list(rank_dict.keys())
    pairs = combinations(tissues, 2)
    n = len(tissues)
    corr_matrix = pd.DataFrame(np.zeros((n, n)), index=tissues, columns=tissues)
    
    # 计算两两组织的Spearman秩相关系数（基于共同节点）
    for i, tissue1 in enumerate(tissues):
        for j, tissue2 in enumerate(tissues):
            if i == j:
                corr_matrix.loc[tissue1, tissue2] = 1.0  # 对角线为1
                continue
            
            # 获取共同节点
            common_nodes = set(rank_dict[tissue1].keys()) & set(rank_dict[tissue2].keys())
            if len(common_nodes) < 2:  # 至少需要2个节点计算相关系数
                corr_matrix.loc[tissue1, tissue2] = np.nan
                continue
            
            # 构建共同节点的排名对
            rank_pairs = pd.DataFrame({
                'rank1': [rank_dict[tissue1][node] for node in common_nodes],
                'rank2': [rank_dict[tissue2][node] for node in common_nodes]
            })
            
            # 计算Spearman相关系数（rho值）
            corr = rank_pairs['rank1'].corr(rank_pairs['rank2'], method='spearman')
            corr_matrix.loc[tissue1, tissue2] = corr
            corr_matrix.loc[tissue2, tissue1] = corr  # 矩阵对称
    
    return corr_matrix


csr_brain_8_dc = process_network_data('data/csr_brain_8_dc.txt')
csr_blood_8_dc = process_network_data('data/csr_blood_8_dc.txt')
csr_breast_8_dc = process_network_data('data/csr_breast_8_dc.txt')
csr_colon_8_dc = process_network_data('data/csr_colon_8_dc.txt')
csr_kidney_8_dc = process_network_data('data/csr_kidney_8_dc.txt')
csr_liver_8_dc = process_network_data('data/csr_liver_8_dc.txt')
csr_lung_8_dc = process_network_data('data/csr_lung_8_dc.txt')
csr_ovary_8_dc = process_network_data('data/csr_ovary_8_dc.txt')
csr_pancreas_8_dc = process_network_data('data/csr_pancreas_8_dc.txt')
csr_stomach_8_dc = process_network_data('data/csr_stomach_8_dc.txt')
csr_throat_8_dc = process_network_data('data/csr_throat_8_dc.txt')

df_dict = {
    'brain': csr_brain_8_dc,
    'blood': csr_blood_8_dc,
    'breast': csr_breast_8_dc,
    'colon': csr_colon_8_dc,
    'kidney': csr_kidney_8_dc,
    'liver': csr_liver_8_dc,
    'lung': csr_lung_8_dc,
    'ovary': csr_ovary_8_dc,
    'pancreas': csr_pancreas_8_dc,
    'stomach': csr_stomach_8_dc,
    'throat': csr_throat_8_dc
}

# 计算排名相关性矩阵
corr_matrix = rank_correlation_matrix(df_dict)

# 保存为CSV文件
corr_matrix.to_csv('correlation_matrix.csv', index=True)