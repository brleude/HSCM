import pandas as pd
import networkx as nx
import sys
from pathlib import Path
python_code_path = str(Path(__file__).parent.parent.parent)
if python_code_path not in sys.path:
    sys.path.append(python_code_path)
from save_load_scores import save_scores
from centrality_improvement import improved_centrality
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 加载保存的edgelist文件
network_path = os.path.join('data', 'throat_8_network.edgelist')  # "brain", "blood", "breast", "colon", "kidney", "liver", "lung", "ovary", "pancreas", "stomach", "throat"
loaded_network = nx.read_edgelist(network_path)
print(f"Loaded network has {loaded_network.number_of_nodes()} nodes and {loaded_network.number_of_edges()} edges.")

# 不优化参数
_, csr_scores = improved_centrality(loaded_network, 'dc', max_clique=5, params=None)

csr_throat_8_dc = sorted(csr_scores.items(), 
                             key=lambda item: item[1], 
                             reverse=True)
# print(cdr_brain_5_dc)
print(csr_throat_8_dc)

# 保存CSR得分
save_scores(csr_throat_8_dc, 'csr_throat_8_dc.txt')


def save_csr_with_rank(csr_data, output_file):
    """
    将排序后的CSR数据保存为带排名的表格
    csr_data: 排序后的CSR数据列表，格式为[(gene, score), ...]
    output_file: 输出文件路径
    """
    # 提取基因和得分
    genes = [item[0] for item in csr_data]
    scores = [item[1] for item in csr_data]
    
    # 计算排名（相同得分共享排名）
    ranks = []
    current_rank = 1
    prev_score = None
    
    for i, score in enumerate(scores):
        if score != prev_score and i > 0:
            current_rank = i + 1
        ranks.append(current_rank)
        prev_score = score
    
    # 创建DataFrame并保存
    df = pd.DataFrame({
        'gene': genes,
        'score': scores,
        'rank': ranks
    })
    df.to_csv(output_file, index=False)
    print(f"表格已保存至 {output_file}")
    return df

# 保存为表格（示例路径，可根据需求修改）
output_path = "csr_throat_8_dc_with_rank1.csv"

save_csr_with_rank(csr_throat_8_dc, output_path)