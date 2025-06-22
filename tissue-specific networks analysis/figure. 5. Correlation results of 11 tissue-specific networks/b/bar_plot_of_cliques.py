import networkx as nx
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)


def plot_clique_sizes_by_tissue(tissue_networks, output_path="clique_sizes_comparison.png"):
    """
    绘制多个组织网络中不同大小团的数量堆叠柱状图
    (图例右侧竖排/增加右侧空间/对数刻度/无边框设计)
    
    参数:
        tissue_networks: 字典，键为组织名称，值为networkx图对象
        output_path: 输出图像路径
    
    返回:
        各组织各大小团数量的统计数据
    """
    # 定义要统计的团大小
    clique_sizes = [2, 3, 4, 5]
    
    # 数据统计
    tissue_data = {}
    for tissue, G in tissue_networks.items():
        clique_counts = {size: 0 for size in clique_sizes}
        max_clique = max(len(c) for c in nx.find_cliques(G)) if G.nodes() else 0
        
        clique_counts[2] = G.number_of_edges()
        for size in range(3, max_clique + 1):
            if size not in clique_sizes:
                continue
            cliques = set()
            for clique in nx.find_cliques(G):
                if len(clique) >= size:
                    for sub_clique in combinations(clique, size):
                        if all(G.has_edge(u, v) for u, v in combinations(sub_clique, 2)):
                            cliques.add(tuple(sorted(sub_clique)))
            clique_counts[size] = len(cliques)
        tissue_data[tissue] = clique_counts
    
    # 准备绘图数据
    tissues = list(tissue_data.keys())
    log_data = np.log10(np.array([[tissue_data[t][size] + 1 for size in clique_sizes] for t in tissues]))
    
    # 创建画布（调整宽度适应右侧图例）
    plt.figure(figsize=(11, 6))  # 增加宽度
    
    # 颜色定义
    colors = {2: "#716767", 3: "#74A485", 4: "#2D6593", 5: "#4A4783"}
    
    # 绘制堆叠柱状图
    bottom = np.zeros(len(tissues))
    for i, size in enumerate(clique_sizes):
        plt.bar(np.arange(len(tissues)), log_data[:, i], 
                width=0.6, bottom=bottom,
                color=colors[size], label=f'{size}-node clique',
                zorder=2, alpha=0.7,
                edgecolor='black', linewidth=0.5)
        bottom += log_data[:, i]
    
    # 坐标轴设置
    ax = plt.gca()
    ax.set_xlabel('Tissue Type', fontsize=21)
    ax.set_ylabel(' log$_{10}$ (# of Cliques)', fontsize=21)
    ax.set_xticks(np.arange(len(tissues)))
    ax.set_xticklabels(tissues, rotation=45, ha='right', fontsize=21)
    ax.tick_params(axis='both', which='major', labelsize=21)
    
    # 移除边框
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    
    # 图例设置（右侧竖排）
    legend = ax.legend(
        loc='center left',          # 定位到左侧中心
        bbox_to_anchor=(1.02, 0.5), # 移到画布右侧
        frameon=False,
        prop={'size': 21},
        # title_fontsize='12',
        handlelength=1,
        markerscale=0.8
    )
    
    # 调整布局（为右侧图例留空间）
    plt.tight_layout(rect=[0, 0, 0.95, 1])  # 右侧保留15%空间
    # plt.show()
    # 保存和返回
    plt.savefig(output_path)
    plt.close()
    return tissue_data


# 组织名称列表
tissues = [
    "brain", "blood", "breast", "colon", 
    "kidney", "liver", "lung", "ovary", 
    "pancreas", "stomach", "throat"
]

# 加载所有组织网络
tissue_networks = {}
for tissue in tissues:
    network_path = os.path.join('data', f'{tissue}_8_network.graphml')
    try:
        G = nx.read_graphml(network_path)
        tissue_networks[tissue] = G
        print(f"已加载 {tissue} 网络: {G.number_of_nodes()} 节点, {G.number_of_edges()} 边")
    except Exception as e:
        print(f"无法加载 {tissue} 网络: {e}")

# 绘制堆叠柱状图
if tissue_networks:
    tissue_data = plot_clique_sizes_by_tissue(tissue_networks, "number of cliques in 11 tissue PPI networks.pdf")
