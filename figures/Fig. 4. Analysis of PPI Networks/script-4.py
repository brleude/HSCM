import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 读取三张图片
p1 = mpimg.imread('Fig 4a Count distributions in S.cerevisiae PPI network/Count distributions in S.cerevisiae PPI network.png')   
p2 = mpimg.imread('Fig 4b Count distributions in E.coli PPI network/Count distributions in E.coli PPI network.png')  
p3 = mpimg.imread('Fig 4c Comparisons of average count ratios/Comparisons of average count ratios.png')  

# 创建一个 3x1 的子图布局
fig = plt.figure(figsize=(11.69, 12))  # 设置大小
gs = plt.GridSpec(3, 1, height_ratios=[1, 1, 1])  # 设置三个子图的高度比例

# 绘制第一张图片
ax1 = fig.add_subplot(gs[0, 0])
ax1.imshow(p1)
ax1.axis('off')
# # 在图片下方添加编号 (a)
# ax1.text(0.5, -0.05, '(a)', transform=ax1.transAxes, ha='center', va='center', fontsize=18)

# 绘制第二张图片
ax2 = fig.add_subplot(gs[1, 0])
ax2.imshow(p2)
ax2.axis('off')
# # 在图片下方添加编号 (b)
# ax2.text(0.5, -0.05, '(b)', transform=ax2.transAxes, ha='center', va='center', fontsize=18)

# 绘制第三张图片
ax3 = fig.add_subplot(gs[2, 0])
ax3.imshow(p3)
ax3.axis('off')
# # 在图片下方添加编号 (c)
# ax3.text(0.5, -0.05, '(c)', transform=ax3.transAxes, ha='center', va='center', fontsize=18)

# 调整子图之间的间距
plt.tight_layout()

# 保存组合后的图片
plt.savefig('Analysis of PPI Networks.png', dpi=300)
plt.savefig('Analysis of PPI Networks.pdf')

# 显示图片
plt.show()