import scanpy as sc
import gc
import sys
import cellanova as cnova
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sea
import harmonypy as hm
import scanpy.external as sce
from scipy import stats

# 设置Scanpy的日志级别为0，即不显示详细信息
sc.settings.verbosity = 0
# 设置绘图参数，dpi为80
sc.settings.set_figure_params(dpi=80)
# 设置Pandas显示所有列
pd.set_option('display.max_columns', None)
# 设置随机种子
seed = 10
np.random.seed(seed)

# 假设从CSV文件中读取数据，文件路径已去除私人路径，需根据实际情况调整路径
lognorm_expr = pd.read_csv('your_data/lognorm_expr.csv', index_col = 0)
meta = pd.read_csv('your_data/meta.csv', index_col = 0)

# 创建AnnData对象
adata_prep = sc.AnnData(np.array(lognorm_expr).T, dtype=np.float64)
# 设置变量名
adata_prep.var_names = lognorm_expr.index.values
# 设置观测数据
adata_prep.obs = meta

# 数据标准化
sc.pp.scale(adata_prep)
# 将标准化后的数据存储在layers中
adata_prep.layers['scale'] = adata_prep.X

# 假设从CSV文件中读取主成分数据，文件路径已去除私人路径，需根据实际情况调整路径
pcs = pd.read_csv('your_data/pcs.csv', index_col = 0)

# 选择与观测数据匹配的主成分数据
Z_corr_combined = pcs.loc[adata_prep.obs_names]
print(np.all(Z_corr_combined.index == adata_prep.obs_names))

# 将主成分数据存储在obsm中
adata_prep.obsm['Cmat'] = Z_corr_combined.to_numpy()

# 检查数据的基本统计信息
print(adata_prep.X.max(), adata_prep.X.min(), adata_prep.X.mean(), adata_prep.X.std())

# 再次标准化数据（如果有必要）
sc.pp.scale(adata_prep)

# 将数据类型转换为np.float64，以避免溢出问题
adata_prep.X = adata_prep.X.astype(np.float64)

# 使用cellanova包计算ME（具体含义根据包的定义）
adata_prep = cnova.model.calc_ME(adata_prep, integrate_key='rep')

# 假设ip是一个Seurat对象
# 从Seurat对象的"data"插槽中获取数据
data <- GetAssayData(ip, slot = "data")
# 将数据写入CSV文件，路径已去除私人路径，需根据实际情况调整路径
write.csv(data, file = "your_data/sobj-WTmESC-data.csv", quote = F)

# 从Seurat对象的另一个插槽中获取数据（slot.data，具体含义根据实际情况）
data <- GetAssayData(ip, slot = "slot.data")
# 将数据写入CSV文件
write.csv((data), file = "your_data/sobj-WTmESC-scaledata.csv", quote = F)

# 获取Seurat对象的元数据
meta <- ip@meta.data
# 将元数据写入CSV文件
write.csv(meta, file = "your_data/sobj-WTmESC-metadata.csv", quote = F)

# 获取Seurat对象的PCA结果的前50个主成分
pc <- ip@reductions$pca[, 1:50]
# 将主成分数据写入CSV文件
write.csv(pc, file = "your_data/sobj-WTmESC-pc.csv", quote = F)
