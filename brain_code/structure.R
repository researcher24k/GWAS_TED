# 加载必要的库
library(ggplot2)

# 读取数据（请根据实际路径调整文件位置）
data <- read.csv("FN1.csv")

# 将数据从长格式转换为宽格式，以便分别绘制两条曲线
data_long <- reshape(data, 
                     varying = c("FN1", "comlex_FN1"), 
                     v.names = "value", 
                     timevar = "type", 
                     times = c("FN1", "complex_FN1"), 
                     direction = "long")
# 计算存在的残基中每隔50的数值作为断点
breaks <- unique(data_long$residue)
breaks <- breaks[breaks %% 50 == 0]

# 创建绘图
# 将残基转换为因子并按顺序排列
data_long$residue_factor <- factor(data_long$residue, levels = unique(data_long$residue))

# 每隔50个残基生成标签（基于顺序而非数值）
breaks <- seq(1, nlevels(data_long$residue_factor), by = 50)
labels <- levels(data_long$residue_factor)[breaks]

# 创建绘图（等距排列）
ggplot(data_long, aes(x = residue_factor, y = value, color = type, group = type)) +
  geom_line(size = 1) +
  scale_x_discrete(
    breaks = labels,  # 显示每隔50个位置的标签
    labels = labels   # 使用实际残基数值作为标签
  ) +
  labs(title = "CABSFlex 2.0", x = "Residue", y = "RMSF [Å]") +
  scale_color_manual(values = c("navy", "darkred")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
