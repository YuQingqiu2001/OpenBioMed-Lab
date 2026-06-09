# 生信绘图 — 全网精选方法与工具

> 搜索时间: 2026-06-06 | 来源: 知乎/CSDN/简书/微信公众号 | 分类整理

## 📊 核心绘图类型与最佳实践

### 1. 火山图 (Volcano Plot)
顶刊级优化：对称渐变配色(`scale_color_gradient2`)、显著性阈值标注、Top基因标注(`ggrepel`)、透明度映射、边缘密度图(`ggExtra`)。关键包: `EnhancedVolcano`, `ggplot2`, `ggrepel`。

### 2. 热图 (Heatmap)
美观要点：高级配色(`viridis`/`ggsci`)、聚类树美化(`ggtree`)、多层注释(`ComplexHeatmap::HeatmapAnnotation`)、环状热图(`circlize`)、交互式(`plotly`/`heatmaply`)。关键包: `ComplexHeatmap`(顶刊首选), `pheatmap`。

### 3. 小提琴图 (Violin/Box Plot)
分半小提琴(`geom_split_violin`/`see`)、叠加散点(`geom_jitter`+alpha)、分组配色(`scale_fill_manual`)、排序优化(`fct_reorder`)。关键包: `ggplot2`, `ggpubr`, `see`, `ggstatsplot`。

### 4. 气泡图 (Dot Plot)
Seurat `DotPlot()` / Scanpy `sc.pl.dotplot()`。大小=表达比例，颜色=平均表达。`scale_color_viridis()`。

### 5. UMAP/t-SNE
`theme_void()`去轴、连续/离散配色、密度等高线(`geom_density_2d`)、箭头标注。

## 🎨 配色方案
| 包名 | 来源 | 用途 |
|------|------|------|
| `ggsci` | Nature/Science/Lancet/NEJM/JAMA | 一键顶刊配色 |
| `viridis` | matplotlib同款 | 色盲友好+连续渐变 |
| `RColorBrewer` | ColorBrewer | 经典调色板 |
| `wesanderson` | Wes Anderson电影 | 复古文艺 |
| `MetBrewer` | 大都会博物馆 | 艺术配色 |
| `paletteer` | 综合 | 统一接口50+调色板 |

**配色黄金法则**: 连续→viridis; 分类→色盲友好(≤8类); 热图→蓝-白-红diverging; 禁止→红绿搭配+jet彩虹色。

## 🧰 生信专用工具
**Python**: `scanpy`, `squidpy`, `matplotlib`+`seaborn`, `plotly`, `pyComplexHeatmap`
**R**: `ggplot2`, `ComplexHeatmap`, `ggtree`, `clusterProfiler`, `enrichplot`, `Seurat`

## 🖼️ 空间转录组专属
`squidpy`, `Giotto`, `Seurat SpatialFeaturePlot`, `SPATA2`, `SpatialDE`, `stLearn`

## 📱 微信公众号精选
- 生信常用分析图形绘制系列(热图/火山图/聚类树/桑基图)
- 从《甄嬛传》学科研绘图配色
- ggplot2扩展包精选集-生信实战指南
- SPATA: 基因集驱动的空间转录组分析框架

## 🎯 发表级Checklist
- [ ] 配色来自顶刊调色板(ggsci/viridis)
- [ ] 字体清晰≥8pt, 无chartjunk
- [ ] 图例不重叠, 轴标签含义明确
- [ ] 统计检验标注(p值/显著性)
- [ ] 导出矢量格式(PDF/SVG > 300dpi TIFF)
- [ ] 色盲友好

## 📚 学习路径
1. 入门: 微信公众号「生信师兄」「生信技能树」
2. 进阶: CSDN搜索"ggplot2顶刊风格"
3. 精通: R Graph Gallery + ggplot2官方文档
4. 空间组学: squidpy教程 + SPATA2 vignettes
5. 配色灵感: ggsci文档+MetBrewer+Nature/Science截图取色
