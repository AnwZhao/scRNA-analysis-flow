import numpy as np
import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as pl
from matplotlib import rcParams
import matplotlib
matplotlib.use('TkAgg')

path=''
adata=''
results_file = 'output.h5ad'

def read_1(path):
    os.chdir(path)
    upDir = os.path.pardir
    os.chdir(upDir)
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3) #设置日志等级
    sc.logging.print_header()  # 打印可能影响数值结果的版本。 Matplotlib和Seaborn被排除在外
    sc.settings.set_figure_params(dpi=80, frameon=True, fontsize=10, facecolor='white')  # 设置输出图形的分辨率/大小、样式和格式


    adata = sc.read_10x_mtx(       #读取 10x-Genomics 格式的 mtx 目录
        path,   #'data',  # the directory with the `.mtx` file
        var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
        cache=True)                              # write a cache file for faster subsequent reading  # If False, read from source, if True, read from fast ‘h5ad’ cache.
    return adata

def highest_2(adata,n_top=20):
    adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
    sc.pl.highest_expr_genes(adata, n_top=n_top, log=False, save='.png')
    # 计算每个基因在细胞内分配给该基因的计数分数。在所有细胞中具有最高平均分数的基因被绘制成箱线图。
    '''看一下top20基因的表达情况'''
    return adata

def filter_3(adata,min_genes=200,min_cells=3):
    '''过滤基因和细胞'''
    sc.pp.filter_cells(adata, min_genes=min_genes)  # 根据表达的基因计数和数量过滤细胞异常值,细胞通过过滤所需的最小表达基因数设为200

    sc.pp.filter_genes(adata, min_cells=min_cells)  # 根据细胞数量或计数过滤基因，基因通过过滤所需的最小表达细胞数为3个

    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    # 计算质量控制指标,percent_top覆盖多少比例的顶级基因，默认为前50个表达最多的基因的累积比例

    # 找出线粒体基因
    mito_genes = adata.var_names.str.startswith('MT-')
    adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    sc.pl.scatter(adata, x='total_counts', y='percent_mito')
    '''线粒体的含量图'''

    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'percent_mito'], jitter=0.4, multi_panel=True,
                 save='QC.png')
    '''用小提琴图将质量信息可视化'''

    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='.png')
    '''细胞和基因的关系，一般是线性关系，斜率越大越好，说明我们可以用较少的细胞测到较多的基因：'''
    return adata

def normlog_4(adata):
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 获取表达基因数在 2500 以下的细胞样本
    # 对数据进行标准化
    sc.pp.normalize_total(adata, target_sum=1e4)  ##标准化
    # target_sum=1e6相当于CPM标准化,如果设置为none，那么会用所有样品的median值来替代

    sc.pp.log1p(adata)  # log标准化后的值  将表达量用对数计算一遍，这样有利于绘图和展示

    sc.pl.highest_expr_genes(adata, n_top=20, log=True, save='-log1p.png')

    return adata

def highest_5(adata):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)  ##鉴定高变基因，参数均为默认值


    adata.raw = adata  # 将标准化后的数值存为.raw属性，方便后续分析

    adata = adata[:, adata.var.highly_variable]  #只留下高变基因进行后续分析
    sc.pp.regress_out(adata, ['total_counts'])  # 回归每个细胞的总计数的影响，将数据缩放到单位方差。

    sc.pp.scale(adata, max_value=10)  # scale数据
    sc.pl.highly_variable_genes(adata, save='.png')  # #可视化高变基因,绘制特异性基因散点图
    '''图黑色的点是高变基因，其他基因是灰色点。'''
    return adata

def pca_6(adata,expid = 'gene-LEMT1'):
    global ExampleID
    ExampleID = expid.replace('\n', '').replace('\r', '')
    # Principal component analysis
    sc.tl.pca(adata, svd_solver='arpack')  # PCA降维

    sc.pl.pca(adata, color=ExampleID, save='.png')  # PCA可视化

    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=40, save='.png')  # 绘制方差比,用碎石图来决定用多少个PC来进行临近细胞的计算
    # 方差贡献率图或者方差值图，便于观察PCA降维最佳值

    adata.write(results_file)
    return adata

def neighbor_7(adata,n_neighbors=30, n_pcs=20):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)  # 计算观察值的邻域图(搜索效率很大程度上依赖于UMAP),参数为默认值
    sc.tl.umap(adata)  # 将邻域图嵌入
    sc.tl.tsne(adata)  # 将邻域图嵌入
    return adata

def umap_all_8(adata):
    sc.pl.umap(adata, color=ExampleID, save='-example.png')
    return adata

def tsne_all_9(adata,resolution=0.4):
    sc.pl.tsne(adata, color=ExampleID, save='_tsne-example.png')

    # Embedding the neighborhood graph
    sc.tl.leiden(adata, resolution=resolution)  # 使用leiden算法对细胞进行聚类
    return adata

def paga_10(adata):
    sc.tl.paga(adata)
    sc.pl.paga(adata, save='.png')  # , plot=False  # remove `plot=False` if you want to see the coarse-grained graph
    return adata

def umap2_11(adata):
    sc.tl.umap(adata, init_pos='paga')
    # Clustering the neighborhood graph
    sc.pl.umap(adata, color=['leiden', ExampleID], save='.png')
    return adata

def tsne2_12(adata):
    sc.pl.tsne(adata, color=['leiden', ExampleID], save='_tsne.png')
    return adata

def rank_13(adata):
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')  # method can be 'wilcoxon','logreg'
    # 识别差异表达的基因，此函数将获取每组细胞，并将每组中每个基因的分布与不在该组中的所有其他细胞中的分布进行比较
    sc.pl.rank_genes_groups(adata, n_genes=15, sharey=False, save='.png')
    adata.write(results_file)

    adata.uns['leiden_colors'] = ['#FF8C00', '#e4a2ff', '#7FFFAA', '#4daf4a',
                                  '#3d3ff0', '#984ea3', '#fa0ff0', '#ff70a0',
                                  '#ffff33', '#e41a1c', '#a6f6f8', '#f781bf',
                                  '#999999', '#01110f']
    return adata

def fig_14(adata):
    from matplotlib.pyplot import rc_context
    ExampleID_list = list(pd.DataFrame(adata.uns['rank_genes_groups']['names']).iloc[0])
    ExampleID_list.append('leiden')

    with rc_context({'figure.figsize': (6, 6)}):
        sc.pl.umap(adata, color=ExampleID_list, s=50, frameon=False, ncols=3, save='_genegroup_umap.png')
        # 为了绘制缩放矫正的基因表达聚类图，需要使用 use_raw=False 参数。

    with rc_context({'figure.figsize': (6, 6)}):
        sc.pl.tsne(adata, color=ExampleID_list, s=50, frameon=False, ncols=3, save='_genegroup_tsne.png')

    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(adata, color='leiden', add_outline=True, legend_loc='on data',
                   legend_fontsize=12, legend_fontoutline=2, frameon=False,
                   title='clustering of cells', palette='Set1')

    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.tsne(adata, color='leiden', add_outline=True, legend_loc='on data',
                   legend_fontsize=12, legend_fontoutline=2, frameon=False,
                   title='clustering of cells', palette='Set1')
    return adata

def table_15(adata):
    # Finding marker genes
    table=pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)  # 之前不同组差异表达的基因的dataframe
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    return table,groups


    #result = adata.uns['rank_genes_groups']

def diff_gene_16(adata,n_genes=8,groups='0'):
    sc.pl.rank_genes_groups_violin(adata,groups=groups, n_genes=n_genes, save='.png')  # 想要某个组的高差异得分的基因更详细的视图  #group指定对比的细胞簇
    return adata

def diff_group_17(adata):
    sc.pl.violin(adata, [ExampleID], groupby='leiden', save='.png')  # 跨组比较某个基因
    return adata

def other_18(adata):
    sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', legend_fontsize='x-large', frameon=False,
               save='-label.png')

    sc.pl.tsne(adata, color='leiden', legend_loc='on data', title='', legend_fontsize='x-large', frameon=False,
               save='_tsne-label.png')

    marker_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).iloc[1,
                   :].tolist()  # 0_n,0_p~7_p的dataframe的第二行

    sc.pl.dotplot(adata, marker_genes, groupby='leiden', dendrogram=True, save='.png')

    sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', dendrogram=True, rotation=90, save='.png')

    adata.write(results_file,
                compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading

    adata.write_csvs(results_file[:-5], )
    return adata
    #return('done')


