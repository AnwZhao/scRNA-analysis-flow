from tkinter import *
from tkinter.ttk import *
import tkinter as tk
import tkinter.ttk
from tkinter import filedialog
from PIL import Image, ImageTk
import os
from tkinter import Label
import scRNA_sep
import tkinter.messagebox


# 打开指定的图片文件，缩放至指定尺寸
def get_image(filename, width, height):
    im = Image.open(filename).resize((width,height))
    return ImageTk.PhotoImage(im)

def txt_read(files):
    txt_dict = {}
    fopen = open(files)
    for line in fopen.readlines():
        line = str(line).replace("\n","")
        txt_dict[line.split(' ',1)[0]] = line.split(' ',1)[1]
        #split（）函数用法，逗号前面是以什么来分割，后面是分割成n+1个部分，且以数组形式从0开始
    fopen.close()
    return txt_dict

txt_dict = txt_read('explain.txt')


class scRNA_App():
    def __init__(self,master):

        self.root = master
        self.root.title("scRNA-seq Task Progress")
        self.root.geometry("800x600")

        # 创建画布，设置要显示的图片，把画布添加至应用程序窗口
        canvas_root = tkinter.Canvas(self.root, width=800, height=600)
        im_root = get_image('窗口图片.jpg', 800, 600)
        canvas_root.create_image(400, 300, image=im_root)
        canvas_root.pack()

        # 初始化任务进度条
        self.progress = tk.DoubleVar()
        self.progress.set(0.0)

        # 创建label
        self.label = tk.Label(self.root, text="")
        #self.label.pack(padx=10, pady=10)
        self.label.place(x=30,y=10)

        global current_path
        current_path = os.getcwd()

        # 创建按钮并添加点击事件
        self.button_select = tk.Button(self.root, text="Select Folder", command=self.select_folder, fg = "Gold", bg = "#006400") #选择文件夹按钮
        #self.button_select.pack(padx=10, pady=10)
        self.button_select.place(x=30,y=40)

        global adata

        self.button1 = tk.Button(self.root, text="1.read data", command=self.do_task1,width=12)
        self.button1.place(x=100,y=100)
        self.button2 = tk.Button(self.root, text="2.find top", command=self.do_task2,width=12)
        self.button2.place(x=200,y=100)
        self.button3 = tk.Button(self.root, text="3.filt data", command=self.do_task3,width=12)
        self.button3.place(x=300,y=100)
        self.button4 = tk.Button(self.root, text="4.norm&log", command=self.do_task4,width=12)
        self.button4.place(x=400, y=100)
        self.button5 = tk.Button(self.root, text="5.highest", command=self.do_task5,width=12)
        self.button5.place(x=500, y=100)
        self.button6 = tk.Button(self.root, text="6.pca process", command=self.do_task6,width=12)
        self.button6.place(x=600, y=100)

        self.button7 = tk.Button(self.root, text="7.find neigh", command=self.do_task7,width=12)
        self.button7.place(x=100, y=160)
        self.button8 = tk.Button(self.root, text="8.umap", command=self.do_task8,width=12)
        self.button8.place(x=200, y=160)
        self.button9 = tk.Button(self.root, text="9.tsne", command=self.do_task9,width=12)
        self.button9.place(x=300, y=160)
        self.button10 = tk.Button(self.root, text="10.paga", command=self.do_task10,width=12)
        self.button10.place(x=400, y=160)
        self.button11 = tk.Button(self.root, text="11.umap2", command=self.do_task11,width=12)
        self.button11.place(x=500, y=160)
        self.button12 = tk.Button(self.root, text="12.tsne2", command=self.do_task12,width=12)
        self.button12.place(x=600, y=160)

        self.button13 = tk.Button(self.root, text="13.rank gene", command=self.do_task13,width=12)
        self.button13.place(x=100, y=220)
        self.button14 = tk.Button(self.root, text="14.figure_all", command=self.do_task14,width=12)
        self.button14.place(x=200, y=220)
        self.button15 = tk.Button(self.root, text="15.marker table", command=self.do_task15,width=12)
        self.button15.place(x=300, y=220)
        self.button16 = tk.Button(self.root, text="16.diff_gene", command=self.do_task16,width=12)
        self.button16.place(x=400, y=220)
        self.button17 = tk.Button(self.root, text="17.diff_group", command=self.do_task17,width=12)
        self.button17.place(x=500, y=220)
        self.button18 = tk.Button(self.root, text="18.show other", command=self.do_task18,width=12)
        self.button18.place(x=600, y=220)

        self.show_select = tk.Button(self.root, text="Slide show result pictures", command=self.display_images, fg = "OrangeRed", bg = "Wheat")  # 幻灯片放映
        #self.show_select.pack(side = 'bottom')
        self.show_select.place(x=320,y=360)

        self.init()

        # 创建进度条并放置在窗口中
        self.progressbar = Progressbar(self.root, orient="horizontal", length=200, mode="determinate", variable=self.progress)
        #self.progressbar.pack(pady=20)
        self.progressbar.place(x=300,y=290)


        self.tip = tk.Label(self.root,
                            text="← Please select the data folder before pressing the data processing buttons.")
        self.tip.place(x=130,y=45)

        self.message = tk.Label(self.root,
                            text="Please follow the instructions.",anchor=SE)
        self.message.place(x=520, y=400)

        self.progress_num = tk.Label(self.root,
                                text="", anchor=SE)
        self.progress_num.place(x=505, y=290)

        self.back_main = tk.Button(self.root, text='back',font=('Helvetica', 15, 'bold'),  command=self.back)
        self.back_main.place(x=20, y=550)
        self.root.mainloop()

    def init(self):
        self.button1.configure(state=DISABLED)
        self.button2.configure(state=DISABLED)
        self.button3.configure(state=DISABLED)
        self.button4.configure(state=DISABLED)
        self.button5.configure(state=DISABLED)
        self.button6.configure(state=DISABLED)
        self.button7.configure(state=DISABLED)
        self.button8.configure(state=DISABLED)
        self.button9.configure(state=DISABLED)
        self.button10.configure(state=DISABLED)
        self.button11.configure(state=DISABLED)
        self.button12.configure(state=DISABLED)
        self.button13.configure(state=DISABLED)
        self.button14.configure(state=DISABLED)
        self.button15.configure(state=DISABLED)
        self.button16.configure(state=DISABLED)
        self.button17.configure(state=DISABLED)
        self.button18.configure(state=DISABLED)
        self.show_select.configure(state=DISABLED)

    def enable(self):
        self.button1.configure(state=NORMAL)
        self.button2.configure(state=NORMAL)
        self.button3.configure(state=NORMAL)
        self.button4.configure(state=NORMAL)
        self.button5.configure(state=NORMAL)
        self.button6.configure(state=NORMAL)
        self.button7.configure(state=NORMAL)
        self.button8.configure(state=NORMAL)
        self.button9.configure(state=NORMAL)
        self.button10.configure(state=NORMAL)
        self.button11.configure(state=NORMAL)
        self.button12.configure(state=NORMAL)
        self.button13.configure(state=NORMAL)
        self.button14.configure(state=NORMAL)
        self.button15.configure(state=NORMAL)
        self.button16.configure(state=NORMAL)
        self.button17.configure(state=NORMAL)
        self.button18.configure(state=NORMAL)
        self.show_select.configure(state=NORMAL)

    def do_task1(self):
        global adata
        self.button1.configure(text="busy...", state=DISABLED)
        adata= scRNA_sep.read_1(folder_path)
        self.progress.set(100/18*1)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button1.configure(text="1.read data\n完成/重启",fg='DarkCyan', state=NORMAL)

        self.message.config(text="After reading the data, press \n2.find topcalculate the counting\n scoreassigned to each gene in the\n cell.The gene with the highest average\n score in all cells is plotted \nas a box graph.",
                            justify='left')
        self.message.update()

        self.progress_num.config(text="1/18",justify='left')
        self.progress_num.update()

        self.n_top_message = tk.Label(self.root, text="show n top genes:", anchor=NW, justify='right')
        self.n_top_message.place(x=50, y=285)

        self.n_top = tk.Text(self.root, height=1.5, width=8)
        self.n_top.place(x=180, y=280)



    def do_task2(self):
        global adata
        self.button2.configure(text="busy...", state=DISABLED)
        if self.n_top.get("1.0", "end") != '\n' :
            n_top = int(self.n_top.get("1.0", "end"))
        else:
            n_top = 20
        adata= scRNA_sep.highest_2(adata,n_top)
        self.progress.set(100/18*2)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button2.configure(text="2.find top\n完成/重启",fg='DarkCyan', state=NORMAL)

        self.n_top.destroy()
        self.n_top_message.destroy()
        self.message.config(text='''According to the number and count 
of expressed genes,the abnormal values
of cells are filtered, and the minimum 
number of expressed genes required for 
cells to pass filtration is set to 200.
According to the number or count of 
cells, the minimum number of expression
cells required for gene filtration is 3.
and press to 3.filt data calculate 
quality control index''',
                            justify='left')
        self.message.update()
        
        self.progress_num.config(text="2/18", justify='left')
        self.progress_num.update()

        self.min_genes_message = tk.Label(self.root, text="filter min genes:", anchor=NW, justify='right')
        self.min_genes_message.place(x=50, y=285)

        self.min_genes = tk.Text(self.root, height=1.5, width=8)
        self.min_genes.place(x=180, y=280)

        self.min_cells_message = tk.Label(self.root, text="filter min cells:", anchor=NW, justify='right')
        self.min_cells_message.place(x=50, y=330)

        self.min_cells = tk.Text(self.root, height=1.5, width=8)
        self.min_cells.place(x=180, y=330)

    def do_task3(self):
        global adata
        self.button3.configure(text="busy...", state=DISABLED)
        if self.min_genes.get("1.0", "end") != '\n' :
            min_genes = int(self.min_genes.get("1.0", "end"))
        else:
            min_genes = 200
        if self.min_cells.get("1.0", "end") != '\n':
            min_cells = int(self.min_cells.get("1.0", "end"))
        else:
            min_cells = 3
        adata= scRNA_sep.filter_3(adata,min_genes,min_cells)
        self.progress.set(100/18*3)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button3.configure(text="3.filter\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(text='Obtain cell samples with the number\nof expressed genes below 2500,\nand standardize and log the data.\npress 4.norm&log',
                            justify='left')
        self.message.update()
        self.progress_num.config(text="3/18", justify='left')
        self.progress_num.update()

        self.min_genes.destroy()
        self.min_cells.destroy()
        self.min_genes_message.destroy()
        self.min_cells_message.destroy()

    def do_task4(self):
        global adata
        self.button4.configure(text="busy...", state=DISABLED)
        adata= scRNA_sep.normlog_4(adata)
        self.progress.set(100/18*4)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button4.configure(text="4.norm&log\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Visualize the hypervariable genes,\ndraw the scatter plot of specificgenes,\nand leave only the hypervariable\ngenes for subsequent analysis,\nand scale the data to unit\nvariance.press 5.highest.',
            justify='left')
        self.message.update()
        self.progress_num.config(text="4/18", justify='left')
        self.progress_num.update()

    def do_task5(self):
        global adata
        self.button5.configure(text="busy...", state=DISABLED)
        adata= scRNA_sep.highest_5(adata)
        self.progress.set(100/18*5)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button5.configure(text="5.highest\n完成/重启",fg='DarkCyan', state=NORMAL)
        print("3")
        self.message.config(
            text='The dimension of PCA is reduced\nand visualized, and the variance\ndiagram is drawn, which is\nconvenient for observing the optimal\nvalue of PCA dimension reduction.\npress 6.pca process.',
            justify='left')
        self.message.update()
        self.progress_num.config(text="5/18", justify='left')
        self.progress_num.update()

        self.expid_message = tk.Label(self.root, text="distribution of which gene:", anchor=NW, justify='right')
        self.expid_message.place(x=50, y=285)

        self.expid = tk.Text(self.root, height=1.5, width=15)
        self.expid.place(x=150, y=320)

    def do_task6(self):
        global adata
        self.button6.configure(text="busy...", state=DISABLED)
        if self.expid.get("1.0", "end") != '\n' :
            expid = self.expid.get("1.0", "end")
        else:
            expid = 'gene-LEMT1'
        adata= scRNA_sep.pca_6(adata,expid)
        self.progress.set(100/18*6)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button6.configure(text="6.pca process\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Embedding neighborhood graph.\nit will take a long time.\npress 7.find neigh',
            justify='left')
        self.message.update()
        self.progress_num.config(text="6/18", justify='left')
        self.progress_num.update()
        self.expid_message.destroy()
        self.expid.destroy()

        self.n_neighbors_message = tk.Label(self.root, text="n_neighbors for FindNeighbors:", anchor=NW, justify='right')
        self.n_neighbors_message.place(x=20, y=285)

        self.n_neighbors = tk.Text(self.root, height=1.5, width=8)
        self.n_neighbors.place(x=220, y=280)

        self.n_pcs_message = tk.Label(self.root, text="n_pcs for FindNeighbors:", anchor=NW, justify='right')
        self.n_pcs_message.place(x=20, y=330)

        self.n_pcs = tk.Text(self.root, height=1.5, width=8)
        self.n_pcs.place(x=220, y=330)

    def do_task7(self):
        global adata
        self.button7.configure(text="busy...", state=DISABLED)
        if self.n_neighbors.get("1.0", "end") != '\n' :
            n_neighbors = int(self.n_neighbors.get("1.0", "end"))
        else:
            n_neighbors = 50
        if self.n_pcs.get("1.0", "end") != '\n':
            n_pcs = int(self.n_pcs.get("1.0", "end"))
        else:
            n_pcs = 50
        adata= scRNA_sep.neighbor_7(adata,n_neighbors,n_pcs)
        self.progress.set(100/18*7)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button7.configure(text="7.find neigh\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Draw umap dimension reduction diagram.\npress 8.umap',
            justify='left')
        self.message.update()
        self.progress_num.config(text="7/18", justify='left')
        self.progress_num.update()

        self.n_neighbors.destroy()
        self.n_neighbors_message.destroy()
        self.n_pcs.destroy()
        self.n_pcs_message.destroy()

    def do_task8(self):
        global adata
        self.button8.configure(text="busy...", state=DISABLED)
        adata= scRNA_sep.umap_all_8(adata)
        self.progress.set(100/18*8)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button8.configure(text="8.umap\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Draw tsne dimension reduction diagram.\nand Using leiden algorithm to cluster\ncells.press 9.tsne',
            justify='left')
        self.message.update()
        self.progress_num.config(text="8/18", justify='left')
        self.progress_num.update()

        self.resolution_message = tk.Label(self.root, text="n_neighbors for leiden:", anchor=NW,
                                            justify='right')
        self.resolution_message.place(x=50, y=285)

        self.resolution = tk.Text(self.root, height=1.5, width=8)
        self.resolution.place(x=220, y=280)

    def do_task9(self):
        global adata
        self.button9.configure(text="busy...", state=DISABLED)

        if self.resolution.get("1.0", "end") != '\n' :
            resolution = int(self.resolution.get("1.0", "end"))
        else:
            resolution = 1.4
        adata= scRNA_sep.tsne_all_9(adata,resolution)
        self.progress.set(100/18*9)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button9.configure(text="9.tsne\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Draw the trace map of paga cells.\npress 10.paga',
            justify='left')
        self.message.update()
        self.progress_num.config(text="9/18", justify='left')
        self.progress_num.update()

        self.resolution_message.destroy()
        self.resolution.destroy()

    def do_task10(self):
        global adata
        self.button10.configure(text="busy...", state=DISABLED)
        adata= scRNA_sep.paga_10(adata)
        self.progress.set(100/18*10)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button10.configure(text="10.paga\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Draw the distribution of a\ngene in leiden dimension reduction\ndiagram in umap form.\npress 11.umap2',
            justify='left')
        self.message.update()
        self.progress_num.config(text="10/18", justify='left')
        self.progress_num.update()

    def do_task11(self):
        global adata
        self.button11.configure(text="busy...", state=DISABLED)
        adata= scRNA_sep.umap2_11(adata)
        self.progress.set(100/18*11)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button11.configure(text="11.umap2\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Draw the distribution of a\ngene in leiden dimension reduction\ndiagram in tsne form.\npress 12.tsne2',
            justify='left')
        self.message.update()
        self.progress_num.config(text="11/18", justify='left')
        self.progress_num.update()

    def do_task12(self):
        global adata
        self.button12.configure(text="busy...", state=DISABLED)
        adata= scRNA_sep.tsne2_12(adata)
        self.progress.set(100/18*12)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button12.configure(text="12.tsne2\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Identify differentially expressed\ngenes. This function will obtain\neach group of cells and compare\nthe distribution of each gene\nin each group with that in\nall other cells that are not in the\ngroup.press 13.rank gene',
            justify='left')
        self.message.update()
        self.progress_num.config(text="12/18", justify='left')
        self.progress_num.update()

    def do_task13(self):
        global adata
        self.button13.configure(text="busy...", state=DISABLED)
        adata= scRNA_sep.rank_13(adata)
        self.progress.set(100/18*13)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button13.configure(text="13.rank gene\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Visualize the distribution\nand expression of the highest\nvariable genes in each group\nin the cluster map.\npress 14.figure_all',
            justify='left')
        self.message.update()
        self.progress_num.config(text="13/18", justify='left')
        self.progress_num.update()

    def do_task14(self):
        global adata
        self.button14.configure(text="busy...", state=DISABLED)
        adata= scRNA_sep.fig_14(adata)
        self.progress.set(100/18*14)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button14.configure(text="14.figure_all\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Dataframe of previously differentially\nexpressed genes in different\ngroups.press 15.marker table',
            justify='left')
        self.message.update()
        self.progress_num.config(text="14/18", justify='left')
        self.progress_num.update()

    def do_task15(self):

        self.button15.configure(text="busy...", state=DISABLED)
        scRNA_sep.table_15(adata)
        self.progress.set(100/18*15)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button15.configure(text="15.marker table\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Dataframe of previously differentially\nexpressed genes in different\ngroups.press 16.diff_gene',
            justify='left')
        self.message.update()
        self.progress_num.config(text="15/18", justify='left')
        self.progress_num.update()

        self.groups_message = tk.Label(self.root, text="which group to show:", anchor=NW, justify='right')
        self.groups_message.place(x=30, y=285)

        self.groups = tk.Text(self.root, height=1.5, width=8)
        self.groups.place(x=190, y=280)

        self.n_genes_message = tk.Label(self.root, text="How many genes to show:", anchor=NW, justify='right')
        self.n_genes_message.place(x=30, y=330)

        self.n_genes = tk.Text(self.root, height=1.5, width=8)
        self.n_genes.place(x=190, y=330)

    def do_task16(self):
        global adata
        self.button16.configure(text="busy...", state=DISABLED)
        if self.groups.get("1.0", "end") != '\n':
            groups = self.groups.get("1.0", "end")
        else:
            groups = '0'

        if self.n_genes.get("1.0", "end") != '\n':
            n_genes = int(self.n_genes.get("1.0", "end"))
        else:
            n_genes = 8

        adata= scRNA_sep.diff_gene_16(adata,n_genes,groups)
        self.progress.set(100/18*16)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button16.configure(text="16.diff_gene\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Want a more detailed view\nof a group of genes with high\ndifferential scores.press\n17.diff_group',
            justify='left')
        self.message.update()
        self.progress_num.config(text="16/18", justify='left')
        self.progress_num.update()

        self.groups_message.destroy()
        self.groups.destroy()
        self.n_genes_message.destroy()
        self.n_genes.destroy()

    def do_task17(self):
        global adata
        self.button17.configure(text="busy...", state=DISABLED)
        adata= scRNA_sep.diff_group_17(adata)
        self.progress.set(100/18*17)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button17.configure(text="17.diff_group\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='View other forms of cluster\nmaps and violin maps of marker\ngene expression in each\ngroup.press 18.show other',
            justify='left')
        self.message.update()
        self.progress_num.config(text="17/18", justify='left')
        self.progress_num.update()

    def do_task18(self):
        global adata
        self.button18.configure(text="busy...", state=DISABLED)
        adata= scRNA_sep.other_18(adata)
        self.progress.set(100/18*18)
        self.root.update()
        tkinter.messagebox.showinfo('System Prompt', 'Task completed.')
        self.button18.configure(text="18.show other\n完成/重启",fg='DarkCyan', state=NORMAL)
        self.message.config(
            text='Finish!\nNow you can press the button\nto play all the generated pictures\nin the form of slides, and\nview all the pictures in the\nfiguers folder of the selected\ndirectory.',
            justify='left')
        self.message.update()
        self.progress_num.config(text="18/18", justify='left')
        self.progress_num.update()


    def select_folder(self):
        global folder_path
        folder_path = filedialog.askdirectory()
        self.label.config(text="Selected Folder: " + folder_path)
        self.tip.config(text="← data folder selected")
        #self.tip = tk.Label(self.root, text="← data folder selected")
        #self.tip.place(x=130, y=45)
        self.tip.update()
        self.enable()
        self.message.config(text="Now you can take the first step:\nread the data of the selected folder.\npress 1.read data",justify='left')
        self.message.update()

    def display_images(self):
        self.top = Toplevel()
        self.top.title('幻灯片展示')
        self.top.geometry("800x600")
        # 创建画布，设置要显示的图片，把画布添加至应用程序窗口
        canvas_top = tkinter.Canvas(self.top, width=800, height=600)
        im_top = get_image(current_path+'/窗口图片.jpg', 800, 600)
        canvas_top.create_image(400, 300, image=im_top)
        canvas_top.pack()
        images = []

        #upDir=os.path.abspath(os.path.join(os.path.dirname(folder_path), os.path.pardir))
        upDir=os.path.abspath (os.path.join (folder_path, "../"))
        print(upDir)
        folder_path_figures = upDir + r'\figures'
        print(folder_path_figures)
        for filename in os.listdir(folder_path_figures):
            if filename.endswith(('.jpg', '.png', '.jpeg')):
                image = Image.open(os.path.join(folder_path_figures, filename))
                image = image.resize((400, 400))
                images.append(ImageTk.PhotoImage(image))

        if images:
            index = 0

            def show_image():
                nonlocal index

                if index >= len(images):
                    #self.top.label.config(text="已经读取完毕")
                    #label(self.top, text="已经读取完毕")
                    widget = Label(self.top, text="已经读取完毕")
                    widget.config(bg='black', fg='yellow')

                    return
                #im_show = Label(self.top)
                #im_show.config(image=images[index])

                label_img = Label(self.top, image=images[index])
                label_img.configure(image=images[index])
                label_img.place(x=10, y=10)
                #label = tk.Label(self.top, image=images[index])
                #label.pack()

                widget = Label(self.top, text="This is the {} th picture.".format(index + 1),font=('Arial', 15))
                widget.config(bg='black', fg='yellow')
                widget.place(x=500, y=100)

                explain = Label(self.top, text=txt_dict[str(index)], height=15, width=35,wraplength=300,
                                font=('Arial', 12), justify='left')
                explain.config(fg="#006400")
                explain.place(x=450, y=150)



                #explain.update()


                index += 1

                self.top.after(4000, show_image)

            show_image()

        self.top.mainloop()

    def back(self):
        self.root.destroy()
        self.root = tk.Tk()
        from main_show import basedesk  #防止互相调用进入死循环
        basedesk(self.root)






#scRNA_App()