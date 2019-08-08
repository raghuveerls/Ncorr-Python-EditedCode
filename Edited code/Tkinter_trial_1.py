#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:13:17 2019

@author: raghuveer
"""
import tkinter as tk
from ncorr_functions import run, set_variables, set_Images, dataprocess, calc_defgrad, calc_strains, dispvid_process, strainvid_process, polar_decomp#, set_action

def printing(var):
    print(var)
  
UI_variables = {'action':'1',                   #1
                'directory':'2',                #2
                'file_type':'3',                #3
                'file_total' : '4',             #4
                'pixelLen' : '1.0',               #5
                'scale_factor' : '3',           #6
                'subregion_radius' : '20',      #7
                'threads' : '4',                #8
                'subregion_shape' : 'circle',   #9
                'analysis_type' : 'small',      #10
                }

   
root = tk.Tk()
root.title('Ncorr for Python')
root.geometry("500x500")
menubar = tk.Menu(root)

filemenu = tk.Menu(menubar, tearoff=0)
menubar.add_cascade(label='Image Processing', menu=filemenu)
filemenu.add_command(label='Print UI variables', command=(lambda: printing(UI_variables)))
filemenu.add_command(label='Set Images',command=(lambda: set_Images(UI_variables)))
filemenu.add_command(label='Set Variables',command=(lambda: set_variables(root, UI_variables)))
filemenu.add_command(label='Process Images',command=(lambda: run(UI_variables, 'calculate')))
filemenu.add_separator()
filemenu.add_command(label="Exit", command=root.destroy)


dataprocmenu = tk.Menu(menubar, tearoff=0)
menubar.add_cascade(label='Data Processing', menu=dataprocmenu)
dispproc = tk.Menu(dataprocmenu, tearoff=0)
dataprocmenu.add_command(label='Process the data',command = (lambda: dataprocess()))
postproc = tk.Menu(dataprocmenu, tearoff =0)
dataprocmenu.add_cascade(label = 'Other', menu=postproc)
postproc.add_command(label='Lagrangian strain',command=(lambda: calc_strains('Lagrangian')))
postproc.add_command(label='Eulerian strain',command=(lambda: calc_strains('Eulerian')))
postproc.add_command(label='Deformation gradient',command=(lambda: calc_defgrad()))
postproc.add_command(label='Polar Decomposition',command=(lambda: polar_decomp()))

videomenu = tk.Menu(menubar, tearoff = 0)
menubar.add_cascade(label='Video Processing', menu=videomenu)
videomenu.add_command(label='Displacements', command = (lambda: dispvid_process()))
videomenu.add_command(label='Lagrangian Strain', command = (lambda: strainvid_process('Lagrangian')))
videomenu.add_command(label='Eulerian Strain', command = (lambda: strainvid_process('Eulerian')))

root.config(menu=menubar)
root.mainloop()