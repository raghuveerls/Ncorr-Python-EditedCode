#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Tue Apr  2 17:35:38 2019

@author: raghuveer
'''
#%%
# =============================================================================
# User chooses the images to perform DIC on
# =============================================================================
def set_Images(UI_var_dict): #UI_var_dict stands for UserInterface_variables_dictionary
    from tkinter.filedialog import askopenfilenames
    import os
        
    pictures = askopenfilenames(title = 'Choose the files') #Get image names of all the images
    
    pic_list = list(pictures)   #Store image names with address in a list
    pic_names = []
    
    for file in pic_list:
       pic_names.append(os.path.basename(file)) #Obtain only the image names
        
    UI_var_dict['directory'] = os.path.dirname(pic_list[0]) #Obtain base directory and store in the dictionary
    (dummy, UI_var_dict['file_type']) = os.path.splitext(pic_names[0])  #Obtain image type and store in the dictionary
    
    i = 0
    for file in pic_names:
        os.rename(os.path.join(UI_var_dict['directory'], file), 
                  os.path.join(UI_var_dict['directory'], 'Image_' + str(i)+UI_var_dict['file_type'])) #Rename all images as "Image_#.file_type"
        i += 1
    
    UI_var_dict['file_total'] = str(len(pic_names)) #Store number of files in the dictionary

#%%
# =============================================================================
# User chooses the DIC parameters to use for analysis
# =============================================================================
def set_variables(root,UI_var_dict):
    import tkinter as tk
    from tkinter import messagebox
    import pickle
    import os
    
    direc = UI_var_dict['directory']    #Obtain directory in which images are stored
    if not os.path.exists(direc + '/DataFiles'):
            os.mkdir(direc + '/DataFiles')  #Create base directory to store the data files
    
    def set_var():
        # Set the variables in the UI variables dictionary
        UI_var_dict['scale_factor'] = str(scale.get()) #scalefactor for the dataset    
        UI_var_dict['subregion_radius'] = str(radius.get()) # subregion radius for the analysis
        UI_var_dict['threads'] = str(threads.get()) # # of threads to use
        UI_var_dict['subregion_shape'] = str(shape.get()) # subergion shape to use
        UI_var_dict['analysis_type'] = str(Type.get()) # small strain, large strain, or discontinuous strain
        pixelLen = pix.get()
        try:
            pixelLen = float(pixelLen)
            UI_var_dict['pixelLen'] = str(pixelLen) #pixels per mm
            setvar.destroy()
            with open(UI_var_dict['directory'] + '/DataFiles/UI_variables.pickle','wb') as handle:
                pickle.dump(UI_var_dict, handle, protocol = pickle.HIGHEST_PROTOCOL)    #Store dictionary variables in file
        except ValueError:  # If user doesn't input a number
            messagebox.showwarning('Error','Please enter numeric for Pixels/mm')
            return
     
    #Create window and set up labels for the different user input options   
    setvar = tk.Toplevel(root)
    tk.Label(setvar,text = 'Pixels/mm').grid(row=0)
    tk.Label(setvar,text = 'Scale Factor').grid(row=1)
    tk.Label(setvar,text = 'Subregion radius').grid(row=2)
    tk.Label(setvar,text = 'Threads').grid(row=3)
    tk.Label(setvar,text = 'Subregion shape').grid(row=4)
    tk.Label(setvar,text = 'Analysis type').grid(row=5)
    
    #Entry widget for inputting pixels/mm in the images
    pix = tk.Entry(setvar)
    pix.insert(1, float(UI_var_dict['pixelLen']))
    pix.grid(row = 0, column = 1)
    
    #Slider widget for inputting the scalefactor to be used
    scale = tk.Scale(setvar, from_ = 1, to=6, tickinterval=1, orient = tk.HORIZONTAL)
    scale.set(int(UI_var_dict['scale_factor']))
    scale.grid(row = 1, column = 1, columnspan = 2)
    
    #Slider widget for inputting the subset radius
    radius = tk.Scale(setvar, from_ = 5, to=50, tickinterval=10, orient = tk.HORIZONTAL)
    radius.set(int(UI_var_dict['subregion_radius']))
    radius.grid(row = 2, column = 1, columnspan = 2)
    
    #Slider widget for inputting the number of threads to use in computations
    threads = tk.Scale(setvar, from_ = 1, to=8, tickinterval=2, orient = tk.HORIZONTAL)
    threads.set(int(UI_var_dict['threads']))
    threads.grid(row = 3, column = 1, columnspan = 2)
    
    #Select subregion shape - Circle or Square 
    shape = tk.StringVar()
    shape.set(UI_var_dict['subregion_shape'])
    tk.Radiobutton(setvar, text = 'Circle', variable = shape, value ='circle').grid(row = 4, column = 1)        
    tk.Radiobutton(setvar, text = 'Square', variable = shape, value ='square').grid(row = 4, column = 2)
    
    #Select analysis type    - Large or Small strain 
    Type = tk.StringVar()
    Type.set(UI_var_dict['analysis_type'])
    tk.Radiobutton(setvar, text = 'Large strain', variable = Type, value ='large').grid(row = 5, column = 1)    
    tk.Radiobutton(setvar, text = 'Small strain', variable = Type, value ='small').grid(row = 5, column = 2)
    
    #Button to complete user input
    done = tk.Button(setvar,text = 'Done', command = set_var)
    done.grid(row = 6, column = 1, columnspan = 2)

#%%
# =============================================================================
# User chooses to process images or calculate displacements and gradients
# =============================================================================
def run(UI_var_dict, action):
    import subprocess as sp
    import tkinter as tk
    from tkinter import messagebox
    
    UI_var_dict['action'] = action  #Process images or calculate displacements?
    
    #Create and set up window to display output messages from Ncorr        
    root = tk.Toplevel()
    root.geometry('1000x400')
    scroll = tk.Scrollbar(root)
    scroll.pack(side = tk.RIGHT, fill = tk.Y)
    text = tk.Text(root)
    text.pack(fill = tk.BOTH, expand = 1)
    scroll.config(command = text.yview)
    text.config(yscrollcommand = scroll.set)
    text.insert(tk.END,'Starting\n')
    text.see(tk.END)
    
    #Run C++ executable
    out = sp.Popen(['./ncorr_test', UI_var_dict['action'], UI_var_dict['directory'], UI_var_dict['file_type'], UI_var_dict['file_total'],
            UI_var_dict['pixelLen'], UI_var_dict['scale_factor'], UI_var_dict['subregion_radius'], UI_var_dict['threads'], 
            UI_var_dict['subregion_shape'], UI_var_dict['analysis_type']], stdout = sp.PIPE, text = True, bufsize = 0)
    temp = out.stdout.readline() # read the first line
    while temp:
        if 'Saving' in temp or 'Processing' in temp or 'Changing perspective' in temp:
            root.title(temp)    #Update window title when necessary
        text.insert(tk.END, temp)
        text.see(tk.END)
        root.update_idletasks()
        temp = out.stdout.readline() # read output

    out.communicate()
    messagebox.showinfo('Complete',UI_var_dict['action'].capitalize()+' Finished!')
    root.destroy()   

#%%
# =============================================================================
# Create the required folders and store displacements and displacements gradients as .npy files
# =============================================================================
def dataprocess():
    import numpy as np
    from tkinter.filedialog import askdirectory
    from tkinter import messagebox
    import pickle
            
    direc = askdirectory()   #Base directory
    
    with open(direc + '/DataFiles/UI_variable.pickle','rb') as handle:
        UI_var_dict = pickle.load(handle)   #Get UI variables
    
    
    framecount = int(UI_var_dict['file_total']) #Total number of images
    
    #Create the required folders
    import os
    datafile_path = direc + '/DataFiles/'    #create file path for new directory for the data files

    def create_folders(perspective):
        if not os.path.exists(datafile_path + perspective):
            os.mkdir(datafile_path + perspective)                                   #Create base directory for each perspective
        if not os.path.exists(datafile_path + perspective + '/v_displacements'):
            os.mkdir(datafile_path + perspective + '/v_displacements')              #Create directory for storing v displacements
        if not os.path.exists(datafile_path + perspective + '/u_displacements'):
            os.mkdir(datafile_path + perspective + '/u_displacements')              #Create directory for storing u displacements
    
        if not os.path.exists(datafile_path + perspective + '/dux_dispgrad'):
            os.mkdir(datafile_path + perspective + '/dux_dispgrad')                 #Create directory for storing du/dx displacement gradient
        if not os.path.exists(datafile_path + perspective + '/duy_dispgrad'):
            os.mkdir(datafile_path + perspective + '/duy_dispgrad')                 #Create directory for storing du/dy displacement gradient
        if not os.path.exists(datafile_path + perspective + '/dvx_dispgrad'):
            os.mkdir(datafile_path + perspective + '/dvx_dispgrad')                 #Create directory for storing dv/dx displacement gradient
        if not os.path.exists(datafile_path + perspective + '/dvy_dispgrad'):
            os.mkdir(datafile_path + perspective + '/dvy_dispgrad')                 #Create directory for storing dv/dy displacement gradient

    
    create_folders('Lagrangian')   #Create folders for Lagrangian data
    create_folders('Eulerian')  #Create folders for Eulerian data
    
    run(UI_var_dict, 'dataprocess') #Run the C++ executable for get the data
    
    def numpy_disps(perspective): #function to get data from .csv files and collect it as .npy file
        
        disp_grad = []  #displacement gradients
        displacements = [] #displacements
        
        for i in range(framecount-1):
            
            #Obtain du/dx(dux), du/dy(duy), dv/dv(dvx), dv/dy(dvy) from their respective .csv files for each frame
            dux = np.genfromtxt(direc + '/DataFiles/' + perspective + '/dux_dispgrad/dux_array_' + perspective + '_frame' + str(i+1) + '.csv', delimiter = ',', dtype = None)
            duy = np.genfromtxt(direc + '/DataFiles/' + perspective + '/duy_dispgrad/duy_array_' + perspective + '_frame' + str(i+1) + '.csv', delimiter = ',', dtype = None)
            dvx = np.genfromtxt(direc + '/DataFiles/' + perspective + '/dvx_dispgrad/dvx_array_' + perspective + '_frame' + str(i+1) + '.csv', delimiter = ',', dtype = None)
            dvy = np.genfromtxt(direc + '/DataFiles/' + perspective + '/dvy_dispgrad/dvy_array_' + perspective + '_frame' + str(i+1) + '.csv', delimiter = ',', dtype = None)
            
            #Obtain u and v displacements from their respective .csv files for each frame
            u = np.genfromtxt(direc + '/DataFiles/' + perspective + '/u_displacements/u_array_' + perspective + '_frame' + str(i+1) + '.csv', delimiter = ',', dtype = None)
            v = np.genfromtxt(direc + '/DataFiles/' + perspective + '/v_displacements/v_array_' + perspective + '_frame' + str(i+1) + '.csv', delimiter = ',', dtype = None)
            
            #Append lists with latest frame information
            
            #displacement gradient[frame][du./dv. - 0/1][d.x/d.y - 0/1][y pixel][x pixel]
            disp_grad.append([[dux]+[duy]] + [[dvx] + [dvy]])
            
            #displacements[frame][u/v - 0/1][y pixel][x pixel]
            displacements.append([u] + [v])
        
        #Convert the displacements list to numpy arrays and save as .npy
        displacements = np.array(displacements)
        np.save(direc + '/DataFiles/' + perspective + '/displacements', displacements) 
        
        #Convert the displacement gradients list to numpy arrays and save as .npy
        disp_grad = np.array(disp_grad)
        np.save(direc + '/DataFiles/' + perspective + '/disp_grad', disp_grad)  
        
        messagebox.showinfo('Done', perspective + ' Data Processed!')
    
    numpy_disps('Lagrangian')   # Create .npy displacement and displacement gradient files for the Lagrangian data
    numpy_disps('Eulerian')     # Create .npy displacement and displacement gradient files for the Eulerian data
 
#%%
# =============================================================================
# Calculate Green or Almansi-Hamel strains
# =============================================================================
def calc_strains(perspective):
    import numpy as np
    from tkinter.filedialog import askdirectory
    from tkinter import messagebox
    import pickle
    
    askcsv = True
    askcsv = messagebox.askyesno('Save options','Avoid saving .csv files for each frame? '
                                 'Saving as .csv takes up a lot of space (.npy file will be saved by default))')
    
    direc = askdirectory()
    
    with open(direc + '/DataFiles/UI_variable.pickle','rb') as handle:
        UI_var_dict = pickle.load(handle)
    
    framecount = int(UI_var_dict['file_total'])
    
    if not askcsv:
        import os
        if not os.path.exists(direc + '/DataFiles/' + perspective + '/xx_Strains'):
            os.mkdir(direc + '/DataFiles/' + perspective + '/xx_Strains')
        if not os.path.exists(direc + '/DataFiles/' + perspective + '/xy_Strains'):
            os.mkdir(direc + '/DataFiles/' + perspective + '/xy_Strains')
        if not os.path.exists(direc + '/DataFiles/' + perspective + '/yy_Strains'):
            os.mkdir(direc + '/DataFiles/' + perspective + '/yy_Strains')
    
    try:
        disp_grad = np.load(direc + '/DataFiles/' + perspective + '/disp_grad.npy')
    except OSError:
        messagebox.showerror('Error','Process the data first!')
          
    Strain = []
    for i in range(framecount-1):
        dux = disp_grad[i][0][0]
        duy = disp_grad[i][0][1]
        dvx = disp_grad[i][1][0]
        dvy = disp_grad[i][1][1]
        if perspective == 'Lagrangian':
            xx = (0.5*(2*dux + dux**2 + dvx**2))
            xy = (0.5*(duy + dvx + dux*duy + dvx*dvy))
            yy = (0.5*(2*dvy + duy**2 + dvy**2))
        elif perspective == 'Eulerian':
            xx = (0.5*(2*dux - dux**2 - dvx**2))
            xy = (0.5*(duy + dvx - dux*duy - dvx*dvy))
            yy = (0.5*(2*dvy - duy**2 - dvy**2))
        
        #E[frame][Ex./Ey.(0/1)][E.x/E.y(0,1)][y pixel position][x pixel position]
        Strain.append([[xx]+[xy]] + [[xy]+[yy]])
        
        if not askcsv:
          #Save the Lagradngian Green Strain as .csv files
          np.savetxt(direc + '/DataFiles/' + perspective + '/xx_Strains/xx_Strain_Frame_' + str(i+1) + '.csv', xx, delimiter = ',')#Save Exx
          np.savetxt(direc + '/DataFiles/' + perspective + '/xy_Strains/xy_Strain_Frame_' + str(i+1) + '.csv', xy, delimiter = ',')#Save Exy
          np.savetxt(direc + '/DataFiles/' + perspective + '/yy_Strains/yy_Strain_Frame_' + str(i+1) + '.csv', yy, delimiter = ',')#Save Eyy
          
        if i == framecount-2:
            messagebox.showinfo('Done', perspective + ' Strains calculated!')
            
    Strain = np.array(Strain)
    np.save(direc + '/DataFiles/' + perspective + '/' + perspective + '_Strains', Strain)

#%%
def calc_defgrad():
    import numpy as np
    from tkinter.filedialog import askdirectory
    from tkinter import messagebox
    import pickle
    

    askcsv = True
    askcsv = messagebox.askyesno('Save options','Avoid saving .csv files for each frame? '
                                 'Saving as .csv takes up a lot of space (.npy file will be saved by default))')
    
    direc = askdirectory()
    with open(direc + '/DataFiles/UI_variable.pickle','rb') as handle:
        UI_var_dict = pickle.load(handle) 
    
    if not askcsv:
        import os
        if not os.path.exists(direc + '/DataFiles/Deformation_gradient'):
            os.mkdir(direc + '/DataFiles/Deformation_gradient')
        if not os.path.exists(direc + '/DataFiles/Deformation_gradient/F11'):
            os.mkdir(direc + '/DataFiles/Deformation_gradient/F11')
        if not os.path.exists(direc + '/DataFiles/Deformation_gradient/F12'):
            os.mkdir(direc + '/DataFiles/Deformation_gradient/F12')
        if not os.path.exists(direc + '/DataFiles/Deformation_gradient/F21'):
            os.mkdir(direc + '/DataFiles/Deformation_gradient/F21')
        if not os.path.exists(direc + '/DataFiles/Deformation_gradient/F22'):
            os.mkdir(direc + '/DataFiles/Deformation_gradient/F22')
    
    try:
        disp_grad_Lag = np.load(direc + '/DataFiles/Lagrangian/disp_grad.npy')
    except OSError:
        messagebox.showerror(title = 'File not found', message = 'Please make sure displacement gradients have been calculated!')
    
        
    
    framecount = int(UI_var_dict['file_total'])
    
    F = []
    for i in range(framecount-1):
        F11 = disp_grad_Lag[i][0][0] + 1    #dU/dX + 1
        F12 = disp_grad_Lag[i][0][1]        #dU/dY
        F21 = disp_grad_Lag[i][1][0]        #dV/dX
        F22 = disp_grad_Lag[i][1][1] + 1    #dV/dY + 1
        
        
        #F[frame][F1./F2.(0/1)[F.1/F.2(0/1)][y pixel position][x pixel position]
        F.append([[F11]+[F12]] + [[F21]+[F22]])
        
        if not askcsv:
            #Save the deformation gradients as .csv
            np.savetxt(direc +'/DataFiles/Deformation_gradient/F11/F11_frame' + str(i+1) + '.csv',F11,delimiter=',') #save F11
            np.savetxt(direc +'/DataFiles/Deformation_gradient/F12/F12_frame' + str(i+1) + '.csv',F12,delimiter=',') #save F11
            np.savetxt(direc +'/DataFiles/Deformation_gradient/F21/F21_frame' + str(i+1) + '.csv',F21,delimiter=',') #save F11
            np.savetxt(direc +'/DataFiles/Deformation_gradient/F22/F22_frame' + str(i+1) + '.csv',F22,delimiter=',') #save F11
        
        if i == framecount-2:
            messagebox.showinfo('Done','Deformation gradient calculated!')
            
    F = np.array(F)
    np.save(direc + '/DataFiles/F', F)
   
#%%
def dispvid_process():
    import numpy as np
    from PIL import Image
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
    from tkinter.filedialog import askdirectory
    import pickle
    from tkinter import messagebox
    
    direc = askdirectory()
    
    with open(direc + '/DataFiles/UI_variable.pickle','rb') as handle:
        UI_var_dict = pickle.load(handle)
        
    framecount = int(UI_var_dict['file_total'])
    file_type = UI_var_dict['file_type']
    scalefactor = int(UI_var_dict['scale_factor'])
    
    img = []
    for i in range(framecount):
        img.append(Image.open(direc + '/Image_' + str(i) + file_type))
    
    imgwidth, imgheight = img[0].size    
    
    def makevideo(perspective):
        displacements = np.load(direc + '/DataFiles/' + perspective + '/displacements.npy')
    
        roi_limits = scalefactor*np.genfromtxt(direc + '/DataFiles/' + perspective + '/roi_limits.csv', delimiter = ',', dtype = None)

        min_u = []    
        max_u = []
        min_v = []
        max_v = []
                
        for i in range(framecount-1):
            
            min_u.append(np.nanmin(displacements[i][0]))
            max_u.append(np.nanmax(displacements[i][0]))
            min_v.append(np.nanmin(displacements[i][1]))
            max_v.append(np.nanmax(displacements[i][1]))
        
        min_u = min(min_u)
        max_u = max(max_u)
        min_v = min(min_v)
        max_v = max(max_v)

        
        def make_animation(vmin, vmax, disp_type):
            if disp_type == 'u':
                array_pos = 0
            elif disp_type == 'v':
                array_pos = 1
                    
            dpi = 96
            fig, ax = plt.subplots(figsize = (imgwidth/dpi, imgheight/dpi), dpi = dpi)
            ax.set_xlim(0,imgwidth)
            ax.set_ylim(imgheight,0)
            if perspective == 'Lagrangian':
                start_frame = 0
            elif perspective == 'Eulerian':
                start_frame = 1
            imRAW = ax.imshow(img[start_frame], extent = [0, imgwidth, imgheight, 0], cmap = 'gray')
            imDIC = ax.imshow(displacements[0][array_pos], extent = [roi_limits[0][2], roi_limits[0][3], roi_limits[0][1], roi_limits[0][0]], 
                              vmin = vmin, vmax = vmax, alpha = 0.75, cmap = 'jet')
            ax.tick_params(axis = 'both', which = 'both', bottom = False, left = False, labelbottom = False, labelleft = False)
            ax.axis('off')
            fig.colorbar(imDIC)
            #fig.subplots_adjust(left = 0.0, right = 1.0, bottom = 0.0, top = 1.0, wspace = 0.0, hspace = 0.0)
            
            def animate(i):
                imDIC.set_data(displacements[i][array_pos][:,:-1])
                imDIC.set_extent([roi_limits[i][2], roi_limits[i][3], roi_limits[i][1], roi_limits[i][0]])
                if perspective == 'Lagrangian':
                    imRAW.set_data(img[0])
                elif perspective == 'Eulerian':
                    imRAW.set_data(img[i+1])
        
            ani = animation.FuncAnimation(fig, animate, frames = framecount-1, interval = 50)
            
            ani.save(direc + '/DataFiles/' + perspective + '/' + disp_type + '_displacement_' + perspective + '.mp4', dpi = dpi)
            
            messagebox.showinfo('Done','Completed ' + disp_type + ' displacement video in ' + perspective + ' coordinates!')
            plt.close(fig)   
        
        make_animation(min_u, max_u, 'u')
        make_animation(min_v, max_v, 'v')
        
    makevideo('Lagrangian')
    makevideo('Eulerian')
    
#%%
def strainvid_process(perspective):
    import numpy as np
    from PIL import Image
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
    from tkinter.filedialog import askdirectory
    import pickle
    from tkinter import messagebox
    
    direc = askdirectory()
    
    with open(direc + '/DataFiles/UI_variable.pickle','rb') as handle:
        UI_var_dict = pickle.load(handle)
        
    framecount = int(UI_var_dict['file_total'])
    file_type = UI_var_dict['file_type']
    scalefactor = int(UI_var_dict['scale_factor'])
    
    img = []
    for i in range(framecount):
        img.append(Image.open(direc + '/Image_' + str(i) + file_type))
    
    imgwidth, imgheight = img[0].size    
    
    strains = np.load(direc + '/DataFiles/' + perspective + '/' + perspective + '_Strains.npy')

    roi_limits = scalefactor*np.genfromtxt(direc + '/DataFiles/' + perspective + '/roi_limits.csv', delimiter = ',', dtype = None)

    min_xx = []    
    max_xx = []
    min_xy = []
    max_xy = []
    min_yy = []
    max_yy = []
            
    for i in range(framecount-1):

        min_xx.append(np.nanmin(strains[i][0][0]))
        max_xx.append(np.nanmax(strains[i][0][0]))
        min_xy.append(np.nanmin(strains[i][0][1]))
        max_xy.append(np.nanmax(strains[i][0][1]))
        min_yy.append(np.nanmin(strains[i][1][1]))
        max_yy.append(np.nanmax(strains[i][1][1]))
        
    min_xx = min(min_xx)
    max_xx = max(max_xx)
    min_xy = min(min_xy)
    max_xy = max(max_xy)
    min_yy = min(min_yy)
    max_yy = max(max_yy)    
    
    def make_animation(vmin, vmax, strain_type):
        if strain_type == 'xx':
            array_pos1 = 0
            array_pos2 = 0
        elif strain_type == 'xy':
            array_pos1 = 0
            array_pos2 = 1
        elif strain_type == 'yy':
            array_pos1 = 1
            array_pos2 = 1
                
        dpi = 96
        fig, ax = plt.subplots(figsize = (imgwidth/dpi, imgheight/dpi), dpi = dpi)
        ax.set_xlim(0,imgwidth)
        ax.set_ylim(imgheight,0)
        if perspective == 'Lagrangian':
            start_frame = 0
        elif perspective == 'Eulerian':
            start_frame = 1
        imRAW = ax.imshow(img[start_frame], extent = [0, imgwidth, imgheight, 0], cmap = 'gray')
        imDIC = ax.imshow(strains[0][array_pos1][array_pos2], extent = [roi_limits[0][2], roi_limits[0][3], roi_limits[0][1], roi_limits[0][0]], 
                          vmin = vmin, vmax = vmax, alpha = 0.75, cmap = 'jet')
        ax.tick_params(axis = 'both', which = 'both', bottom = False, left = False, labelbottom = False, labelleft = False)
        ax.axis('off')
        fig.colorbar(imDIC)
        
        def animate(i):
            imDIC.set_data(strains[i][array_pos1][array_pos2][:,:-1])
            imDIC.set_extent([roi_limits[i][2], roi_limits[i][3], roi_limits[i][1], roi_limits[i][0]])
            if perspective == 'Lagrangian':
                imRAW.set_data(img[0])
            elif perspective == 'Eulerian':
                imRAW.set_data(img[i+1])
    
        ani = animation.FuncAnimation(fig, animate, frames = framecount-1, interval = 50)
        
        ani.save(direc + '/DataFiles/' + perspective + '/' + strain_type + '_strain_' + perspective + '.mp4', dpi = dpi)
        
        messagebox.showinfo('Done','Completed ' + strain_type + ' strain video in ' + perspective + ' coordinates!')
        plt.close(fig)    
    
    make_animation(min_xx, max_xx, 'xx')
    make_animation(min_xy, max_xy, 'xy')
    make_animation(min_yy, max_yy, 'yy')

#%%
def polar_decomp():
    import numpy as np
    from numpy.linalg import qr
    from scipy.linalg import polar
    import pickle
    from tkinter.filedialog import askdirectory
    from tkinter import messagebox
    
    direc = askdirectory()
    
    with open(direc + '/DataFiles/UI_variable.pickle','rb') as handle:
        UI_var_dict = pickle.load(handle)
        
    F = np.load(direc + '/DataFiles/F.npy')
    
    framecount = int(UI_var_dict['file_total'])
    

    Ru = []
    Rv = []
    U = []
    V = []
    Q = []
    Rq = []
    rows, cols = F[0][0][0].shape
    
    for i in range(framecount-1):
        Ru_frame = []
        U_frame = []
        Rv_frame = []
        V_frame = []
        Q_frame = []
        Rq_frame = []
        
        for r in range(rows):
            Ru_row = []
            U_row = []
            Rv_row = []
            V_row = []
            Q_row = []
            Rq_row = []
            
            for c in range(cols):
                F_frame_pixel = np.array([[F[i][0][0][r][c], F[i][0][1][r][c]], 
                                          [F[i][1][0][r][c], F[i][1][1][r][c]]])
                
                Ru_pixel, U_pixel = polar(F_frame_pixel, side = 'right')        
                Ru_row = Ru_row + [Ru_pixel]
                U_row = U_row + [U_pixel]
                
                Rv_pixel, V_pixel = polar(F_frame_pixel, side = 'left')
                Rv_row = Rv_row + [Rv_pixel]
                V_row = V_row + [V_pixel]
                
                Q_pixel, Rq_pixel = qr(F_frame_pixel, mode = 'complete')
                Q_row = Q_row + [Q_pixel]
                Rq_row = Rq_row + [Rq_pixel]
            
            Ru_frame.append([Ru_row])
            U_frame.append([U_row])
            Rv_frame.append([Rv_row])
            V_frame.append([V_row])
            Q_frame.append([Q_row])
            Rq_frame.append([Rq_row])
            
        Ru.append([Ru_frame])
        U.append([U_frame])
        Rv.append([Rv_frame])
        V.append([V_frame])
        Q.append([Q_frame])
        Rq.append([Rq_frame])
    
    left_polar = [Ru] + [U]
    right_polar = [V] + [Rv]
    qr_decomp = [Q] + [Rq]
    
    left_polar = np.array(left_polar)
    np.save(direc + '/DataFiles/left_polar', left_polar)
    
    right_polar = np.array(right_polar)
    np.save(direc + '/DataFiles/right_polar', right_polar)
    
    qr_decomp = np.array(qr_decomp)
    np.save(direc + '/DataFiles/qr_decomp', qr_decomp)
    
    messagebox.showinfo('Done', 'Completed the decompositions!')
        
    
#%%    
def set_action(UI_var_dict, action = None):
    pass

