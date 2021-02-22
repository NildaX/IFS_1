from tkinter import filedialog
from tkinter import *
from tkinter import ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np
from astropy.io import fits
import matplotlib.colors
import cv2
from astropy.visualization.mpl_normalize import simple_norm
import math
from skimage import exposure
from astropy.wcs import WCS
import os
import warnings
from tkinter import messagebox as MessageBox
from astropy.stats import sigma_clip
from scipy import interpolate
import ctypes
from astropy.coordinates import SkyCoord
from PIL import Image
from PIL import ImageTk
import matplotlib.patches as patches
from tkinter import scrolledtext
import tkinter.font as tkFont
import platform

warnings.filterwarnings("ignore")


def onclick_(event):  
    global flag_explorer,flag_integrate_region2
    if flag_file==1:
        if flag_integrate_region == 0:
            if flag_explorer == 1:
                flag_explorer = 0
                varla10.set("Explorer OFF")
            else:
                flag_explorer = 1
                varla10.set("Explorer ON")
        else:
            if flag_integrate_region2 == 0:
                try:
                    cord_x = int(round(event.xdata))
                    cord_y = int(round(event.ydata))
                    draw_circle(int(round(event.xdata)), int(round(event.ydata)))
                    
                except Exception as e:
                    print(e)
               

    
def coordinates_(cord_x,cord_y):
    global tam_x,tam_y,ax0,canvas, arrlambda,header_file,wcs,spectrum,arr_w,arr_f,red_marks
    ID = (cord_y*tam_x)+cord_x
    ax0.cla()
    ax0.set_xlabel( 'Wavelength' )
    ax0.set_ylabel( 'Flux' )
    spectrum=data[:,cord_y,cord_x]
    spectrum=np.nan_to_num(spectrum) 
    for i in range(0,pixels):
       if (abs(spectrum[i])>1e30):
           spectrum[i]=0
    ax0.plot(arrlambda,spectrum,color='blue')
    if flag_wave == 1:
       ax0.set_xlim(xmin=varla5.get(),xmax=varla6.get())
    else:
       ax0.set_xlim(xmin=min_value_la,xmax=max_value_la)
    if flag_flux==1:
       ax0.set_ylim(ymin=varla3.get(),ymax=varla4.get())
       ax0.plot(arrlambda,res*varla4.get(),'--',color = 'orange')
       for i in red_marks:
           ax0.axvline(int(i),varla3.get(),varla4.get(),color='red')
       
    else:
       ax0.plot(arrlambda,res*np.amax(spectrum)*1.2,'--',color = 'orange')
       ax0.set_ylim(ymin=np.amin(spectrum),ymax=np.amax(spectrum)*1.2)
       for i in red_marks:
           ax0.axvline(int(i),min_value_da,max_value_da,color='red')
    canvas.draw()
    if wcs == 0:
       varla8.set("   [ %d , %d ]       %d        "%(cord_x,cord_y,ID))
    else:
       try: 
           celestial, spectral = wcs_header.pixel_to_world([cord_x], [cord_y], [1])  
           ra_1=celestial.ra.hms
           c = str(celestial.dec)
           c=c.replace("d"," d ")
           c=c.replace("m"," m ")
           c=c.replace("s"," s ")
           c=c.replace("["," ")
           c=c.replace("]"," ")
           degr = str(celestial.to_string('decimal'))
           degr = degr.replace(" ","              ")
           degr = degr.replace("["," ")
           degr = degr.replace("]"," ")
           degr = degr.replace("'"," ")
           varla8.set("   [ %d , %d ]       %d        %d h  %d m  %f s     %s     %s "%(cord_x,cord_y,ID,ra_1[0],ra_1[1],ra_1[2],c,degr))
       except Exception as e:
           celestial = wcs_header.pixel_to_world([cord_x], [cord_y], [1])    
           coo = SkyCoord(ra=celestial[0], dec=celestial[1],unit="deg")
           ra_1=coo.ra.hms
           c = str(coo.dec)
           c=c.replace("d"," d ")
           c=c.replace("m"," m ")
           c=c.replace("s"," s ")
           c=c.replace("["," ")
           c=c.replace("]"," ")
           completa= "   [ %d , %d ]       %d        %d h  %d m  %f s   %s  "%(cord_x,cord_y,ID,ra_1[0],ra_1[1],ra_1[2],c)
           completa_1 = completa+str(celestial[0])+str(celestial[1])
           completa_1 = completa_1.replace("deg"," ")
           completa_1 = completa_1.replace("["," ")
           completa_1 = completa_1.replace("]"," ")
           varla8.set(completa_1)
           
def create_txt():
    try:
        name_text = varla1.get()
        name_text = name_text.split('.')
        name_text = name_text[0]+'.spectrum_'+str(cir_x)+'_'+str(cir_y)+'.text'
        final_dir = file_dir.replace(varla1.get(),name_text)
        
        file = open(final_dir, "w")
        for i in range(len(arrlambda)):
            line = "   %d     %f     %f  \n"%(i+1,arrlambda[i],Integrated_spectrum[i])
            file.write(line)
        file.close()
        text_1 = (u" File created: %s \n "%(name_text))
        box_entry.configure(state='normal')
        box_entry.insert(INSERT, text_1)
        box_entry.see(END)
        box_entry.configure(state='disabled')
    except Exception as e:
        print(e)

def create_spectrum_fits():
    try:
        name_spectrum = varla1.get()
        name_spectrum = name_spectrum.split('.')
        name_spectrum = name_spectrum[0]+'.spectrum_'+str(cir_x)+'_'+str(cir_y)+'.fits'
        fits_s = fits.PrimaryHDU(Integrated_spectrum)
        final_dir = file_dir.replace(varla1.get(),name_spectrum)
        fits_s.writeto(final_dir)  
        text_1 = (u" File created: %s \n "%(name_spectrum))
        box_entry.configure(state='normal')
        box_entry.insert(INSERT, text_1)
        box_entry.see(END)
        box_entry.configure(state='disabled')
    except Exception as e:
        print(e)
    
def circle_fits():
    var = (int(radius_.get())*2)+1
    new_data = np.zeros([pixels,var,var])
    k = 0
    lim = 1
    cont = 0
    cont2 = int(radius_.get())
    n = var-1
    for i in range(int(var/2)+1):
        cont2 = int(radius_.get())-i
        cont=0
        for j in range(var):
            if j >=cont2 and cont<lim and k < len(integrated_x):
                if 0<=integrated_y[k]<tam_y and 0<=integrated_x[k]<tam_x:
                    new_data[:,n,j] = data[:,integrated_y[k],integrated_x[k]]
                k = k+1
                cont = cont+1    
        n = n-1
        lim = lim+2
    lim = lim-4
    v = 1
    while i<var-1:
        cont2 = int(radius_.get())-i+v
        cont=0
        for j in range(var):
            if j >= cont2 and k < len(integrated_x) and cont<lim:
                if 0<=integrated_y[k]<tam_y and 0<=integrated_x[k]<tam_x:
                    new_data[:,n,j] = data[:,integrated_y[k],integrated_x[k]]
                k = k+1
                cont = cont+1
        n = n-1
        i = i+1
        lim = lim-2
        v = v+2
    return new_data

def create_circle_fits():
    try:
        new_data = circle_fits()
        name_circle = varla1.get()
        name_circle = name_circle.split('.')
        name_circle = name_circle[0]+'.cirle.rscube_'+str(cir_x)+'_'+str(cir_y)+'.fits'
        final_dir = file_dir.replace(varla1.get(),name_circle) 
        r = int(radius_.get())
        v = (r*2)+1
        header_file['NAXIS1'] = v
        header_file['NAXIS2'] = v
        header_file['CRPIX1'] = int(v/2)
        header_file['CRPIX2'] = int(v/2)
        fits.writeto(final_dir,data=new_data,header=header_file)
        text_1 = (u" File created: %s \n "%(name_circle))
        box_entry.configure(state='normal')
        box_entry.insert(INSERT, text_1)
        box_entry.see(END)
        box_entry.configure(state='disabled')
    except Exception as e:
        print(e)
    
def create_files():
    global flag_create_fits
    if flag_file==1:
        if flag_integrate_region == 1 and flag_integrate_region2==1:
            if flag_create_fits == 0:
                create_txt()
                create_spectrum_fits()
                create_circle_fits()
                flag_create_fits = 1
            else:
                MessageBox.showerror("Error!","Please, select other region to create files")
            
            
        else:
            MessageBox.showerror("Error!","Please, first integer region")#revisar
    else:
        MessageBox.showerror("Error!","Please, first choose a file")
           
def get_coor(x,y,radius):
    pix_x = []
    pix_y = []   
    k = 1
    r = 1     
    num = ((radius*2)+1)
    for i in range(int(num/2)+1):
        k = 1
        while k <= r:
            dif = (k-1) - i
          #  if 0 < (x+dif) < tam_x and 0 < (y+radius-i) < tam_y:
            pix_x.append(x+dif)
            pix_y.append(y+radius-i)
            k = k+1
        r += 2
        
    r = ((radius*2)+1)-2
    i = i+1
    l = 0
    for j in range(int(num/2)):
        k = 1
        while k <= r:
       #     if 0 < (x-radius+k+l) < tam_x and 0 < (y+radius-i) < tam_y:
           pix_x.append(x-radius+k+l)
           pix_y.append(y+radius-i)
           k = k+1
            
        r = r-2  
        i = i+1
        l = l+1
        
    return pix_x,pix_y
    
def get_integrated_spectrum(pix_x,pix_y):
    global Integrated_spectrum 
    Integrated_spectrum = np.zeros(pixels)   
    for j in range(len(pix_x)):
        if 0<=pix_y[j]<tam_y and 0<=pix_x[j]<tam_x:
            esp=data[:,pix_y[j],pix_x[j]]
            esp=np.nan_to_num(esp)   
            for i in range(0,pixels):
                if (abs(esp[i])>1e30):
                    esp[i]=0
            Integrated_spectrum = Integrated_spectrum + esp
    
   
    
    
def draw_circle(cord_x,cord_y):
    global flag_integrate_region,flag_explorer,ax1,canvas2,flag_integrate_region2,cir_x,cir_y,integrated_x,integrated_y
    try:
        if int(radius_.get()) <= 0 or int(radius_.get()) >=tam_x/2:
            MessageBox.showwarning("Warning!","Please, enter a radius positive and logic")
            flag_integrate_region = 0
            flag_integrate_region2 = 0
            integrated_x = []
            integrated_y = []
            radius_.set(0)
            flag_explorer = 0
            varla10.set("Explorer OFF")
            label29.config(state=NORMAL)
            
        else:
             cir = patches.Circle((cord_x,cord_y),
                      int(radius_.get()),
                      edgecolor='red',
                      fill = False)
             ax1.add_patch(cir)
             canvas2.draw()
             
             cir_x=cord_x
             cir_y=cord_y
             flag_explorer = 0
             arr_x,arr_y=get_coor(cord_x,cord_y,int(radius_.get()))
             get_integrated_spectrum(arr_x,arr_y)
             varla10.set("Show integrated flux")
             text_1 = (u"\n Drawing circle with center %d,%d and radius of %d \n"%(cord_x,cord_y,radius_.get()))
             integrated_x = arr_x
             integrated_y = arr_y
             box_entry.configure(state='normal')
             box_entry.insert(INSERT, text_1)
             box_entry.see(END)
             box_entry.configure(state='disabled')     
             flag_integrate_region2 = 1        
             filters_(name_f)
    except Exception as e: 
           print(e)
           MessageBox.showerror("Error!","Please, enter numbers")  
           flag_integrate_region = 0
           flag_integrate_region2 = 0
           radius_.set(0)

 
def move_mouse(event):
    if flag_explorer == 1:
        try:
            cord_x = int(round(event.xdata))
            cord_y = int(round(event.ydata))
            coordinates_(int(round(event.xdata)), int(round(event.ydata)))
        except Exception as e:
            var = "not pixel in graph"
    
 
def create_label_offset():
    global x_ticks,y_ticks,x_ticks_l,y_ticks_l
    try:
        x_ticks=[]
        y_ticks=[]
        x_ticks_l=[]
        y_ticks_l=[]
        cdelt_1 = header_file['CDELT1']
        cdelt_2 = header_file['CDELT2']
        
            
            
        if cdelt_1 <= 0: 
            cdelt_1 = cdelt_1*-1
        if cdelt_2 <= 0: 
            cdelt_1 = cdelt_2*-1
            
        if cdelt_1 <= 0.01: 
            cdelt_1 = cdelt_1*100
            
        if cdelt_2 <= 0.01: 
            cdelt_2 = cdelt_2*100
            
        text_1 = (u"CDELT1 = %f  \n CDELT2 = %f \n"%(cdelt_1,cdelt_2) )
        box_entry.configure(state='normal')
        box_entry.insert(INSERT, text_1)
        box_entry.see(END)
        box_entry.configure(state='disabled')
        crpix_1 = header_file['CRPIX1']
        crpix_2 = header_file['CRPIX2']
        if crpix_1 <= 0 or crpix_1 > tam_x:
            crpix_1 = int(tam_x/2)
            text_1 = (u"  CRPIX1 entry was negative or greater than NAXIS1, using CRPIX1= NAXIS1/2 \n" )
            box_entry.configure(state='normal')
            box_entry.insert(INSERT, text_1)
            box_entry.see(END)
            box_entry.configure(state='disabled')
        else:
            text_1 = (u"  CRPIX1 = %f  \n "%(crpix_1) )
            box_entry.configure(state='normal')
            box_entry.insert(INSERT, text_1)
            box_entry.see(END)
            box_entry.configure(state='disabled')
            
        if crpix_2 <= 0 or crpix_2 > tam_y:
            crpix_2 = int(tam_y/2)
            text_1 = (u"CRPIX2 entry was negative or greater than NAXIS2, using CRPIX2= NAXIS2/2 \n" )
            box_entry.configure(state='normal')
            box_entry.insert(INSERT, text_1)
            box_entry.see(END)
            box_entry.configure(state='disabled')
        else:
            text_1 = (u"CRPIX2 = %f  \n "%(crpix_2) )
            box_entry.configure(state='normal')
            box_entry.insert(INSERT, text_1)
            box_entry.see(END)
            box_entry.configure(state='disabled')
            
        x_ticks= np.linspace(0,tam_x-1,tam_x)
        y_ticks= np.linspace(0,tam_y-1,tam_y)
        for s in range(tam_x):
            x_ticks_l.append(" ")
        for s in range(tam_y):
            y_ticks_l.append(" ")
                
        x_ticks_l[int(crpix_1)] = 0    
        y_ticks_l[int(crpix_2)] = 0

        
        
        num_a=(tam_x * cdelt_1)/tam_x
        num_b=(tam_y * cdelt_2)/tam_y
            
        s_1 = int(crpix_1+1)
        inter_1=int(tam_x/10)
        cont_1=0
        while s_1 < tam_x:
            if cont_1%inter_1 == 0 and cont_1 != 0:
                x_ticks_l[s_1]=int(cont_1*-1*num_a)
            cont_1=cont_1+1
            s_1 = s_1+1
        s_1 = int(crpix_1-1)
        cont_1=0;
        while s_1 > 0:
            if cont_1%inter_1 == 0 and cont_1 != 0:
                x_ticks_l[s_1]=int(cont_1*num_a)
            cont_1=cont_1+1
            s_1 = s_1-1    
        s_1 = int(crpix_2+1)
        inter_1=int(tam_y/10)
        cont_1=0
        while s_1 < tam_y:
            if cont_1%inter_1 == 0 and cont_1 != 0 :
                y_ticks_l[s_1]=int(cont_1*num_b)
            cont_1=cont_1+1
            s_1 = s_1+1
        s_1 = int(crpix_2-1)
        cont_1=0
        while s_1 > 0:
            if cont_1%inter_1 == 0 and cont_1 != 0:
                y_ticks_l[s_1]=int(cont_1*-1*num_b)
            cont_1=cont_1+1
            s_1 = s_1-1
        
    except Exception as e:    
        print(e)
        
        
def new_file():
    global flag_create_fits,flag_integrate_region,flag_integrate_region2,file_dir,box_entry,Integrated_spectrum,band_sticks,combo1,combo2,canvas2,f2,ax1,wcs_header,wcs,name,hi_data,data,header_file,datosflux,pixels,arrlambda,spectrum,crval,cdelt,crpix,tam_x,tam_y,imagen_final,imagen_final,frame7,varla5,varla6,flag_wave,flag_flux,red_marks,band,name_f,min_value_la,max_value_la,flag_band,array_data,dband,bar_,min_value_da,max_value_da,flag_file,combo1,combo2,integrated_x, integrated_y
    
    
    if flag_file == 1:
        hi_data.close()
        
 #   folder_selected = filedialog.askopenfile(mode="r")#, filetypes = (("fits files","*.fits"),("all files","*.*")))
    py_exts = r"*.fits  *.fits.gz *.fits.rar"
    folder_selected = filedialog.askopenfile(mode="r", filetypes = (("fits files",py_exts),("all files","*.*")))
    file_dir = os.path.abspath(folder_selected.name)
    name=folder_selected.name
    nombre2=name.split('/')
    n=nombre2[len(nombre2)-1]
    varla1.set(n)
    name = n
    print(file_dir)
    hi_data = fits.open(file_dir)
    data = hi_data[0].data
    
    min_value_da=np.amin(data)
    max_value_da=np.amax(data)
    hi_data.info()
    header_file = hi_data[0].header
    
    band_sticks.set(0)
    integrated_x = []
    integrated_y = []
    if flag_file == 1:
        box_entry.configure(state='normal')
        box_entry.delete("1.0","end")
        box_entry.configure(state='disabled')  
        
    text_1 = (u" \n IFS Explorer 3D cube spectra viewer\n" 
                  u"Move mouse over the mosaic to plot the spectra \n"
                  u"Click for ON/OFF the explorer  \n")
    box_entry.configure(state='normal')
    box_entry.insert(INSERT, text_1)
    box_entry.see(END)
    box_entry.configure(state='disabled')
    
    
    #---
      
    flag_file=1
    try:
        crval  = header_file['CRVAL3']
        cdelt  = header_file['CDELT3']
        crpix  = header_file['CRPIX3']
        tam_x  = header_file['NAXIS1']
        tam_y  = header_file['NAXIS2']
        pixels = header_file['NAXIS3']
        text_1 = (u"  The dimensions of the cube are: %d x %d x %d \n"%(tam_x,tam_y,pixels) )
        box_entry.configure(state='normal')
        box_entry.insert(INSERT, text_1)
        box_entry.see(END)
        box_entry.configure(state='disabled')
            
            
        create_label_offset()        
        wcs_header = WCS(header_file)
        Dec = 0
        RA = 0
        try:
            Dec = header_file['CRVAL2']
            RA = header_file['CRVAL1']
        except KeyError as e:
            box_entry.configure(state='normal')
            text_1 = (u"Warning: No WCS founder in header. \n")
            box_entry.insert(INSERT, text_1)
            box_entry.see(END)
            box_entry.configure(state='disabled')
        if RA <= 0: 
            wcs = 0 
            varla7.set("     Spaxel        ID           ")
        else: 
            wcs = 1
            varla7.set("      Spaxel             ID                        RA                                    DEC                      RA-deg                DEC-deg")
        try:
            varla9.set("Object: %s"%(header_file['OBJECT']))
        except Exception as e:  
            try:
                varla9.set("Object: %s"%(header_file['OBJNAME']))
            except Exception as e: 
                varla9.set("Object:  %s "%(name))
            
        arrlambda = np.zeros(pixels) 
        if cdelt == 0:
            cdelt = 1
            text_1 = (u"  CDELT entry not found, using CDELT=1 \n" )
            box_entry.configure(state='normal')
            box_entry.insert(INSERT, text_1)
            box_entry.see(END)
            box_entry.configure(state='disabled')
        for x in range(pixels):
            arrlambda[x] = (crval + x*cdelt)-(crpix-1)*cdelt
        promedio = np.mean(arrlambda)
        if promedio < 100:
            arrlambda = arrlambda*10000
        
   #     if flag_file==1:
    #        ax1.cla()
            
        Integrated_spectrum = np.zeros(pixels)       
   #     flag_file=1
        f2.clf()
        ax1 = f2.add_subplot(projection=wcs_header, slices=('x', 'y', 2))
        dband=np.zeros(tam_x*tam_y)
        imagen_final = np.zeros((tam_y,tam_x))
        array_data=np.reshape(data,(pixels,tam_x*tam_y))
        min_value_la=np.amin(arrlambda)
        max_value_la=np.amax(arrlambda)
        bar_ = Scale(frame3, orient = HORIZONTAL, showvalue= 0, from_=min_value_la+1, to=max_value_la-1, sliderlength = 30, length=480, command = set_bar)
        bar_.grid(row=2,column=2,sticky="nsew",pady=3,padx=5)
        flag_wave = 0
        flag_flux = 0
        flag_band = 0
        flag_explorer=0
        flag_integrate_region=0
        flag_integrate_region2=0
        flag_create_fits = 0
        radius_.set(0)
        varla10.set("Explorer OFF")
        sp1.set(0)
        red_marks = []
        band = 0
        combo1.current(0)
        combo2.current(0)
        update_graph()     
        var.set(1)
        var3.set("")   
        varla5.set(min_value_la)
        varla6.set(max_value_la)
        aux_56 = "%d - %d "%(varla5.get(),varla6.get())
        varlawave.set(aux_56)
        varlaflux.set(" ")
        name_f = "Halpha-KPN0 6547-80A.txt"
        filters_(name_f)
        ax0.set_xlim(xmin=min_value_la,xmax=max_value_la)
        coordinates_(int(tam_x/2),int(tam_y/2))
        label29.config(state=NORMAL)
        label14.config(state=NORMAL)
        label10.config(state=NORMAL)
        label19.config(state=NORMAL)
        
    except Exception as e:    
           print(e)
        #   MessageBox.showerror("Error!","The selected object is not a cube FITS") 
           MessageBox.showerror("Error!",e) 
           button_quit()
          
           

def update_graph():
    global ax0,toolbar,canvas,f,spectrum
    ax0.cla()
    ax0.set_xlabel( 'Wavelength' )
    ax0.set_ylabel( 'Flux' )
    spectrum=data[:,3,3]
    spectrum=np.nan_to_num(spectrum)   
    for i in range(0,pixels):
        if (abs(spectrum[i])>1e30):
            spectrum[i]=0
    ax0.plot(arrlambda,spectrum,color='blue')
    canvas.draw()
    
def set_flux_range():
    global varla3,varla4,flag_flux,red_marks
    if flag_file==1:
        try:
            aux_34 = varlaflux.get()
            aux_34 = aux_34.split('-')
       #     print(aux_34)
            if len(aux_34) == 4:
                aux_34[0] = float (aux_34[1]) *-1
                aux_34[1] = float (aux_34[3]) *-1
            else:
                if len(aux_34) == 3:
                    if aux_34[0]=='' or aux_34[0]==' ':
                        aux_34[0]= float (aux_34[1]) *-1
                        aux_34[1] = aux_34[2]
      #      print(aux_34)           
            varla3.set(aux_34[0])
            varla4.set(aux_34[1])
            if varla3.get() < varla4.get():
                ax0.set_ylim(ymin=varla3.get(),ymax=varla4.get())
             #   canvas.draw()
                flag_flux=1
                for i in red_marks:
                    ax0.axvline(int(i),varla3.get(),varla4.get(),color='red')
                canvas.draw()
                    
            else:
                MessageBox.showerror("Error!","The minimum value should be minimun that the maximus value")
        except Exception as e:   
               print(e)
               MessageBox.showerror("Error!","Please, enter numbers separater by a -")
               reset_flux_range()
               
               
    else:
        MessageBox.showerror("Error!","Please, first choose a file") 
        
def reset_flux_range():
    global flag_flux,varla3,varla4,spectrum
    if flag_file==1:
        flag_flux=0
        varla3.set(0)
        varla4.set(0)
        varlaflux.set("")
        ax0.set_ylim(ymin=np.amin(spectrum),ymax=np.amax(spectrum)*1.2)
        for i in red_marks:
            ax0.axvline(int(i),min_value_da,max_value_da,color='red')
        canvas.draw()
    else:
        MessageBox.showerror("Error!","Please, first choose a file") 
    

def set_wavelength_range():
    global varla5,varla6,arrlambda,flag_wave
    if flag_file==1:
        try:
            aux_56 = varlawave.get()
            aux_56 = aux_56.split('-')
            
            if int(aux_56[0]) >= min_value_la and int(aux_56[1]) <= max_value_la:
                varla5.set(aux_56[0])
                varla6.set(aux_56[1])
                ax0.set_xlim(xmin=varla5.get(),xmax=varla6.get())
                flag_wave = 1
                
                canvas.draw()
            else:
                if int(aux_56[0]) >= min_value_la and int(aux_56[1]) > max_value_la:
                    MessageBox.showwarning("Warning!","The maximun value is %d"%(np.amax(arrlambda)))
                    varla5.set(aux_56[0])
                    varla6.set(aux_56[1])
                    ax0.set_xlim(xmin=varla5.get(),xmax=max_value_la)
                    varla6.set((ax0.get_xlim())[1])
                    aux_56 = "%d - %d "%(varla5.get(),varla6.get())
                    varlawave.set(aux_56)
                    flag_wave = 1
                
                    canvas.draw()
                else:
                    if int(aux_56[0]) < min_value_la and int(aux_56[1]) <= max_value_la:
                        MessageBox.showwarning("Warning!","The minimun value is %d"%(min_value_la))
                        varla5.set(aux_56[0])
                        varla6.set(aux_56[1])
                        ax0.set_xlim(xmin=min_value_la,xmax=varla6.get())
                        varla5.set((ax0.get_xlim())[0])
                        aux_56 = "%d - %d "%(varla5.get(),varla6.get())
                        varlawave.set(aux_56)
                        flag_wave = 1
                        canvas.draw()
                    else:
                        MessageBox.showwarning("Warning!","The minimun value is %d and the maximun value is %d"%(min_value_la,max_value_la))
                        varla5.set((ax0.get_xlim())[0])
                        varla6.set((ax0.get_xlim())[1])
                        aux_56 = "%d - %d "%(varla5.get(),varla6.get())
                        varlawave.set(aux_56)
        except Exception as e:    
               MessageBox.showerror("Error!","Please, enter numbers")
               varla5.set(min_value_la)
               varla6.set(max_value_la)
               aux_56 = "%d - %d "%(varla5.get(),varla6.get())
               varlawave.set(aux_56)
        
    else:
        MessageBox.showerror("Error!","Please, first choose a file") 
    
  
def reset_wavelength_range():
    global varla5,varla6,arrlambda,flag_wave,canvas, ax0
    if flag_file==1:
        if flag_wave == 1:
            flag_wave = 0
            varla5.set(min_value_la) 
            varla6.set(max_value_la)
            aux_56 = "%d - %d "%(varla5.get(),varla6.get())
            varlawave.set(aux_56)
            ax0.set_xlim(xmin=min_value_la,xmax=max_value_la)
            canvas.draw()
    else:
        MessageBox.showerror("Error!","Please, first choose a file") 
    
    
def integrated_region():
    global flag_integrate_region,flag_explorer,label29
    if flag_file == 1:
        if flag_integrate_region == 1 or flag_integrate_region2 ==1:
            MessageBox.showerror("Error!","Please, first reset integrate region") 
        else:
            flag_integrate_region = 1
            label29.config(state=DISABLED)
            if flag_explorer == 0:
                flag_explorer = 1
            varla10.set("Select a pixel")
    else:
        MessageBox.showerror("Error!","Please, first choose a file") 
    
   
def set_offsets():
    
    global band_sticks,f2,ax1
    if flag_file==1:
        if band_sticks.get() == 1:
         #   band_sticks = 1
            f2.clf()
            ax1 = f2.add_axes( (0.05, .15, .90, .80), frameon=False)
            filters_(name_f)
        else:
            f2.clf()
            ax1 = f2.add_subplot(projection=wcs_header, slices=('x', 'y', 2))
            filters_(name_f)
    else:
        MessageBox.showerror("Error!","Please, first choose a file")
        band_sticks.set(0)
    
  

def reset_integrated_region():
    global flag_integrate_region2, flag_integrate_region,integrated_x,integrated_y,flag_create_fits
    if flag_file == 1:
        if flag_integrate_region2 == 1 and flag_integrate_region == 1:
            flag_integrate_region2 = 0
            flag_integrate_region = 0
            text_1 = (u"    \n Circle erased    \n")
            box_entry.configure(state='normal')
            box_entry.insert(INSERT, text_1)
            box_entry.see(END)
            box_entry.configure(state='disabled')   
            varla10.set("Explorer OFF")
            integrated_x = []
            integrated_y = []
            filters_(name_f)
            radius_.set(0)
            label29.config(state=NORMAL)
            flag_create_fits = 0
    else:
        MessageBox.showerror("Error!","Please, first choose a file") 
       
def set_bar(bar_1):
    varla11.set(bar_1)
    

def mark_wavelength():
    
    global band,ax0,arrlambda,canvas,spectrum,red_marks
    if flag_file==1:
        mark = var3.get()
        red_marks = mark.split(',')
        b = 0
        try:
            for i in red_marks:
                if int(i) < min_value_la or int(i) > max_value_la:
                    b = 1
                else:
                    ax0.axvline(int(i),min_value_da,max_value_da,color='red')
            if b == 1:
                MessageBox.showwarning("Warning!","The minimun value is %d and the maximun value is %d"%(min_value_la,max_value_la))
                red_marks = []
                var3.set("")
            else:            
              canvas.draw()
        except Exception as e: 
               MessageBox.showerror("Error!","Please, enter numbers") 
               red_marks = []
               var3.set("")
    else:
        MessageBox.showerror("Error!","Please, first choose a file") 
      
def set_mark_wavelength():
    global red_marks,var3
    if flag_file==1:
        red_marks = []
        var3.set("")
        ax0.cla()
        ax0.set_xlabel( 'Wavelength' )
        ax0.set_ylabel( 'Flux' )
        ax0.plot(arrlambda,spectrum,color='blue')
        if flag_flux==1:
            ax0.set_ylim(ymin=varla3.get(),ymax=varla4.get())
            ax0.plot(arrlambda,res*varla4.get(),'--',color = 'orange')
        else:
            ax0.plot(arrlambda,res*np.amax(spectrum)*1.2,'--',color = 'orange')
            ax0.set_ylim(ymin=np.amin(spectrum),ymax=np.amax(spectrum)*1.2)
        if flag_wave == 1:
            ax0.set_xlim(xmin=varla5.get(),xmax=varla6.get())
        else:
            ax0.set_xlim(xmin=min_value_la,xmax=max_value_la)
        canvas.draw()
    else:
        MessageBox.showerror("Error!","Please, first choose a file") 
    
    
           
def set_band():
    
    global varla11,band,name_f,flag_band
    if flag_file==1:
        try:
            band = varla11.get()
            flag_band=1
            filters_(name_f)
            set_scaling() 
        except Exception as e: 
               MessageBox.showerror("Error!","Please, enter numbers") 
    else:
        MessageBox.showerror("Error!","Please, first choose a file")
    

def set_scaling():
    global saved_image,canvas2,ax1,f2,canvas2
    if flag_file==1:
        scaling = var.get()
        if sp1.get()==1:
            cmap_1=maps_array_inv[combo1.current()]
        else:
            cmap_1=maps_array[combo1.current()]
        if band_sticks.get() == 0:
            f2.clf()
            ax1 = f2.add_subplot(projection=wcs_header, slices=('x', 'y', 2))
        else:
            ax1.cla()
            ax1.set_xlabel( 'RA (arcsec)' )
            ax1.set_ylabel( 'DEC (arcsec)' )
            ax1.set_xticks(x_ticks)
            ax1.set_yticks(y_ticks)
            ax1.set_xticklabels(x_ticks_l)
            ax1.set_yticklabels(y_ticks_l)      
        if scaling == 1:
            saved_image=ax1.imshow(imagen_final,cmap=cmap_1,interpolation='nearest',origin='lower' )
            lineal = simple_norm(imagen_final, stretch='linear')
            saved_image.set_norm(lineal)
        else:
            if scaling == 2:
                total = sigma_clip(imagen_final, sigma=2)
                saved_image=ax1.imshow(total,cmap=cmap_1,interpolation='nearest',origin='lower' )
            else:
                if scaling == 3:
                    saved_image=ax1.imshow(imagen_final,cmap=cmap_1,interpolation='nearest',origin='lower' )
                    asin_h = simple_norm(imagen_final, stretch='asinh')
                    saved_image.set_norm(asin_h)
                else:
                    if scaling == 4:
                        saved_image=ax1.imshow(imagen_final,cmap=cmap_1,interpolation='nearest',origin='lower' )
                        power = 2.0
                        power_l = simple_norm(imagen_final, stretch='power', power=power)
                        saved_image.set_norm(power_l)
                    else:
                        if scaling == 5:
                            saved_image=ax1.imshow(imagen_final,cmap=cmap_1,interpolation='nearest',origin='lower')
                            raiz_c = simple_norm(imagen_final, stretch='sqrt')
                            saved_image.set_norm(raiz_c)
                        else:
                            if scaling == 6:
                                saved_image=ax1.imshow(imagen_final,cmap=cmap_1,interpolation='nearest',origin='lower' )
                                img_cdf, bin_centers = exposure.cumulative_distribution(imagen_final)
                                final = np.interp(imagen_final,bin_centers,img_cdf)
                        #        ax1.cla()
                                saved_image=ax1.imshow(final,cmap=cmap_1,interpolation='nearest',origin='lower' )
                            else:
                                if scaling == 7:
                                    saved_image=ax1.imshow(imagen_final,cmap=cmap_1,interpolation='nearest',origin='lower' )
                                    sigma=1
                                    norm_img = np.zeros((36,36))
                                    imagen_pi = cv2.normalize(imagen_final,norm_img, -math.pi, math.pi, cv2.NORM_MINMAX)
                                    una=1/(sigma*math.sqrt(2*math.pi))
                                    cuadrado=np.power(imagen_pi,2)
                                    division= cuadrado/(2*(sigma**2))
                                    dos = np.exp(-division)
                                    total=una*dos
                         #           ax1.cla()
                                    saved_image=ax1.imshow(total,cmap=cmap_1,interpolation='nearest',origin='lower' )
                                else:
                                    saved_image=ax1.imshow(imagen_final,cmap=cmap_1,interpolation='nearest',origin='lower' )
                                    log_a = 100
                                    loga = simple_norm(imagen_final, stretch='log', log_a=log_a)
                                    saved_image.set_norm(loga)
                                    
        if flag_integrate_region == 1 and flag_integrate_region2 == 1:
                cir = patches.Circle((cir_x,cir_y),int(radius_.get()),edgecolor='red',fill = False)
                ax1.add_patch(cir)
                
        canvas2.draw()
    else:
        MessageBox.showerror("Error!","Please, first choose a file") 
        var.set(1)

    
def set_filter(event):
   
    global name_f,flag_band,combo2
    if flag_file==1:
        combo_2 = event.widget.get()
        name_f = combo_2 + '.txt'
        flag_band=0
        filters_(name_f)
    else:
        MessageBox.showerror("Error!","Please, first choose a file")
        combo2.current(0)
    
def filters_(name_1):
    global name_f,ax0,canvas,arrlambda,spectrum,band,arr_w,arr_f,flag_flux,varla3,varla4,varla5,varla6,imagen_final,dband,res,bar_
    lineas = []
    cwl = 0
    filter_dir=os.getcwd()
    if flag_system == 1:
        filter_dir = filter_dir+'/Filters/'
    else:
        filter_dir = filter_dir+'\\Filters\\'
    nuev_=filter_dir+name_1
    if os.path.isfile(nuev_):
        archivo = open(nuev_, 'r')
        for linea in archivo.readlines():
            if linea.startswith("# CWL:")==True:
                cwl = float(linea[7:len(linea)-1])
            if linea.startswith("#")==False and linea.isalnum()==False and linea!= '':
                linea = linea.rstrip('\r\n')
                lineas.append(linea)
        lineas = list(filter(bool,lineas))
        arr_w = np.zeros(len(lineas))
        arr_f =  np.zeros(len(lineas))
        j= 0
        for i in lineas:
            primeralinea=i.split()
            arr_w[j] = primeralinea[0]
            arr_f[j] = primeralinea[1]
            j = j+1
        if (arr_w[0] < min_value_la and arr_w[j-1] < min_value_la) or (max_value_la < arr_w[0] and max_value_la < arr_w[j-1]):
            text_1 = (u" WARNING: Using a V-band filter within the data wavelength \n limits: %d - %d \n "%(min_value_la,max_value_la))
            box_entry.configure(state='normal')
            box_entry.insert(INSERT, text_1)
            box_entry.see(END)
            box_entry.configure(state='disabled')
            archivo.close()
            lineas = []
            name_f = "V-Johnson.txt"
            nuev_=filter_dir+"V-Johnson.txt"
            combo2.current(3)
            archivo = open(nuev_, 'r')
            for linea in archivo.readlines():
                if linea.startswith("# CWL:")==True:
                    cwl = float(linea[7:len(linea)-1])
                if linea.startswith("#")==False and linea.isalnum()==False and linea!= '':
                    linea = linea.rstrip('\r\n')
                    lineas.append(linea)
            
            lineas = list(filter(bool,lineas))
            arr_w = np.zeros(len(lineas))
            arr_f =  np.zeros(len(lineas))
            j= 0
            for i in lineas:
                primeralinea=i.split()
                arr_w[j] = primeralinea[0]
                arr_f[j] = primeralinea[1]
                j = j+1
            ##-hasta aqui
            if flag_band == 0:
                band= np.mean(arrlambda)  
                varla11.set(band)
                bar_.set(band)
           
        
        else:
           if flag_band==0:
               band=cwl
               varla11.set(band)
               bar_.set(band)
               
        shift = cwl/band
        new1 = arr_w/shift
        res=interpolate.InterpolatedUnivariateSpline(new1, arr_f)(arrlambda)
        posiciones = []
        posiciones2 = []
        for i in range(len(arrlambda)):
            if arrlambda[i] <= np.amin(new1) or arrlambda[i] >= np.amax(new1):
                posiciones2.append(i)
            if arrlambda[i] >= np.amin(new1) and arrlambda[i] <= np.amax(new1):
                posiciones.append(i)
        res[posiciones2]=0
        convolve = res[posiciones]
        for i in range(len(dband)):
            dband[i] = sum(convolve*array_data[posiciones,i])
        imagen_final=np.reshape(dband,(tam_y,tam_x))
        ax0.clear()
        if flag_integrate_region == 1 and flag_integrate_region2 == 1:
            if flag_flux==1:
                ax0.set_ylim(ymin=varla3.get(),ymax=varla4.get())
                ax0.plot(arrlambda,res*varla4.get(),'--',color = 'orange')
                for i in red_marks:
                        ax0.axvline(int(i),varla3.get(),varla4.get(),color='red')
            else:
             #   ax0.plot(arrlambda,res*np.amax(Integrated_spectrum)*1.2,'--',color = 'orange')
                ax0.set_ylim(ymin=np.amin(Integrated_spectrum),ymax=np.amax(Integrated_spectrum)*1.2)
                for i in red_marks:
                        ax0.axvline(int(i),min_value_da,max_value_da,color='red')
            
            ax0.plot(arrlambda,Integrated_spectrum,color='red')
            if flag_wave == 1:
                    ax0.set_xlim(xmin=varla5.get(),xmax=varla6.get())
            else:
                    ax0.set_xlim(xmin=min_value_la,xmax=max_value_la)
           
        else:
            if flag_flux==1:
                ax0.set_ylim(ymin=varla3.get(),ymax=varla4.get())
                ax0.plot(arrlambda,res*varla4.get(),'--',color = 'orange')
                for i in red_marks:
                        ax0.axvline(int(i),varla3.get(),varla4.get(),color='red')
            else:
                ax0.plot(arrlambda,res*np.amax(spectrum)*1.2,'--',color = 'orange')
                ax0.set_ylim(ymin=np.amin(spectrum),ymax=np.amax(spectrum)*1.2)
                for i in red_marks:
                        ax0.axvline(int(i),min_value_da,max_value_da,color='red')
            
            ax0.plot(arrlambda,spectrum,color='blue')
            if flag_wave == 1:
                    ax0.set_xlim(xmin=varla5.get(),xmax=varla6.get())
            else:
                    ax0.set_xlim(xmin=min_value_la,xmax=max_value_la)
                    
        nomb_a = name_f.split(".")
        text_1 = (u" Filter used %s with band %d \n"%(nomb_a[0],band))
        box_entry.configure(state='normal')
        box_entry.insert(INSERT, text_1)
        box_entry.see(END)
        box_entry.configure(state='disabled')    
        canvas.draw()
        archivo.close()
        set_scaling()

        
    else:
        MessageBox.showerror("Error!","Not exists the folder of filters") 
    
    directorio= os.listdir("./")
    
    
def set_color_map(event=''):
  
    global combo1,saved_image,canvas2
    if flag_file==1:
        if sp1.get()==1:
            saved_image.set_cmap(maps_array_inv[combo1.current()])
        else:
            saved_image.set_cmap(maps_array[combo1.current()])
        canvas2.draw()
    else:
        sp1.set(0)
        combo1.current(0)
        MessageBox.showerror("Error!","Please, first choose a file") 
    
def button_quit_destroy():
    window_.destroy()
    
def button_quit():
    global flag_file,ax0,f2,canvas,canvas2,hi_data
    if flag_file == 1:
        hi_data.close()
        flag_file=0
        ax0.cla()
        f2.clf()
        canvas.draw()
        canvas2.draw()
        var.set(1)
        varla8.set("")
        varla10.set("")
        varla9.set("")
        varla7.set("")
        varla5.set(0)
        varla6.set(0)
        varlaflux.set("")
        varlawave.set("")
        radius_.set(0)
        combo2.current(0)
        combo1.current(0)
        var3.set("")
        varla11.set(0)
        varla1.set("")
        box_entry.configure(state='normal')
        box_entry.delete("1.0","end")
        box_entry.configure(state='disabled')
        label29.config(state=DISABLED)
        label14.config(state=DISABLED)
        label10.config(state=DISABLED)
        label19.config(state=DISABLED)
        
        
    else:
        MessageBox.showerror("Error!","Please, first choose a file") 
    
    
def get_resolution():
  #  user32 = ctypes.windll.user32
  #  user32.SetProcessDPIAware()
  #  width, height = user32.GetSystemMetrics(0), user32.GetSystemMetrics(1)
    return 1440-50,760-50
  #  width = 1410 #1350x780
  #  height = 780
   # return width-50,height-50
    
   # user32 = ctypes.windll.user32
   # user32.SetProcessDPIAware()
   # width, height = user32.GetSystemMetrics(0), user32.GetSystemMetrics(1)
  #  return width-50,height-50
    
def operative_system():
    global flag_system
    sistema = platform.system()
    if sistema == "Windows":
        flag_system = 0
    else:
        flag_system = 1
        
if __name__=="__main__":    
   
    file_dir = ''
    band = 0
    arr_w  = np.zeros(2)
    arr_f = np.zeros(2)
    name = ''
    hi_data = ''
    data = ''
    header_file = ''
    wcs_header = ''
    crval  = 0
    cdelt  = 0
    crpix  = 0
    tam_x  = 0
    tam_y  = 0
    pixels = 1
    x_ticks = []
    y_ticks = []
    x_ticks_l = []
    y_ticks_l = []
    Integrated_spectrum =[]
    integrated_x = []
    integrated_y = []
    
    Dec = 0
    RA = 0
    wcs = 0 
    arrlambda = np.zeros(pixels) 
    array_data = 0
    dband=0
    res = []
    cir_x = 0
    cir_y = 0
    name_f = "" 
    
    
    sout = 0
    infinite= 0
    spectrum= 0    
    min_value_da= 0
    max_value_da= 0
    
    
    #--------banderas
    
    flag_explorer = 0
    flag_flux=0
    flag_wave=0
    flag_band=0
    flag_file=0
  #  band_sticks = 0
    flag_integrate_region = 0
    flag_integrate_region2 = 0
    flag_create_fits = 0
    flag_system = 0 # 0 windows 1 ubuntu/mac
    operative_system()
    
    #parte auxiliar
    
    red_marks = []
    maps_array = []
    maps_array_inv = []
    ax1 = 0
    saved_image = 0
    
    #parte para crear los mapas solo una vez:
    prism = matplotlib.colors.LinearSegmentedColormap.from_list('custom prism', [(0,    "white"),(0.2, '#000000'),(0.4, '#8b0000'),(0.6, '#f63f2b'),(0.8, '#15E818'),(1, '#1139d9'  )], N=256)
    stern = matplotlib.colors.LinearSegmentedColormap.from_list('custom stern',[(0,    "white"),(0.2, '#8b0000'),(0.3, '#e42121'),(0.4, '#252850'),(0.6, '#0588EF'), (0.8, '#3b83bd'),(1, '#c6ce00'  )], N=256)
    std = matplotlib.colors.LinearSegmentedColormap.from_list('custom Std-Gamma', [(0,    "white"),(0.2, '#0000ff'),(0.4, '#2178E4'),(0.6, '#ff0000'),(0.8, '#ff8000'),(1, '#ffff00'  )], N=256)
    BGRY = matplotlib.colors.LinearSegmentedColormap.from_list('custom BGRY', [(0,    "white"),(0.2, '#ff8000'),(0.4, '#EFEE05'),(0.6, '#EF5A05'),(0.8, '#51EF05'),(1, '#0000ff'  )], N=256)
    califa = matplotlib.colors.LinearSegmentedColormap.from_list('custom CALIFA special', [(0,    "white"),(0.25, '#00008B'),(0.5, '#B2FFFF'),(0.62, '#B2FFFF'),(0.75, '#ff4000'),(1, '#008f39'  )], N=256)
    ping = matplotlib.colors.LinearSegmentedColormap.from_list('custom Pingsoft-special', [(0,    "white"),(0.25, '#00008B'),(0.5, '#3b83bd'),(0.75, '#ff8000'),(1, '#ffff00'  )], N=256)
    
    
    prism_r= matplotlib.colors.LinearSegmentedColormap.from_list('custom prism inv', [(0,    '#1139d9'),(0.2, '#15E818'),(0.4, '#f63f2b'),(0.6, '#8b0000'),(0.8, '#000000'),(1, "white"  )], N=256)
    stern_r= matplotlib.colors.LinearSegmentedColormap.from_list('custom strn inv', [(0,    '#c6ce00'),(0.2, '#3b83bd'),(0.3, '#0588EF'),(0.4, '#252850'),(0.6, '#e42121'),(0.8, '#8b0000'),(1, "white"  )], N=256)
    std_r= matplotlib.colors.LinearSegmentedColormap.from_list('custom std inv', [(0,    '#ffff00'),(0.2, '#ff8000'),(0.4, '#ff0000'),(0.6, '#2178E4'),(0.8, '#0000ff'),(1, "white"  )], N=256)
    BGRY_r= matplotlib.colors.LinearSegmentedColormap.from_list('custom BGRY inv', [(0,    '#0000ff'),(0.2, '#51EF05'),(0.4, '#EF5A05'),(0.6, '#EFEE05'),(0.8, '#ff8000'),(1, "white"  )], N=256)
    califa_r= matplotlib.colors.LinearSegmentedColormap.from_list('custom CALIFA inv', [(0,    '#008f39'),(0.25, '#ff4000'),(0.5, '#B2FFFF'),(0.62, '#B2FFFF'),(0.75, '#00008B'),(1, "white"  )], N=256)
    ping_r= matplotlib.colors.LinearSegmentedColormap.from_list('custom Pingsoft inv', [(0,    '#ffff00'),(0.25, '#ff8000'),(0.5, '#3b83bd'),(0.75, '#00008B'), (1, "white"  )], N=256)
    
    maps_array.append('Blues')
    maps_array.append('Reds')
    maps_array.append('Greens')
    maps_array.append('Greys')
    maps_array.append(ping)
    maps_array.append(califa)
    maps_array.append('rainbow')
    maps_array.append(BGRY)
    maps_array.append(prism)
    maps_array.append(stern)
    maps_array.append(std)
    maps_array_inv.append('Blues_r')
    maps_array_inv.append('Reds_r')
    maps_array_inv.append('Greens_r')
    maps_array_inv.append('Greys_r')
    maps_array_inv.append(ping_r)
    maps_array_inv.append(califa_r)
    maps_array_inv.append('rainbow_r')
    maps_array_inv.append(BGRY_r)
    maps_array_inv.append(prism_r)
    maps_array_inv.append(stern_r)
    maps_array_inv.append(std_r)
    imagen_final = 0
    cmap = maps_array[0]
    
    
    #####################------------------------
    window_=Tk()
    window_.title("IFS Explorer")
    if flag_system == 0:
        fontStyle = tkFont.Font(family="Arial", size=8,slant="roman")
    else:
        fontStyle = tkFont.Font(family="Arial", size=12,slant="roman")
    
    res_width,res_height=get_resolution()
    valor_res="%dx%d"%(res_width,res_height)
    window_.geometry(valor_res)
    
    
    
    
    
    ###------------variables ()
    radius_ = IntVar()
    band_sticks = IntVar()
    min_value_la = 0
    max_value_la = 0
    varla1 = StringVar()
    varla2 = StringVar()
    varlaflux = StringVar() 
    varla3 = DoubleVar()
    varla4 = DoubleVar()
    varlawave = StringVar()
    varla5 = IntVar()
    varla6 = IntVar()
    varla7 = StringVar()
    varla8 = StringVar()
    varla9 = StringVar()
    varla10 = StringVar()
    varla11 = DoubleVar()
    varla12 = IntVar()
    var = IntVar() 
    var.set(1)
    var3 = StringVar()
    varsc = DoubleVar()
    var3 = StringVar()
    sp = IntVar()
    sp1 = IntVar()
    varsc = DoubleVar()
    
    
    
    
    #image_ de logotipo
    image_dir=os.getcwd()
    if flag_system == 1:
        image_dir = image_dir+'/Img/logoIFSexplorer.png'
    else:
        image_dir = image_dir+'\\Img\\logoIFSexplorer.png'
    image_=Image.open(image_dir)
    image_=image_.resize((250,250),Image.ANTIALIAS)
    pImage = ImageTk.PhotoImage(image_)
    window_.configure(background='light gray')

    background_=Label(window_,image=pImage).grid(row = 0, column = 0, pady = 2) 
    
    #frame para navegador
    frame = Frame()
    frame.config(bg="light grey",bd=1, relief="groove",width=605,height=275)
    frame.grid(row = 0, column = 1, pady = 2,padx=2)
    frame.grid_propagate(False)
    
    
    Label(frame, text="FITS",width="10", height="2",font=fontStyle).grid(row=0,column=0,sticky = W,pady=2,padx = 2)  
    Entry(frame,state=DISABLED,textvariable=varla1,font=fontStyle).grid(row=0,column=1, columnspan=5,sticky="ew",padx = 2) 
    Button(frame, font=fontStyle,text="Browse", command=new_file,width="6").grid(row=0,column=6,pady = 2,padx = 2)
    Button(frame, text="Quit",command =button_quit_destroy,font=fontStyle,width="6").grid(row=0,column=7,pady = 2)
    Entry(frame,state=DISABLED,font=fontStyle, textvariable=varla9, width="53").grid(row=1,column=0,columnspan=8)
    Entry(frame,state=DISABLED, font=fontStyle,textvariable=varla10, width="53").grid(row=2,column=0,columnspan=8)
    
    box_entry = scrolledtext.ScrolledText(frame)
    box_entry.config(width = 50, height=10,font=fontStyle,yscrollcommand=TRUE)
    box_entry.grid(row = 3, column = 0, columnspan=8)  
    box_entry.configure(state='disabled')
    Label(frame, text="Flux range",width="12", height="2",font=fontStyle).grid(row=4,column=0,pady = 5)
    label14 = Entry(frame,textvariable=varlaflux, width="10",font=fontStyle)
    label14.grid(row=4,column=1,pady = 5,padx=10)
    label14.config(state=DISABLED)
    varlaflux.set(" ")
    Button(frame, text="Set", command=set_flux_range,font=fontStyle,width="6").grid(row=4,column=2,padx = 2,pady = 5)
    Button(frame, text="Reset", command=reset_flux_range,font=fontStyle,width="6").grid(row=4,column=3,padx=10,pady = 5)
    
    
    Label(frame, text="Wavelenth range",width="15", height="2",font=fontStyle).grid(row=4,column=4,pady = 5)
    label19 = Entry(frame,textvariable=varlawave, width="10",font=fontStyle)
    label19.grid(row=4,column=5,padx=10,pady = 5)
    varlawave.set(" ")
    label19.config(state=DISABLED)
    Button(frame, text="Set", command=set_wavelength_range,font=fontStyle,width="6").grid(row=4,column=6,pady = 5)
    Button(frame, text="Reset", command=reset_wavelength_range,font=fontStyle,width="6").grid(row=4,column=7,pady=5,padx=10)
    
    #Frame de Integrated region
    frame9 = Frame()
    frame9.config(bg="light grey", bd=1, relief="groove",padx=10,width=110, height=275)
    frame9.grid(row = 0, column = 2, pady = 2,sticky="nsew") 
    frame9.grid_propagate(False)
    
    Button(frame9, text="Integrated region",font=fontStyle,command=integrated_region).grid(columnspan=2,row=0,column=0,pady=3)
    Label(frame9, text="Radius",font=fontStyle).grid(columnspan=2,row=1,column=0,pady=2,padx=2,sticky="nsew")
    label29=Entry(frame9, width="8",font=fontStyle,textvariable=radius_)
    label29.grid(row=2,column=0,padx=3)    
  #  label29.config(justify=RIGHT) #Dirección del texto
    label29.config(state=DISABLED)
    
    
    Button(frame9, text="Reset",command=reset_integrated_region,font=fontStyle,width="6").grid(row=2,column=1)
    Button(frame9, text="Create Files",command=create_files,font=fontStyle,width="10").grid(row=3,column=0,columnspan=2)
    
    
    
    Label(frame9, text="Display axis",font=fontStyle).grid(columnspan=2,row=4,column=0,pady=6,sticky="nsew") 
    Radiobutton(frame9, text="RA-DEC", variable=band_sticks, value=0, command=set_offsets,font=fontStyle).grid(columnspan=2,row=5,column=0,sticky="nsew")
    Radiobutton(frame9, text="Offset", variable=band_sticks, value=1, command=set_offsets,font=fontStyle).grid(columnspan=2,row=6,column=0,sticky="nsew")
    
    Label(frame9, text="Mark wavelength",font=fontStyle).grid(columnspan=2,row=7,column=0,sticky="nsew",pady=6)
    label10 =Entry(frame9,textvariable=var3, width="10",font=fontStyle)
    label10.grid(columnspan=2,row=8,column=0,pady=3)
    label10.config(state=DISABLED)
    Button(frame9, text="Set", command=mark_wavelength,font=fontStyle,width="6").grid(row=9,column=0,pady=3)
    Button(frame9, text="Reset", command=set_mark_wavelength,font=fontStyle,width="6").grid(row=9,column=1,pady=3)
    
    #Frame para display
    frame10 = Frame()
 #  frame10.config(bg="light grey", bd=1, relief="groove",padx=10,width=395, height=200)
   #frame10.config(bg="light grey", bd=1, relief="groove",padx=1,width=380,height=200)
    frame10.config(bg="light grey", bd=1, relief="groove",padx=1,width=372,height=200)
    frame10.grid(row = 0, column = 3, pady = 2) 
    frame10.grid_propagate(False)
  #  frame10.config(width=480,height=320)
    
    Label(frame10, text="Display",font=fontStyle).grid(columnspan=2,row=0,column=0,sticky="nsew",pady=10,padx=5)
    Label(frame10, text="Color Map",font=fontStyle).grid(row=1,column=0,columnspan=2,sticky=W)
    combo1 = ttk.Combobox(frame10, width = 25,state="readonly",font=fontStyle) 
    combo1['values'] = ( 'Blue scaling', 
                                        'Red scaling',
                                        'Green scaling',
                                        'Grayscale',
                                        'PINGSoft special',
                                        'CALIFA-special',
                                        'Rainbow',
                                        'BGRY',
                                        'Prism',
                                        'Stern',
                                        'Std-Gamma')   
    combo1.current(0)
    combo1.grid(column =0, row = 3,sticky = W,columnspan=2) 
    combo1.bind("<<ComboboxSelected>>", set_color_map)  
    Label(frame10, text="Filter",font=fontStyle).grid(row=4,column=0,sticky=W)
    combo2 = ttk.Combobox(frame10, width = 25,state="readonly",font=fontStyle) 
    combo2['values'] = ( 'Halpha-KPN0 6547-80A', 
                                        'HALPHA-CTI0 6586-20A',
                                        'B-Johnson (1965)',
                                        'V-Johnson',
                                        'u-SDSS-III',
                                        'g-SDSS-III',
                                        'r-SDSS-III',
                                        'i-SDSS-III',
                                        'B-Bessell (1990)',
                                        'V-Bessell',
                                        'R-Bessell',
                                        'B-KPN0-Harris',
                                        'V-KPN0-Harris',
                                        'R-KPN0-Harris')   
    combo2.current(0)
    combo2.grid(column = 0, row = 5,sticky = W,columnspan=2) 
    combo2.bind("<<ComboboxSelected>>", set_filter)
    
    Checkbutton(frame10,text="Invert color map", variable=sp1, onvalue=1, offvalue=0, command=set_color_map,font=fontStyle).grid(row=6,column=0,sticky = W,pady=6)
    
    Label(frame10, text="Scaling",font=fontStyle).grid(row=0,column=2,columnspan=2,sticky="nsew",pady=10,padx=5)
    
    #opciones de boton
    Radiobutton(frame10, text="Linear", variable=var, value=1, command=set_scaling,font=fontStyle).grid(sticky=W,row=3,column=2,pady=5,padx=5)
    Radiobutton(frame10, text="2% clipping", variable=var, value=2, command=set_scaling,font=fontStyle).grid(sticky=W,row=4,column=2,pady=3,padx=5)
    Radiobutton(frame10, text="Asinh", variable=var, value=3, command=set_scaling,font=fontStyle).grid(sticky=W,row=5,column=2,pady=3,padx=5)
    Radiobutton(frame10, text="Power-Law", variable=var, value=4, command=set_scaling,font=fontStyle).grid(sticky=W,row=6,column=2,pady=6,padx=5)
    Radiobutton(frame10, text="Square Root", variable=var, value=5, command=set_scaling,font=fontStyle).grid(sticky=W,row=3,column=3,pady=5,padx=5)
    Radiobutton(frame10, text="Hist Equal", variable=var, value=6, command=set_scaling,font=fontStyle).grid(sticky=W,row=4,column=3,pady=3,padx=5)
    Radiobutton(frame10, text="Gaussian", variable=var, value=7, command=set_scaling,font=fontStyle).grid(sticky=W,row=5,column=3,pady=3,padx=5)
    Radiobutton(frame10, text="Logarithmic", variable=var, value=8, command=set_scaling,font=fontStyle).grid(sticky=W,row=6,column=3,pady=6,padx=5)
 #   frame10.grid_propagate(False)
    frame3 = Frame()
    frame3.config(bg="light grey", bd=1, relief="groove", pady=10,padx=10,width=850, height=80)
    frame3.grid(row = 1, column = 0, pady = 2, columnspan=2) 
    frame3.grid_propagate(False)
    
    Entry(frame3,state=DISABLED,width=137, textvariable=varla7,font=fontStyle).grid(row=0,column=0,columnspan=4,sticky="nsew")
    Entry(frame3,state=DISABLED, textvariable=varla8,font=fontStyle).grid(row=1,column=0,columnspan=4,sticky="nsew")
    Label(frame3, text="Shift Filter",font=fontStyle).grid(row=2,column=0,sticky="nsew",pady=3,padx=5)
    Entry(frame3,state=DISABLED, textvariable=varla11, width="10",font=fontStyle).grid(row=2,column=1,pady=5,padx=3,sticky="nsew")
    bar_ = Scale(frame3, orient = HORIZONTAL, showvalue= 0, from_=min_value_la+1, to=max_value_la-1, sliderlength = 30, length=480,command = set_bar )
    bar_.grid(row=2,column=2,sticky="nsew",pady=3,padx=5)
    Button(frame3, text="Apply", command=set_band,font=fontStyle,width="6").grid(row=2,column=3,pady=3,padx=5,sticky="nsew")
  
    #frame para image_
    frame6 = Frame()
  #  frame6.config(bg="light grey", bd=1, relief="groove", pady=1,padx=1,width=521, height=419)
    frame6.config(bg="light grey", bd=1, relief="groove", pady=1,padx=1)
    frame6.grid(row = 1, column = 2, pady = 2, columnspan=2,rowspan=3) 
 #   frame6.grid_propagate(False)
    
    f2 = Figure( figsize=(6.5, 4.8), dpi=80 )
    
    saved_image = 0 
    canvas2 = FigureCanvasTkAgg(f2, master=frame6)
    canvas2.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    canvas2.mpl_connect("motion_notify_event", move_mouse)
    canvas2.mpl_connect("button_press_event", onclick_)
    canvas2.draw()
    toolbar2 = NavigationToolbar2Tk(canvas2,frame6)
    toolbar2.pack()
    toolbar2.update() 
    
    #frame para grafica
    f = Figure( figsize=(10.3, 3.7), dpi=80 )
    ax0 = f.add_axes( (0.088, .15, .90, .80), frameon=False)
    frame5 = Frame() #width=610, height=180
    frame5.config(bg="light grey", bd=1, relief="groove", pady=2,padx=2)
    frame5.grid(row = 2, column = 0, pady = 2, columnspan=2) 
 #   frame5.grid_propagate(False)
    
    canvas = FigureCanvasTkAgg(f, master=frame5)
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    canvas.draw()
    toolbar = NavigationToolbar2Tk(canvas,frame5 )
    toolbar.pack()
    toolbar.update()
   
    
    window_.resizable(0,0)
    window_.mainloop()
