import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import CheckButtons
from matplotlib.figure import Figure

class Main(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        container = tk.Frame(self)
        container.pack(fill="both", expand=True)
        pages = (HomePage, IgnitionACOMP, RockwellACOMP, MWUV, Others)
        self.frames = {}
        for page in pages:
            self.frames[page] = page(container, self)
        self.show_frame(HomePage)
        self.wm_title("R&D Data Analysis V1.02")
    def show_frame(self, page):
        frame = self.frames[page]
        frame.tkraise()
class HomePage(tk.Frame):
    def __init__(self, container, controller):
        tk.Frame.__init__(self, container)
        self.grid(row=0, column=0, sticky='nsew')
        button1 = ttk.Button(self, text="Ignition ACOMP", command=lambda: controller.show_frame(IgnitionACOMP), width=20)
        button1.grid(row=0, column=0, sticky='nsew')
        button2 = ttk.Button(self, text="Rockwell ACOMP", command=lambda: controller.show_frame(RockwellACOMP), width=20)
        button2.grid(row=1, column=0, sticky='nsew')
        button3 = ttk.Button(self, text="Multi-Wavelength UV", command=lambda: controller.show_frame(MWUV), width=20)
        button3.grid(row=2, column=0, sticky='nsew')
        button4 = ttk.Button(self, text="Customize System", command=lambda: controller.show_frame(Others), width=20)
        button4.grid(row=3, column=0, sticky='nsew')
        button5 = tk.Button(self, text="Quit", bg="red", command=quit, width=20)
        button5.grid(row=4, column=0, sticky='nsew')
class IgnitionACOMP(tk.Frame):
    def __init__(self, container, controller):
        tk.Frame.__init__(self, container)
        self.grid(row=0, column=0, stick='nsew')
        self.value = 0
        self.button_frame = ttk.Frame(self)
        self.button_frame.grid(row=0, column=0, sticky='nsew')
        self.canvas_frame = ttk.Frame(self)
        self.canvas_frame.grid(row=0, column=1, sticky='nsew')
        self.loadbutton = ttk.Button(self.button_frame, text='Load File',
                                    command=self.load_file, width=20)
        self.loadbutton.grid(row=0, column=0, stick='nsew')
        self.plotbutton = ttk.Button(self.button_frame, text='Generate Plots',
                                   command=self.data_plot, width=20)
        self.plotbutton.grid(row=1, column=0, sticky='nsew')
        self.back = ttk.Button(self.button_frame, text='Back the Home',
                              command=lambda: controller.show_frame(HomePage), width=20)
        self.back.grid(row=2, column=0, sticky='nsew')
        self.LF1 = ttk.LabelFrame(self.button_frame, text='Zimm Plot Angle')
        self.LF1.grid(row=4, column=0, columnspan=1, sticky='nsew')
        self.bv_0 = tk.BooleanVar()
        self.bv_0.set(True)
        self.bv_1 = tk.BooleanVar()
        self.bv_1.set(True)
        self.bv_2 = tk.BooleanVar()
        self.bv_2.set(True)
        self.bv_3 = tk.BooleanVar()
        self.bv_3.set(True)
        self.bv_4 = tk.BooleanVar()
        self.bv_4.set(True)
        self.cb0 = tk.Checkbutton(self.LF1, text='45deg', variable=self.bv_0,
                             command=self.bv_0.get)
        self.cb0.grid(row=1, column=0, sticky='nsew')
        self.cb0.bind("<ButtonRelease>", self.bv_get)
        self.cb1 = tk.Checkbutton(self.LF1, text='65deg', variable=self.bv_1,
                             command=self.bv_1.get)
        self.cb1.grid(row=2, column=0, sticky='nsew')
        self.cb1.bind("<ButtonRelease>", self.bv_get)
        self.cb2 = tk.Checkbutton(self.LF1, text='90deg', variable=self.bv_2,
                             command=self.bv_2.get)
        self.cb2.grid(row=3, column=0, sticky='nsew')
        self.cb2.bind("<ButtonRelease>", self.bv_get)
        self.cb3 = tk.Checkbutton(self.LF1, text='115deg', variable=self.bv_3,
                             command=self.bv_3.get)
        self.cb3.grid(row=4, column=0, sticky='nsew')
        self.cb3.bind("<ButtonRelease>", self.bv_get)
        self.cb4 = tk.Checkbutton(self.LF1, text='135deg', variable=self.bv_4,
                             command=self.bv_4.get)
        self.cb4.grid(row=5, column=0, sticky='nsew')
        self.cb4.bind("<ButtonRelease>", self.bv_get)
        self.v = tk.StringVar()
        self.v.set(self.value)
        self.LF2 = ttk.LabelFrame(self.button_frame, text='Reaction Time')
        self.scale = ttk.Scale(self.LF2, from_=0, to=1, command=lambda s: self.v.set('Time: %0.0d s' % float(s)))
        self.scale.bind("<ButtonRelease-1>", self.update_value)
        self.LF2.grid(row=5, column=0, columnspan=1, sticky='nsew')
        self.scale.grid(row=0, column=0, sticky='nsew')
        self.val_label = ttk.Label(self.LF2, textvariable=self.v)
        self.val_label.grid(row=1, column=0, sticky='nsew')
        self.LF3 = ttk.LabelFrame(self.button_frame, text='Reaction Details')
        self.LF3.grid(row=6, column=0, columnspan=1, sticky='nsew')

        self.v1 = tk.StringVar()
        self.v1.set("")
        self.v2 = tk.StringVar()
        self.v2.set("")
        self.v3 = tk.StringVar()
        self.v3.set("")
        self.v4 = tk.StringVar()
        self.v4.set("")
        self.v5 = tk.StringVar()
        self.v5.set("")
        self.v6 = tk.StringVar()
        self.v6.set("")
        self.v7 = tk.StringVar()
        self.v7.set("")
        self.v8 = tk.StringVar()
        self.v8.set("")
        self.v9 = tk.StringVar()
        self.v9.set("")
        self.v10 = tk.StringVar()
        self.v10.set("")

        self.CONC = ttk.Label(self.LF3, text='Rector Concentration:').grid(row=0, column=0, sticky='nsew')
        self.CONC_v = ttk.Label(self.LF3, textvariable=self.v1)
        self.CONC_v.grid(row=1, column=0, sticky='nsew')
        self.DF = ttk.Label(self.LF3, text='Dilution Factor:').grid(row=2, column=0, sticky='nsew')
        self.DF_v = ttk.Label(self.LF3, textvariable=self.v2)
        self.DF_v.grid(row=3, column=0, sticky='nsew')
        self.OC = ttk.Label(self.LF3, text='Optical Constant:').grid(row=4, column=0, sticky='nsew')
        self.OC_v = ttk.Label(self.LF3, textvariable=self.v3)
        self.OC_v.grid(row=5, column=0, sticky='nsew')
        self.Dark = ttk.Label(self.LF3, text='Toluene_Dark:').grid(row=6, column=0, sticky='nsew')
        self.Dark_v = ttk.Label(self.LF3, textvariable=self.v4)
        self.Dark_v.grid(row=7, column=0, sticky='nsew')
        self.Toluene = ttk.Label(self.LF3, text='Toluene_90deg:').grid(row=8, column=0, sticky='nsew')
        self.Toluene_v = ttk.Label(self.LF3, textvariable=self.v5)
        self.Toluene_v.grid(row=9, column=0, sticky='nsew')
        self.NF = ttk.Label(self.LF3, text='Normalization Factor:').grid(row=10, column=0, sticky='nsew')
        self.NF_v1 = ttk.Label(self.LF3, textvariable=self.v6)
        self.NF_v1.grid(row=11, column=0, sticky='nsew')
        self.NF_v2 = ttk.Label(self.LF3, textvariable=self.v7)
        self.NF_v2.grid(row=12, column=0, sticky='nsew')
        self.NF_v3 = ttk.Label(self.LF3, textvariable=self.v8)
        self.NF_v3.grid(row=13, column=0, sticky='nsew')
        self.NF_v4 = ttk.Label(self.LF3, textvariable=self.v9)
        self.NF_v4.grid(row=14, column=0, sticky='nsew')
        self.NF_v5 = ttk.Label(self.LF3, textvariable=self.v10)
        self.NF_v5.grid(row=15, column=0, sticky='nsew')
        self.export = ttk.Button(self.button_frame, text='Export Data', command=self.export_data, width=20)
        self.export.grid(row=7, column=0, sticky='nsew')
        self.f = Figure(figsize=(10, 6))
        self.gs = plt.GridSpec(20, 20)
        self.zimm = self.f.add_subplot(self.gs[1:10, 0:9])
        self.sls = self.f.add_subplot(self.gs[1:10, 10:19])
        self.visc = self.f.add_subplot(self.gs[10:19, 0:9])
        self.UV = self.f.add_subplot(self.gs[10:19, 10:19])
        self.checkbox_1 = self.f.add_subplot(self.gs[1:8, 10:14])
        self.checkbox_2 = self.f.add_subplot(self.gs[10:17, 15:19])
        self.checkbox_1.axis('off')
        self.checkbox_2.axis('off')
        self.checkbox_1.patch.set_alpha(0.4)
        self.checkbox_2.patch.set_alpha(0.4)
        self.checkbox_1.tick_params(axis='x', which='both', bottom=False, left=False, labelbottom=False)
        self.checkbox_1.tick_params(axis='y', which='both', bottom=False, left=False, labelleft=False)
        self.checkbox_2.tick_params(axis='x', which='both', bottom=False, left=False, labelbottom=False)
        self.checkbox_2.tick_params(axis='y', which='both', bottom=False, left=False, labelleft=False)
        self.zimm.set_autoscale_on(True)
        self.sls.set_autoscale_on(True)
        self.visc.set_autoscale_on(True)
        self.UV.set_autoscale_on(True)
        self.zimm.relim()
        self.sls.relim()
        self.visc.relim()
        self.UV.relim()
        self.zimm.autoscale_view(True, True, True)
        self.sls.autoscale_view(True, True, True)
        self.visc.autoscale_view(True, True, True)
        self.UV.autoscale_view(True, True, True)
        self.zimm_plot, = self.zimm.plot([], [])
        self.zimm.set_title('Zimm plot')
        self.zimm.set_xlabel('q^2')
        self.zimm.set_ylabel('k/I')
        self.sls_plot, = self.sls.plot([], [])
        self.sls.set_title('Light Scattering Data')
        self.sls.set_xlabel('Time(s)')
        self.sls.set_ylabel('Intensity')
        self.visc_plot, = self.visc.plot([], [])
        self.visc.set_title('Specific Viscosity Data')
        self.visc.set_xlabel('Time(s)')
        self.visc.set_ylabel('Voltage')
        self.UV_plot, = self.UV.plot([], [])
        self.UV.set_title('UV Data')
        self.UV.set_xlabel('Time(s)')
        self.UV.set_ylabel('Absorption')
        self.canvas = FigureCanvasTkAgg(self.f, self.canvas_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.canvas_frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.LEFT, fill="both", expand=True)
        self.canvas.draw()
        self.f.tight_layout()
        self.canvas.get_tk_widget().pack(side=tk.LEFT, fill="both", expand=True)

    def load_file(self):
        self.filename = filedialog.askopenfilename(initialdir='C:/Users/Desktop',
                                              title='Please select a ACOMP-Ignition Data file')
        df = pd.read_csv(self.filename, header=0, low_memory=False) #pandas closes file after reading
        dm = df.fillna(0)
        self.dm = dm
        self.v1.set('%0.4f g/ml'% float(dm.at[0, 'init_batch_mono_conc']))
        self.v2.set('%0.2f times'% float(dm.at[0, 'sls_dilution_factor']))
        self.v3.set('%0.5e'% float(dm.at[0, 'sls_optical_factor']))
        self.v4.set('%0.0f'% float(dm.at[0, 'sls_dark_0']))
        self.v5.set('%0.0f'% float(dm.at[0, 'sls_toluene_0']))
        self.v6.set('sls_angle_0: %0.3f'% float(dm.at[0, 'angle_normalization_constant_0']))
        self.v7.set('sls_angle_1: %0.3f'% float(dm.at[0, 'angle_normalization_constant_1']))
        self.v8.set('sls_angle_2: %0.3f'% float(dm.at[0, 'angle_normalization_constant_2']))
        self.v9.set('sls_angle_3: %0.3f'% float(dm.at[0, 'angle_normalization_constant_3']))
        self.v10.set('sls_angle_4: %0.3f'% float(dm.at[0, 'angle_normalization_constant_4']))
        sls_col = [f'sls_midpoint_{i}' for i in range(5)]
        df_nominator = dm.at[0, 'sls_toluene_0'] - dm.at[0, 'sls_dark_0']
        solvent = dm[dm['reaction_trigger']==2].index[0]
        if solvent >= 500:
            s = solvent-500
            e = solvent
        else:
            s = 0
            e = solvent
        self.sls_process = (dm[sls_col].rolling(window=29).median() - dm[sls_col][s:e].mean())*1.19e-5 * 0.98 / df_nominator  # this needs to be modified eventually (automatic selecting)
        for i in range(5):
            self.sls_process[sls_col[i]] *= (dm[f'angle_normalization_constant_{i}'][0])
        self.q_sqr = [(4 * np.pi * 1.33 * np.sin(np.radians(i/2)) / 6.6e-7) ** 2 for i in (45, 65, 90, 115, 135)]
        self.lambda_0 = str(int(dm.at[0, 'uv_wavelength_0_pv']))+'nm'
        self.lambda_1 = str(int(dm.at[0, 'uv_wavelength_1_pv']))+'nm'
        self.lambda_2 = str(int(dm.at[0, 'uv_wavelength_2_pv']))+'nm'
        self.lambda_3 = str(int(dm.at[0, 'uv_wavelength_3_pv']))+'nm'
        uv_col = [f'uv_absorbance_{i}' for i in range(4)]
        self.uv_process = (df[uv_col].rolling(window=29).median() - df[uv_col][s:e].mean()) * 1e-6
        self.uv_max = self.uv_process.max()
        self.conversion = 1 - self.uv_process/self.uv_max
        visc_col =[f'visc_dpt']
        self.visc_process = (df[visc_col].rolling(window=29).median() - df[visc_col][s:e].mean()) / df[visc_col][s:e].mean()

        self.f.suptitle(self.filename, x=0.2, y=1, horizontalalignment='left', verticalalignment='top', fontsize=10)
        self.sls.cla()
        self.sls.set_title('Light Scattering Data')
        self.sls.set_xlabel('Time(s)')
        self.sls.set_ylabel('Intensity')
        self.checkbox_1.cla()
        self.checkbox_1.patch.set_alpha(0.0)
        self.checkbox_1.axis('off')
        self.UV.cla()
        self.UV.set_title('UV Data')
        self.UV.set_xlabel('Time(s)')
        self.UV.set_ylabel('Absorption')
        self.checkbox_2.cla()
        self.checkbox_2.patch.set_alpha(0.0)
        self.checkbox_2.axis('off')
        self.visc.cla()
        self.visc.set_title('Specific Viscosity Data')
        self.visc.set_xlabel('Time(s)')
        self.visc.set_ylabel('')
        self.zimm.cla()
        self.zimm.set_title('Zimm plot')
        self.zimm.set_xlabel('q^2')
        self.zimm.set_ylabel('k/I')
        self.canvas.draw()

    def data_plot(self):
        (self.m, self.n) = self.dm.shape
        self.scale = ttk.Scale(self.LF2, from_=0, to=self.m-1, command=lambda s: self.v.set('Time: %0.0d s' % float(s)))
        self.scale.bind("<ButtonRelease-1>", self.update_value)
        self.scale.grid(row=0, column=0, sticky='nsew')
        self.val_label = ttk.Label(self.LF2, textvariable=self.v)
        self.val_label.grid(row=1, column=0, sticky='nsew')
        self.sls.cla()
        self.checkbox_1.cla()
        self.sls.set_title('Light Scattering Data')
        ls0, = self.sls.plot(self.dm['Index'], self.sls_process['sls_midpoint_0'])
        ls1, = self.sls.plot(self.dm['Index'], self.sls_process['sls_midpoint_1'])
        ls2, = self.sls.plot(self.dm['Index'], self.sls_process['sls_midpoint_2'])
        ls3, = self.sls.plot(self.dm['Index'], self.sls_process['sls_midpoint_3'])
        ls4, = self.sls.plot(self.dm['Index'], self.sls_process['sls_midpoint_4'])
        self.check_1 = CheckButtons(self.checkbox_1, ('45$^\circ$', '65$^\circ$', '90$^\circ$', '115$^\circ$', '135$^\circ$'), (True, True, True, True, True))
        self.checkbox_1.patch.set_alpha(0.0)
        self.checkbox_1.axis('off')

        def func_1(label):
            if label == '45$^\circ$':
                ls0.set_visible(not ls0.get_visible())
            elif label == '65$^\circ$':
                ls1.set_visible(not ls1.get_visible())
            elif label == '90$^\circ$':
                ls2.set_visible(not ls2.get_visible())
            elif label == '115$^\circ$':
                ls3.set_visible(not ls3.get_visible())
            elif label == '135$^\circ$':
                ls4.set_visible(not ls4.get_visible())
            a = self.sls.axvline(x=self.value)
            self.canvas.draw()
            a.remove()
        self.check_1.on_clicked(func_1)
        self.sls.set_title('Light Scattering Data')
        self.sls.set_xlabel('Time(s)')
        self.sls.set_ylabel('Intensity')
        self.update()

        self.UV.cla()
        self.checkbox_2.cla()
        self.UV.set_title('UV Data')
        uv0, = self.UV.plot(self.dm['Index'], self.uv_process['uv_absorbance_0'])
        uv1, = self.UV.plot(self.dm['Index'], self.uv_process['uv_absorbance_1'])
        uv2, = self.UV.plot(self.dm['Index'], self.uv_process['uv_absorbance_2'])
        uv3, = self.UV.plot(self.dm['Index'], self.uv_process['uv_absorbance_3'])
        self.check_2 = CheckButtons(self.checkbox_2, (self.lambda_0, self.lambda_1, self.lambda_2, self.lambda_3), (True, True, True, True))
        self.checkbox_2.patch.set_alpha(0.0)
        self.checkbox_2.axis('off')

        def func_2(label):
            if label == self.lambda_0:
                uv0.set_visible(not uv0.get_visible())
            elif label == self.lambda_1:
                uv1.set_visible(not uv1.get_visible())
            elif label == self.lambda_2:
                uv2.set_visible(not uv2.get_visible())
            elif label == self.lambda_3:
                uv3.set_visible(not uv3.get_visible())
            a = self.sls.axvline(x=self.value)
            self.canvas.draw()
            a.remove()
        self.check_2.on_clicked(func_2)
        self.UV.set_title('UV Data')
        self.UV.set_xlabel('Time(s)')
        self.UV.set_ylabel('Absorption')

        self.visc.cla()
        self.visc.set_title('Specific Viscosity Data')
        self.visc.set_xlabel('Time(s)')
        self.visc.set_ylabel('')
        self.visc.plot(self.dm['Index'], self.visc_process['visc_dpt'])
        self.visc.set_ylim(0)

        self.canvas.draw()
        self.update()

    def update_value(self, event):
        self.value = int(self.scale.get())
        a = self.sls.axvline(x=self.value)
        NF= self.sls_process.iloc[self.value, 2]/self.sls_process.iloc[self.value]
        print(NF)
        self.I_theta = self.dm.at[0, 'sls_optical_factor'] / self.sls_process.iloc[self.value]
        self.angle = [self.bv_0.get(), self.bv_1.get(), self.bv_2.get(), self.bv_3.get(), self.bv_4.get()]
        self.dot = [(self.q_sqr[i], self.I_theta[i]) for i in range(5) if self.angle[i] == True]
        xdata = [self.dot[i][0] for i in range(len(self.dot))]
        ydata = [self.dot[i][1] for i in range(len(self.dot))]
        self.zimm.cla()
        self.zimm.set_title('Zimm plot')
        self.zimm.set_xlabel('q^2')
        self.zimm.set_ylabel('k/I')
        self.zimm.scatter(xdata, ydata, c='r')
        if len(self.dot) > 1:
            k, b = np.polyfit(xdata, ydata, 1)
            xdata1 = np.array(xdata)
            ydata1 = np.array(k * xdata1 + b)
            self.zimm.plot(xdata1, ydata1, c='g')
        self.canvas.draw()
        a.remove()
        self.update()

    def bv_get(self, event):
        self.after(100, self.update_zimm_plot)

    def update_zimm_plot(self):
        self.value = int(self.scale.get())
        a = self.sls.axvline(x=self.value)
        self.I_theta = [f'I{i}' for i in range(5)]
        self.I_theta = self.dm.at[0, 'sls_optical_factor'] / self.sls_process.iloc[self.value]
        self.angle = [self.bv_0.get(), self.bv_1.get(), self.bv_2.get(), self.bv_3.get(), self.bv_4.get()]
        self.dot = [(self.q_sqr[i], self.I_theta[i]) for i in range(5) if self.angle[i] == True]
        xdata = [self.dot[i][0] for i in range(len(self.dot))]
        ydata = [self.dot[i][1] for i in range(len(self.dot))]
        self.zimm.cla()
        self.zimm.set_title('Zimm plot')
        self.zimm.set_xlabel('q^2')
        self.zimm.set_ylabel('k/I')
        self.zimm.scatter(xdata, ydata, c='r')
        if len(self.dot) > 1:
            k, b = np.polyfit(xdata, ydata, 1)
            xdata1 = np.array(xdata)
            ydata1 = np.array(k * xdata1 + b)
            self.zimm.plot(xdata1, ydata1, c='g')
        self.canvas.draw()
        a.remove()
        self.update()
    def export_data(self):
        self.export_file = pd.concat([self.dm['Index'], self.sls_process, self.uv_process, self.visc_process], axis=1)
        self.export_file.to_csv(self.filename+'-processed.csv', index=False, header=1)
        messagebox.showinfo("Message", "Data file has been created")

class RockwellACOMP(tk.Frame):
    def __init__(self, container, controller):
        tk.Frame.__init__(self, container)
        self.grid(row=0, column=0, sticky='nsew')
        self.t = 1000
        self.button_frame = ttk.Frame(self)
        self.button_frame.grid(row=0, column=0, sticky='nsew')
        self.loadbutton = ttk.Button(self.button_frame, text='Load File',
                                    command=self.load_file, width=15)
        self.loadbutton.grid(row=0, column=0, sticky='nsew')
        self.plotbutton = ttk.Button(self.button_frame, text='Start Zimm Plot',
                                    command=self.update_plot, width=15)
        self.plotbutton.grid(row=1, column=0, sticky='nsew')
        self.plotrestbutton = ttk.Button(self.button_frame, text='Reset Zimm Plot',
                                        command=self.reset_plot, width=15)
        self.plotrestbutton.grid(row=2, column=0, sticky='nsew')
        self.slsbutton = ttk.Button(self.button_frame, text='Plot MALS data',
                                   command=self.MALS_plot, width=15)
        self.slsbutton.grid(row=3, column=0, sticky='nsew')
        self.continue_plotting = False
        self.NF1 = ttk.Label(self.button_frame, text="Normalization Factor_45deg").grid(row=5, column=0, sticky='nsew')
        self.NF2 = ttk.Label(self.button_frame, text="Normalization Factor_65deg").grid(row=6, column=0, sticky='nsew')
        self.NF3 = ttk.Label(self.button_frame, text="Normalization Factor_90deg").grid(row=7, column=0, sticky='nsew')
        self.NF4 = ttk.Label(self.button_frame, text="Normalization Factor_15deg").grid(row=8, column=0, sticky='nsew')
        self.NF5 = ttk.Label(self.button_frame, text="Normalization Factor_135deg").grid(row=9, column=0, sticky='nsew')
        self.NF5 = ttk.Label(self.button_frame, text="Toluene_90deg Value").grid(row=10, column=0, sticky='nsew')
        self.NF1_text = tk.StringVar()
        self.NF1 = tk.Entry(self.button_frame, textvariable=self.NF1_text, width=20)
        self.NF1_text.set('')
        self.NF1.grid(row=5, column=1, stick='nsew')
        self.NF2_text = tk.StringVar()
        self.NF2 = tk.Entry(self.button_frame, textvariable=self.NF2_text, width=20)
        self.NF2_text.set('')
        self.NF2.grid(row=6, column=1, stick='nsew')
        self.NF3_text = tk.StringVar()
        self.NF3 = tk.Entry(self.button_frame, textvariable=self.NF3_text, width=20)
        self.NF3_text.set('')
        self.NF3.grid(row=7, column=1, stick='nsew')
        self.NF4_text = tk.StringVar()
        self.NF4 = tk.Entry(self.button_frame, textvariable=self.NF4_text, width=20)
        self.NF4_text.set('')
        self.NF4.grid(row=8, column=1, stick='nsew')
        self.NF5_text = tk.StringVar()
        self.NF5 = tk.Entry(self.button_frame, textvariable=self.NF5_text, width=20)
        self.NF5_text.set('')
        self.NF5.grid(row=9, column=1, stick='nsew')
        self.TV90_text = tk.StringVar()
        self.TV90 = tk.Entry(self.button_frame, textvariable=self.TV90_text, width=20)
        self.TV90_text.set('')
        self.TV90.grid(row=10, column=1, stick='nsew')
        self.confirmbutton = ttk.Button(self.button_frame, text='confirm',
                                        command=self.on_click, width=15)
        self.confirmbutton.grid(row=11, column=1, sticky='nsew')
        self.back = ttk.Button(self.button_frame, text='Back the Home',
                               command=lambda: controller.show_frame(HomePage), width=15)
        self.back.grid(row=17, column=0, sticky='nsew')
        self.bv_0 = tk.BooleanVar()
        self.bv_0.set(True)
        self.bv_1 = tk.BooleanVar()
        self.bv_1.set(True)
        self.bv_2 = tk.BooleanVar()
        self.bv_2.set(True)
        self.bv_3 = tk.BooleanVar()
        self.bv_3.set(True)
        self.bv_4 = tk.BooleanVar()
        self.bv_4.set(True)
        self.cb0 = ttk.Checkbutton(self.button_frame, text='45deg', variable=self.bv_0,
                                  command=self.bv_0.get).grid(row=0, column=1, sticky='nsew')
        self.cb1 = ttk.Checkbutton(self.button_frame, text='65deg', variable=self.bv_1,
                                  command=self.bv_1.get).grid(row=1, column=1, sticky='nsww')
        self.cb2 = ttk.Checkbutton(self.button_frame, text='90deg', variable=self.bv_2,
                                  command=self.bv_2.get).grid(row=2, column=1, sticky='nsew')
        self.cb3 = ttk.Checkbutton(self.button_frame, text='115deg', variable=self.bv_3,
                                  command=self.bv_3.get).grid(row=3, column=1, sticky='nsew')
        self.cb4 = ttk.Checkbutton(self.button_frame, text='135deg', variable=self.bv_4,
                                  command=self.bv_4.get).grid(row=4, column=1, sticky='nsew')
        self.f = Figure()
        self.a = self.f.add_subplot(121)
        self.sls = self.f.add_subplot(122)
        self.a.set_autoscale_on(True)
        self.sls.set_autoscale_on(True)
        self.a.relim()
        self.sls.relim()
        self.a.autoscale_view(True, True, True)
        self.sls.autoscale_view(True, True, True)
        self.zimm_plot, = self.a.plot([], [])
        self.a.set_title('Zimm plot')
        self.a.set_xlabel('q^2')
        self.a.set_ylabel('1/I')
        self.sls_plot, = self.sls.plot([], [])
        self.sls.set_title('Light Scattering Data')
        self.canvas = FigureCanvasTkAgg(self.f, self)
        self.canvas.draw()
        self.f.tight_layout()
        self.canvas.get_tk_widget().grid(row=0, column=1, sticky='nsew')

    def on_click(self):
        from tkinter import messagebox
        NF1 = self.NF1_text.get()
        NF2 = self.NF2_text.get()
        NF3 = self.NF3_text.get()
        NF4 = self.NF4_text.get()
        NF5 = self.NF5_text.get()
        TV90 = self.TV90_text.get()
        Input="Normalization Factor_45deg: %s\n Normalization Factor_65deg: %s\n Normalization Factor_90deg: %s\n " \
              "Normalization Factor_115deg: %s\n Normalization Factor_135deg: %s\n Toluene 90deg Value: %s" \
              %(NF1, NF2, NF3, NF4, NF5, TV90)
        messagebox.showinfo("Make sure the normalization factor is correct", Input)
    def load_file(self):
        filename = filedialog.askopenfilename(initialdir='C:/Users/Desktop',
                                              title='Please select a ACOMP-Ignition Data file')
        df = pd.read_csv(filename, header=0)  # pandas closes file after reading
        dm = df.fillna(0)
        self.dm = dm
        sls_col = [f'D3_LS_Region_midpoint{i}' for i in range(5)]
        df_nominator = dm.at[0, 'sls_toluene_0'] - float(self.TV90_text.get())
        # baseline=
        self.sls_process = (dm[sls_col].rolling(window=29).median() - dm[sls_col][
                                                                      1000:1500].mean()) * 1.19e-5 * 0.98 / df_nominator
        # this needs to be modified eventually (automatic selecting)
        for i in range(5):
            self.sls_process[sls_col[i]] *= float(f'self.NF{i}.get()')
        self.q_sqr = [(4 * np.pi * 1.33 * np.sin(np.radians(i / 2)) / 6.6e-7) ** 2 for i in (45, 65, 90, 115, 135)]
        uv_col = [f'D1_UV_abs{i}' for i in range(4)]
        self.uv_process = (df[uv_col].rolling(window=29).median() - df[uv_col][1000:1500].mean()) * 1e-6
        self.t = 0
        self.plotbutton.config(text='Start Zimm plot')
    def reset_plot(self):
        self.t = 500
        self.a.cla()
    def MALS_plot(self):
        self.sls.cla()
        self.sls.set_title('Light Scattering Data')
        ls0, = self.sls.plot(self.dm['Index'], self.sls_process['D3_LS_Region_midpoint0'])
        ls1, = self.sls.plot(self.dm['Index'], self.sls_process['D3_LS_Region_midpoint1'])
        ls2, = self.sls.plot(self.dm['Index'], self.sls_process['D3_LS_Region_midpoint2'])
        ls3, = self.sls.plot(self.dm['Index'], self.sls_process['D3_LS_Region_midpoint3'])
        ls4, = self.sls.plot(self.dm['Index'], self.sls_process['D3_LS_Region_midpoint4'])
        self.check = CheckButtons(self.sls, ('45deg', '65deg', '90deg', '115deg', '135deg'),
                                  (True, True, True, True, True))
        def func(label):
            if label == '45deg':
                ls0.set_visible(not ls0.get_visible())
            elif label == '65deg':
                ls1.set_visible(not ls1.get_visible())
            elif label == '90deg':
                ls2.set_visible(not ls2.get_visible())
            elif label == '115deg':
                ls3.set_visible(not ls3.get_visible())
            elif label == '135deg':
                ls4.set_visible(not ls4.get_visible())
            self.canvas.draw()
        self.check.on_clicked(func)
        self.canvas.draw()
    def update_data(self):
        self.I_theta = [f'I{i}' for i in range(5)]
        self.I_theta = 1 / self.sls_process.iloc[self.t]
        self.angle = [self.bv_0.get(), self.bv_1.get(), self.bv_2.get(), self.bv_3.get(), self.bv_4.get()]
        self.dot = [(self.q_sqr[i], self.I_theta[i]) for i in range(5) if self.angle[i] == True]
        self.t += 1
    def update_plot(self):
        (m, n) = self.sls_process.shape
        self.continue_plotting = not self.continue_plotting
        if self.continue_plotting:
            self.plotbutton.config(text='Pause Zimm plot')
        else:
            self.plotbutton.config(text='Continue Zimm plot')
        # while self.continue_plotting == True and self.t < m:
            self.after(50, self.update_data())
            xdata = [self.dot[i][0] for i in range(len(self.dot))]
            ydata = [self.dot[i][1] for i in range(len(self.dot))]
            self.a.cla()
            self.a.set_title('Zimm plot')
            self.a.set_xlabel('q^2')
            self.a.set_ylabel('1/I')
            self.a.scatter(xdata, ydata, c='r')
            if len(self.dot) > 1:
                k, b = np.polyfit(xdata, ydata, 1)
                xdata1 = np.array(xdata)
                ydata1 = np.array(k * xdata1 + b)
                self.a.plot(xdata1, ydata1, c='g')
            self.canvas.draw()
            self.update()
class MWUV(tk.Frame):
    def __init__(self, container, controller):
        tk.Frame.__init__(self, container)
        self.grid(row=0, column=0, sticky='nsew')
        self.button_frame = tk.Frame(self)
        self.button_frame.grid(row=0, column=0, sticky='nsew')

        self.loadECF = ttk.Button(self.button_frame, text='Load Extinction Coefficient Spectra File',
                                     command=self.load_extinction_coefficient_file, width=15)
        self.loadECF.grid(row=8, column=0, sticky='nsew')

        self.loadRF = ttk.Button(self.button_frame, text='Load Reaction File',
                                  command=self.load_reaction_file, width=15)
        self.loadRF.grid(row=9, column=0, sticky='nsew')

        self.l1 = ttk.Label(self.button_frame, text='Mass of Styrene Sulfonate(g)')
        self.l1.grid(row=0, column=0, sticky='nsew')
        self.SS_text = tk.StringVar()
        self.SS = tk.Entry(self.button_frame, textvariable=self.SS_text)
        self.SS_text.set('')
        self.SS.grid(row=0, column=1, sticky='nsew')

        self.l2 = ttk.Label(self.button_frame, text='Mass of Acrylamide(g)')
        self.l2.grid(row=1, column=0, sticky='nsew')
        self.AM_text = tk.StringVar()
        self.AM = tk.Entry(self.button_frame, textvariable=self.AM_text)
        self.AM_text.set('')
        self.AM.grid(row=1, column=1, sticky='nsew')

        self.l3 = ttk.Label(self.button_frame, text='Total Volume(ml)')
        self.l3.grid(row=2, column=0, sticky='nsew')
        self.V_text = tk.StringVar()
        self.V = tk.Entry(self.button_frame, textvariable=self.V_text)
        self.V_text.set('')
        self.V.grid(row=2, column=1, sticky='nsew')

        self.l4 = ttk.Label(self.button_frame, text='Dilution Factor')
        self.l4.grid(row=3, column=0, sticky='nsew')
        self.DF_text = tk.StringVar()
        self.DF = tk.Entry(self.button_frame, textvariable=self.DF_text)
        self.DF_text.set('')
        self.DF.grid(row=3, column=1, sticky='nsew')

        self.l5 = ttk.Label(self.button_frame, text='Starting Wavelength(nm)')
        self.l5.grid(row=4, column=0, sticky='nsew')
        self.L1_text = tk.StringVar()
        self.L1 = tk.Entry(self.button_frame, textvariable=self.L1_text)
        self.L1_text.set('')
        self.L1.grid(row=4, column=1, sticky='nsew')

        self.l6 = ttk.Label(self.button_frame, text='Ending Wavelength(nm)')
        self.l6.grid(row=5, column=0, sticky='nsew')
        self.L2_text = tk.StringVar()
        self.L2 = tk.Entry(self.button_frame, textvariable=self.L2_text)
        self.L2_text.set('')
        self.L2.grid(row=5, column=1, sticky='nsew')

        self.l7 = ttk.Label(self.button_frame, text="Styrene Sulfonate Concentration")
        self.l7.grid(row=3, column=2, sticky='nsew')
        self.ss_answerbox = tk.Entry(self.button_frame)
        self.ss_answerbox.grid(row=3, column=3, sticky='nsew')

        self.l7 = ttk.Label(self.button_frame, text="Styrene Sulfonate Concentration")
        self.l7.grid(row=4, column=2, sticky='nsew')
        self.am_answerbox = tk.Entry(self.button_frame)
        self.am_answerbox.grid(row=4, column=3, sticky='nsew')

        self.confirmbutton = ttk.Button(self.button_frame, text='Input', command=self.on_click)
        self.confirmbutton.grid(row=6, column=0, sticky='nsew')

        self.calculate_button = ttk.Button(self.button_frame, text='Calculate', command=self.calcualtion)
        self.calculate_button.grid(row=7, column=0, sticky='nsew')

        self.back = ttk.Button(self.button_frame, text='Back the Home',
                              command=lambda: controller.show_frame(HomePage), width=15)
        self.back.grid(row=10, column=0, sticky='nsew')

    def on_click(self):
        from tkinter import messagebox
        SS_mass = self.SS_text.get()
        AM_mass = self.AM_text.get()
        Volume = self.V_text.get()
        Dilution = self.DF_text.get()
        Lambda_s = self.L1_text.get()
        Lambda_e = self.L2_text.get()
        Input = "Mss: %s grams\n Mam: %s grams\n Volume: %s ml\n DF=%s\n Wavelength from %s nm to %s nm" % (
        SS_mass, AM_mass, Volume, Dilution, Lambda_s, Lambda_e)
        messagebox.showinfo("Make sure the reaction information is correct", Input)

    def calcualtion(self):
        self.mss = float(self.SS_text.get())
        self.mam = float(self.AM_text.get())
        self.v = float(self.V_text.get())
        self.dfx = float(self.DF_text.get())
        self.l1 = int(self.L1_text.get())
        self.l2 = int(self.L2_text.get())
        self.conc_ss = self.mss / (self.v * self.dfx)
        self.conc_am = self.mam / (self.v * self.dfx)

        self.ss_answerbox.insert(0, self.conc_ss)
        self.am_answerbox.insert(0, self.conc_am)

    def load_extinction_coefficient_file(self):
        filename1 = filedialog.askopenfilename(initialdir='C:/Users/Desktop',
                                          title='Please Select Extinction Coefficient Spectra File')
        self.DataFrame = pd.read_csv(filename1, header=0)  # Starting line needs to change
        self.DF1 = self.DataFrame.fillna(0)
        (m1, n1) = self.DataFrame.shape
    # Find Corresponding Wavelength Columns
        self.x1 = self.l1+4
        self.x2 = self.l2+5
    # Extinction of Polymers for Selected Range
        self.ep1 = self.DF1.iloc[2, self.x1:self.x2]
        self.ep2 = self.DF1.iloc[3, self.x1:self.x2]
    # Extinction Coefficient Difference Between Polymers and Monomers
        self.delta1 = self.DF1.iloc[0, self.x1:self.x2] - self.DF1.iloc[2, self.x1:self.x2]  # monomer1 subtract polymer1
        self.delta2 = self.DF1.iloc[1, self.x1:self.x2] - self.DF1.iloc[3, self.x1:self.x2]  # monomer2 subtract polymer2
        self.a11 = self.delta1.dot(self.delta1.T)
        self.a12 = self.delta1.dot(self.delta2.T)
        self.a22 = self.delta2.dot(self.delta2.T)
        self.Aa1 = self.ep1 * self.conc_ss
        self.Aa2 = self.ep2 * self.conc_am
        self.AA = self.Aa1 + self.Aa2
        self.B = np.array([[self.a11, self.a12], [self.a12, self.a22]])
        print(self.B)
        self.B_invert = np.linalg.inv(self.B)
    def load_reaction_file(self):
        filename2 = filedialog.askopenfilename(initialdir='C:/Users/Desktop', title='Please Select Reaction Data File')
        self.DataFrame2 = pd.read_csv(filename2, header=9)
        self.DF2 = self.DataFrame2.fillna(0)
        (m2, n2) = self.DF2.shape
    # Calculate the Sigma Absorption of Defined UV range
        self.sum_array = np.zeros((m2, 1))
        self.A = np.zeros((m2, 2))
        self.Conc = np.zeros((m2, 2))
        for i in range(0, m2):
            self.sum_array[i] = self.DF2.iloc[i, self.x1:self.x2].sum()
            self.A[i, 0] = self.DF2.iloc[i, self.x1:self.x2].dot(self.delta1) - self.AA.dot(self.delta1)
            self.A[i, 1] = self.DF2.iloc[i, self.x1:self.x2].dot(self.delta2) - self.AA.dot(self.delta2)
            self. Conc[i] = np.array([self.A[i, 0], self.A[i, 1]]).dot(self.B_invert)
class Others(tk.Frame):
    def __init__(self, container, controller):
        tk.Frame.__init__(self, container)
        # self.root = root
        self.grid(row=0, column=0, sticky='nsew')
        self.t = 1000
        self.button_frame = tk.Frame(self)
        self.button_frame.grid(row=0, column=0, sticky='nsew')
        self.loadbutton1 = ttk.Button(self.button_frame, text='Light Scattering Data File',
                                      command=self.bokeh, width=20)
        self.loadbutton1.grid(row=0, column=0, sticky='nsew')
        self.loadbutton2 = ttk.Button(self.button_frame, text='Load MC-DAQ File',
                                    command=self.mc_daq, width=20)
        self.loadbutton2.grid(row=1, column=0, sticky='nsew')
        self.loadbutton3 = ttk.Button(self.button_frame, text='Viscosity/RI',
                                    command=self.rv, width=20)
        self.loadbutton3.grid(row=2, column=0, sticky='nsew')
        self.save_button = ttk.Button(self.button_frame, text='Merge and Save Date Files',
                                      command=self.merge, width=20)
        self.save_button.grid(row=3, column=0, sticky='nsew')
        self.back = tk.Button(self.button_frame, text='Back the Home',
                              command=lambda: controller.show_frame(HomePage), width=20)
        self.back.grid(row=4, column=0, sticky='nsew')
        self.bv_0 = tk.BooleanVar()
        self.bv_0.set(True)
        self.bv_1 = tk.BooleanVar()
        self.bv_1.set(True)
        self.bv_2 = tk.BooleanVar()
        self.bv_2.set(True)
        self.bv_3 = tk.BooleanVar()
        self.bv_3.set(True)
        self.bv_4 = tk.BooleanVar()
        self.bv_4.set(True)
        self.cb0 = ttk.Checkbutton(self.button_frame, text='45deg', variable=self.bv_0,
                                  command=self.bv_0.get).grid(row=0, column=1, sticky='nsew')
        self.cb1 = ttk.Checkbutton(self.button_frame, text='65deg', variable=self.bv_1,
                                  command=self.bv_1.get).grid(row=1, column=1, sticky='nsww')
        self.cb2 = ttk.Checkbutton(self.button_frame, text='90deg', variable=self.bv_2,
                                  command=self.bv_2.get).grid(row=2, column=1, sticky='nsew')
        self.cb3 = ttk.Checkbutton(self.button_frame, text='115deg', variable=self.bv_3,
                                  command=self.bv_3.get).grid(row=3, column=1, sticky='nsew')
        self.cb4 = ttk.Checkbutton(self.button_frame, text='135deg', variable=self.bv_4,
                                  command=self.bv_4.get).grid(row=4, column=1, sticky='nsew')
        f = plt.figure()
        a = f.add_subplot(121)
        sls = f.add_subplot(122)
        a.set_autoscale_on(True)
        sls.set_autoscale_on(True)
        a.relim()
        sls.relim()
        a.autoscale_view(True, True, True)
        sls.autoscale_view(True, True, True)
        self.a = a
        self.sls = sls
        self.zimm_plot, = a.plot([], [])
        self.a.set_title('UV absorption plot')
        self.a.set_xlabel('time(s)')
        self.a.set_ylabel('UV absorption')
        self.sls_plot, = sls.plot([], [])
        self.sls.set_title('Light Scattering Data')
        self.sls.set_xlabel('time(s)')
        self.sls.set_ylabel('Intensity')
        self.canvas = FigureCanvasTkAgg(f, self)
        self.canvas.draw()
        f.tight_layout()
        self.canvas.get_tk_widget().grid(row=0, column=2, sticky='nsew')
    def bokeh(self):
        filename1 = filedialog.askopenfilename(initialdir='C:/Users/Desktop',
                                              title='Please select a Light Scattering Data file')
        df = pd.read_csv(filename1, header=0)  # pandas closes file after reading
        self.dm_LS= df.fillna(0)
    def mc_daq(self):
        filename2 = filedialog.askopenfilename(initialdir='C:/Users/Desktop',
                                              title='Please select a MC-DAQ/UV Data file')
        df_UV = pd.read_table(filename2, header=0)
        self.dm_UV = df_UV.fillna(0)
    def rv(self):
        filename3 = filedialog.askopenfilename(initialdir='C:/Users/Desktop',
                                              title='Please select a Viscosity/RI Data file')
        df = pd.read_csv(filename3, header=6, sep=',', encoding='utf-8')  # pandas closes file after reading
        self.dm_VR = df.fillna(0)
    def merge(self):
        UV = self.dm_UV.iloc[:, 2:6]
        MALS = self.dm_LS.iloc[:, 4:7]
        VR = self.dm_VR.iloc[:, 2:4]
        self.MERGE = pd.concat([UV, MALS, VR], axis=1)
        self.MERGE.columns = ['230nm', '240nm', '260nm', '280nm', '45deg', '90deg', '135deg', 'Viscosity', 'RI_voltage']
        filename = filedialog.asksaveasfilename(initialdir="C:/Users/Desktop", title="create a file", filetypes=[('csv files', ".csv")])
        self.MERGE.to_csv(self.filename+'.csv', index=False, header=1)

if __name__ == "__main__":
    app = Main()
    app.mainloop()

