import matplotlib, sys
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from .ftan_call import pickgroupdispcurv
from msnoise.api import *
import shutil
import glob
import random

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk as NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from tkinter import *
from tkinter import ttk
from msnoise.api import *

try:
    from tkFileDialog import askopenfilename, askdirectory
    from tkFont import Font
except:
    from tkinter.filedialog import askopenfilename, askdirectory
from tkinter.font import Font




def main():
    # PLOTDIAGR = show
    PLOTRAWDISP = 0
    PLOTDISPALL = 0
    SAVEFILES = 0

    if not os.path.isdir("TOMO_DISP"):
        os.makedirs("TOMO_DISP")
    db = connect()
    PER = get_config(db, "ftan_periods", plugin="Tomo")
    PER = np.array([float(pi) for pi in PER.split(',')])
    fmin = float(get_config(db, "ftan_fmin", plugin="Tomo"))
    fmax = float(get_config(db, "ftan_fmax", plugin="Tomo"))
    vgmin = float(get_config(db, "ftan_vgmin", plugin="Tomo"))
    vgmax = float(get_config(db, "ftan_vgmax", plugin="Tomo"))
    bmin = float(get_config(db, "ftan_bmin", plugin="Tomo"))
    bmax = float(get_config(db, "ftan_bmax", plugin="Tomo"))

    diagramtype = get_config(db, "ftan_diagramtype", plugin="Tomo")

    nfreq = int(get_config(db, "ftan_nfreq", plugin="Tomo"))
    ampmin = float(get_config(db, "ftan_ampmin", plugin="Tomo"))
    global data
    data = pd.DataFrame()
    # nfreq = 100
    db = connect()

    def load_dir():
        folder = askdirectory(parent=root)
        files = sorted(glob.glob(os.path.join(os.path.realpath(folder),"*_MEAN.*")))
        cb['values'] = files
        cb_val.set(files[0])
    
    def save():
        global data
        print(data.head())
        filename = cb_val.get()
        filename = filename.replace("TOMO_SAC", "TOMO_DISP").replace(".sac",".csv").replace(".SAC",".csv")
        if not os.path.isdir(os.path.split(filename)[0]):
            os.makedirs(os.path.split(filename)[0])
        data.to_csv(filename)
        
    
    def process(e=None, xdata=0, ydata=0):
        filename = cb_val.get()
        NET1, STA1, NET2, STA2, crap = os.path.split(filename)[1].split('_')

        st = read(filename)
        dist = st[0].stats.sac.dist
        dt = st[0].stats.delta

        p = ccf.gca()
        p.cla()
        taxis = np.arange(st[0].stats.npts)*st[0].stats.delta
        vaxis = dist/taxis[1:]
        zone = np.where((vaxis >= float(_vgmin.get())) & (vaxis <= float(_vgmax.get())))[0]

        p.plot(taxis, st[0].data)
        p.plot(taxis[zone], st[0].data[zone], c='r')
        ccf.subplots_adjust(bottom=0.25)
        ccfcanvas.draw()

        per, disper, seeds = pickgroupdispcurv(filename,
                                               _fmin.get(),
                                               _fmax.get(),
                                               _vgmin.get(),
                                               _vgmax.get(),
                                               _bmin.get(),
                                               _bmax.get(),
                                               diagramtype,
                                               nfreq,
                                               _ampmin.get(),
                                               dist,
                                               xdata,
                                               ydata)
        basename = "%s.%s_%s.%s_%s" % (NET1, STA1, NET2, STA2, crap)
        basename = basename.replace(".SAC", "")
        for _ in ["write_amp.txt",
                  "write_disp.txt",
                  "write_FP.txt",
                  "write_ph.txt",
                  "write_TV.txt",
                  ]:
            shutil.move(_, _.replace("write", basename))

        U = np.loadtxt('%s_TV.txt' % basename)
        P = np.loadtxt('%s_FP.txt' % basename)
        xmin = min(P)
        xmax = max(P)
        ymin = min(U)
        ymax = max(U)
        amp = np.loadtxt('%s_amp.txt' % basename).T
        iu = np.where((disper >= ymin) & (disper <= ymax))
        per = per[iu]
        disper = disper[iu]
        Per, Vitg = np.meshgrid(P, U)
        f.clf()
        p = f.gca()
        p.cla()
        if int(_normed.get()):
            for i in range(amp.shape[1]):
                amp[:,i] /= amp[:,i].max() 
        
        c = p.contourf(Per, Vitg, amp, 35, cmap=cm_val.get())
        f.colorbar(c)
        
        idxx = np.where(float(_minWL.get())*per*disper < dist)[0]
        p.plot(per, disper,'-ok', lw=1.5)
        p.scatter(per[idxx], disper[idxx], marker='o', s=100, c="green", lw=1.5)
        p.scatter(seeds[:,0], seeds[:,1], s=10, marker='d', facecolor='w',
                  edgecolor="k")
        global data
        data = pd.DataFrame(disper[idxx], index=per[idxx], columns=["Velocity"])
        data.index.name = "Period"
        
        # TODO add axes labels depending on diagramtype
        p.set_xlabel("Period (s)", fontsize=12)
        p.set_ylabel('Velocity (km/s)', fontsize=12)
        p.set_xlim(xmin, xmax)
        p.set_ylim(ymin, ymax)
        p.set_title("%s.%s - %s.%s (%.2f km)" % (NET1, STA1, NET2, STA2, dist))
        canvas.draw()

    def previous_file(e=None):
        idx = cb['values'].index(cb.get())
        if idx > 0:
            idx -= 1
        cb.set(cb['values'][idx])
        process()

    def next_file(e=None):
        idx = cb['values'].index(cb.get())
        if idx < len(cb['values']):
            idx += 1
        cb.set(cb['values'][idx])
        process()

    #root
    root = Tk()
    root.title("MSNoise-TOMO: Time-Frequency & Dispersion Curve Picking Tool")
    # root.resizable(True, True)
    # menubar = Menu(root)
    # filemenu = Menu(menubar, tearoff=0)
    # filemenu.add_command(label="Open", command=openfile)
    # filemenu.add_separator()
    # filemenu.add_command(label="Exit", command=root.quit)
    # menubar.add_cascade(label="File", menu=filemenu)

    # root.config(menu=menubar)

    #mainframe
    gui_style = ttk.Style()
    gui_style.configure('My.TLabel', background='#ffffff')
    gui_style.configure('My.TButton', background='#ff0011')
    gui_style.configure('M2.TButton', background='#B7E2F0')
    gui_style.configure('My.TFrame', background='#ffffff')

    mainframe = ttk.Frame(root, padding = "3 3 12 12", style='My.TFrame')
    mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
    mainframe.columnconfigure(0, weight=1)
    mainframe.rowconfigure(0, weight=1)

    _folder = StringVar(root, value=os.path.join(os.getcwd(), "TOMO_SAC"))
    ttk.Button(mainframe, text="Select Folder", command=load_dir,
               style='M2.TButton').grid(column=1, row=1, sticky=W)

    cb_val = StringVar()
    files = sorted(
        glob.glob(os.path.join(os.path.realpath(_folder.get()), "*_MEAN.*")))
    cb = ttk.Combobox(mainframe, width=55, textvariable=cb_val, height=4)
    cb.grid(column=2, row=1, sticky=(W, E), columnspan=2)
    if len(files):
        cb['values'] = files
        cb_val.set(files[0])
    else:
        cb_val.set("<-- SELECT A FOLDER BY CLICKING ON THE BUTTON")
    # PARAMS PANEL
    # myFont = Font(family="Helvetica", size=12)

    _vgmin = StringVar(root, value=vgmin)
    vgminSize = ttk.Entry(mainframe, width=7, textvariable=_vgmin)
    vgminSize.grid(column=2, row=2, sticky=(W, E))
    tmp = ttk.Label(mainframe, text="Vg Min ", style='My.TLabel').grid(column=1, row=2, sticky=W)

    _vgmax = StringVar(root, value=vgmax)
    vgmaxSize = ttk.Entry(mainframe, width=7, textvariable=_vgmax)
    vgmaxSize.grid(column=2, row=3, sticky=(W, E))
    ttk.Label(mainframe, text="Vg Max ", style='My.TLabel').grid(column=1, row=3, sticky=W)

    _minSNR = StringVar(root, value=0.0)
    minSNRSize = ttk.Entry(mainframe, width=7, textvariable=_minSNR)
    minSNRSize.grid(column=2, row=4, sticky=(W, E))
    ttk.Label(mainframe, text="Min SNR ", style='My.TLabel').grid(column=1, row=4, sticky=W)

    _minWL = StringVar(root, value=1.0)
    minWLSize = ttk.Entry(mainframe, width=7, textvariable=_minWL)
    minWLSize.grid(column=2, row=5, sticky=(W, E))
    ttk.Label(mainframe, text="Min Wavelength ", style='My.TLabel').grid(column=1, row=5, sticky=W)

    _fmin = StringVar(root, value=fmin)
    fminSize = ttk.Entry(mainframe, width=7, textvariable=_fmin)
    fminSize.grid(column=2, row=6, sticky=(W, E))
    ttk.Label(mainframe, text="Min Frequency ", style='My.TLabel').grid(column=1, row=6, sticky=W)

    _fmax = StringVar(root, value=fmax)
    fmaxSize = ttk.Entry(mainframe, width=7, textvariable=_fmax)
    fmaxSize.grid(column=2, row=7, sticky=(W, E))
    ttk.Label(mainframe, text="Max Frequency ", style='My.TLabel').grid(column=1, row=7, sticky=W)

    _diagType = StringVar(root, value="PV")
    diagTypeSize = ttk.Entry(mainframe, width=7, textvariable=_diagType)
    diagTypeSize.grid(column=2, row=8, sticky=(W, E))
    ttk.Label(mainframe, text="Diagram (...) ", style='My.TLabel').grid(column=1, row=8, sticky=W)


    _bmin = StringVar(root, value=bmin)
    bminSize = ttk.Entry(mainframe, width = 7, textvariable = _bmin)
    bminSize.grid(column = 2, row = 9, sticky =(W, E))
    ttk.Label(mainframe, text="Bmin ", style='My.TLabel').grid(column=1, row=9, sticky=W)

    _bmax = StringVar(root, value=bmax)
    bminSize = ttk.Entry(mainframe, width = 7, textvariable = _bmax)
    bminSize.grid(column = 2, row = 10, sticky =(W, E))
    ttk.Label(mainframe, text="Bmax ", style='My.TLabel').grid(column=1, row=10, sticky=W)

    _ampmin = StringVar(root, value=ampmin)
    ampminSize = ttk.Entry(mainframe, width = 7, textvariable = _ampmin)
    ampminSize.grid(column = 2, row = 11, sticky =(W, E))
    ttk.Label(mainframe, text="Amp Min ", style='My.TLabel').grid(column=1, row=11, sticky=W)

    maps = sorted(m for m in plt.cm.datad)
    cm_val = StringVar()
    cm = ttk.Combobox(mainframe, width=7, textvariable=cm_val, height=4)
    cm.grid(column=2, row=12, sticky=(W, E))
    ttk.Label(mainframe, text="Cmap ", style='My.TLabel').grid(column=1, row=12,
                                                               sticky=W)
    cm['values'] = maps
    cm_val.set("hot_r")


    ccf = Figure(figsize=(5, 1), dpi=100)
    ccfcanvas = FigureCanvasTkAgg(ccf, master=mainframe)
    ccfcanvas.get_tk_widget().grid(row=2, column=3, rowspan=2)

    f = Figure(dpi=100)
    canvas = FigureCanvasTkAgg(f, master=mainframe)
    canvas.get_tk_widget().grid(row=4, column=3, rowspan=9)

    _normed = IntVar()
    chh = ttk.Checkbutton(mainframe, text="Normed", variable=_normed, \
                     onvalue=1, offvalue=0).grid(column=1, row=13, sticky=W)

    ttk.Button(mainframe, text="Compute", command=process,
               style='My.TButton').grid(column=2, row=13, sticky=W)

    ttk.Button(mainframe, text="Save", command=save,
               style='My.TButton').grid(column=2, row=14, sticky=W)

    # toolbar = NavigationToolbar2TkAgg(canvas, mainframe)
    # toolbar.update()
    # toolbar.grid(column=3, row=14, columnspan=1, sticky=W)
    # canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=0)

    def onclick(event):
        print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              (event.button, event.x, event.y, event.xdata, event.ydata))
        process(event, xdata=event.xdata, ydata=event.ydata)

    cid = f.canvas.mpl_connect('button_press_event', onclick)

    for child in mainframe.winfo_children():
        child.grid_configure(padx=0.5, pady=0.5)

    root.bind('<Return>', process)
    root.bind('<KP_Enter>', process)
    root.bind('<Control-Key-Left>', previous_file)
    root.bind('<Control-Key-Right>', next_file)
    icon = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'img',
                        'msnoise.gif')
    img = PhotoImage(file=icon)

    root.tk.call('wm', 'iconphoto', root._w, img)
    # print(os.path.isfile(icon))
    # root.iconbitmap(icon)

    root.mainloop()
