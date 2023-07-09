from tkinter import *
from tkinter import filedialog
from widgets import *
from colors import *
from gecatsim.pyfiles.CommonTools import *


class XcistScanFrame:
    def __init__(self, window):
        self.frame = XcistFrame(window, "SCAN PROTOCOL", 700, 250)
        self.frame.frame.grid(column=1, columnspan=2, row=1, padx=5, pady=5)
        self.frame.load_btn.configure(command=self.load_protocol)

        self.contents = Frame(self.frame.frame, bg=BG_PRIMARY)
        self.contents.grid(sticky='NSWE', column=0, row=1)
        self.contents.grid_columnconfigure(0, weight=1)
        self.contents.grid_columnconfigure(2, weight=1)
        self.contents.grid_columnconfigure(4, weight=1)
        self.contents.grid_rowconfigure(5, weight=1)

        self.filename_lbl = XcistLabel(self.contents, text="Filename")
        self.filename = XcistNameLabel(self.contents, "\"example\"")
        self.filename_lbl.grid(sticky='W', column=0, columnspan=2, row=0, pady=4)
        self.filename.grid(sticky='W', column=2, columnspan=4, row=0, pady=4)

        self.kvp_label = XcistLabel(self.contents, "kVp")
        self.kvp_field = Entry(self.contents, width=10)
        self.ma_label = XcistLabel(self.contents, "mA")
        self.ma_field = Entry(self.contents, width=10)
        self.cu_label = XcistLabel(self.contents, "Cu (mm)")
        self.cu_field = Entry(self.contents, width=10)
        self.kvp_label.grid(sticky='W', column=0, row=1, pady=4)
        self.kvp_field.grid(column=1, row=1, padx=10, pady=4)
        self.ma_label.grid(sticky='W', column=2, row=1, pady=4)
        self.ma_field.grid(column=3, row=1, padx=10, pady=4)
        self.cu_label.grid(sticky='W', column=4, row=1, pady=4)
        self.cu_field.grid(column=5, row=1, padx=(10,0), pady=4)

        self.rp_label = XcistLabel(self.contents, "Rotation period (s)")
        self.rp_field = Entry(self.contents, width=10)
        self.bowtie_label = XcistLabel(self.contents, "Bowtie")
        self.bowtie_field = Entry(self.contents, width=10)
        self.rp_label.grid(sticky='W', column=0, columnspan=2, row=2, pady=4)
        self.rp_field.grid(sticky='WE', column=2, columnspan=2, row=2, padx=10, pady=4)
        self.bowtie_label.grid(sticky='W', column=4, row=2, pady=4)
        self.bowtie_field.grid(column=5, row=2, padx=(10,0), pady=4)

        self.ah_label = XcistLabel(self.contents, "Axial or helical")
        self.ah_field = Entry(self.contents, width=10)
        self.ma_mod_label = XcistLabel(self.contents, "mA-mod")
        self.ma_mod_field = Entry(self.contents, width=10)
        self.ah_label.grid(sticky='W', column=0, columnspan=2, row=3, pady=4)
        self.ah_field.grid(sticky='WE', column=2, columnspan=2, row=3, padx=10, pady=4)
        self.ma_mod_label.grid(sticky='W', column=4, row=3, pady=4)
        self.ma_mod_field.grid(column=5, row=3, padx=(10,0), pady=4)

        self.noise_lbl = XcistLabel(self.contents, "Noise")
        self.noise_box = Checkbutton(self.contents, bg=BG_PRIMARY, activebackground=BG_PRIMARY)
        self.mono_lbl = XcistLabel(self.contents, "Mono")
        self.mono_box = Checkbutton(self.contents, bg=BG_PRIMARY, activebackground=BG_PRIMARY)
        self.noise_lbl.grid(sticky='W', column=0, row=4, pady=4)
        self.noise_box.grid(sticky='W', column=1, row=4, padx=(10,0), pady=4)
        self.mono_lbl.grid(sticky='W', column=2, row=4, pady=4)
        self.mono_box.grid(sticky='W', column=3, row=4, padx=(10,0), pady=4)

        self.scan_btn = Button(self.contents, width=10, text="SCAN", bg='red', fg=FG_PRIMARY, font=('Arial 12 bold'))
        self.scan_btn.grid(sticky='SE', column=6, row=4, rowspan=2, padx=(20,0), pady=(20,0))

        self.image = XcistImage(self.frame.frame, "SINOGRAM")
        self.image.frame.grid(column=1, row=0, rowspan=2, padx=(20,0))

    def load_protocol(self):
        filename = filedialog.askopenfilename(title="Open a protocol file", filetypes=[("Configuration Source Files", "*.cfg")])
        self.protocol_file = filename
        self.filename.configure(text="\"" + filename[filename.rindex('/')+1:] + "\"")
        cfg = CFG(filename)
        self.ma_field.delete(0, END)
        self.ma_field.insert(0, cfg.protocol.mA)
        self.rp_field.delete(0, END)
        self.rp_field.insert(0, cfg.protocol.rotationTime)
        self.bowtie_field.delete(0, END)
        self.bowtie_field.insert(0, cfg.protocol.bowtie)
        self.ah_field.delete(0, END)
        self.ah_field.insert(0, cfg.protocol.scanTrajectory)
