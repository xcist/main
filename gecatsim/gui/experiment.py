from tkinter import *
from tkinter import filedialog
from widgets import *
from colors import *


class XcistExperimentFrame:
    def __init__(self, window):
        self.frame = XcistFrame(window, "EXPERIMENT", 300, 250)
        self.frame.frame.grid(column=0, row=1, padx=5, pady=5)
        self.frame.load_btn.configure(command=self.load_experiment)
        self.frame.save_btn.configure(command=self.save_experiment)

        self.contents = Frame(self.frame.frame, bg=BG_PRIMARY)
        self.contents.grid(sticky='WE', column=0, row=1)
        self.contents.grid_columnconfigure(1, weight=1)

        self.name_lbl = XcistLabel(self.contents, "Name")
        self.name_field = Entry(self.contents)
        self.name_lbl.grid(sticky='W', column=0, row=0, pady=4)
        self.name_field.grid(sticky='WE', column=1, row=0, padx=8, pady=4)

        self.phantom_lbl = XcistLabel(self.contents, "Phantom")
        self.phantom_name = XcistNameLabel(self.contents, "\"example\"")
        self.phantom_btn = Button(self.contents, text="EDIT", command=self.edit_phantom)
        self.phantom_lbl.grid(sticky='W', column=0, row=1, pady=4)
        self.phantom_name.grid(sticky='W', column=1, row=1, padx=(4,0), pady=4)
        self.phantom_btn.grid(sticky='E', column=2, row=1, padx=(4,0), pady=4)

        self.scanner_lbl = XcistLabel(self.contents, "Scanner")
        self.scanner_name = XcistNameLabel(self.contents, "\"example\"")
        self.scanner_btn = Button(self.contents, text="EDIT", command=self.edit_scanner)
        self.scanner_lbl.grid(sticky='W', column=0, row=2, pady=4)
        self.scanner_name.grid(sticky='W', column=1, row=2, padx=(4,0), pady=4)
        self.scanner_btn.grid(sticky='E', column=2, row=2, padx=(4,0), pady=4)

        self.physics_lbl = XcistLabel(self.contents, "Physics")
        self.physics_name = XcistNameLabel(self.contents, "\"example\"")
        self.physics_btn = Button(self.contents, text="EDIT", command=self.edit_physics)
        self.physics_lbl.grid(sticky='W', column=0, row=3, pady=4)
        self.physics_name.grid(sticky='W', column=1, row=3, padx=(4,0), pady=4)
        self.physics_btn.grid(sticky='E', column=2, row=3, padx=(4,0), pady=4)

    def load_experiment(self):
        filename = filedialog.askopenfilename()

    def save_experiment(self):
        pass

    def edit_phantom(self):
        filename = filedialog.askopenfilename(title="Open a phantom file", filetypes=[("Configuration Source Files", "*.cfg")])
        self.phantom_file = filename
        self.phantom_name.configure(text="\"" + filename[filename.rindex('/')+1:] + "\"")

    def edit_scanner(self):
        filename = filedialog.askopenfilename(title="Open a scanner file", filetypes=[("Configuration Source Files", "*.cfg")])
        self.scanner_file = filename
        self.scanner_name.configure(text="\"" + filename[filename.rindex('/')+1:] + "\"")

    def edit_physics(self):
        filename = filedialog.askopenfilename(title="Open a physics file", filetypes=[("Configuration Source Files", "*.cfg")])
        self.physics_file = filename
        self.physics_name.configure(text="\"" + filename[filename.rindex('/')+1:] + "\"")
