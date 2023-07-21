from tkinter import *
from widgets import *
from colors import *


class XcistDoseFrame:
    def __init__(self, window):
        self.frame = XcistFrame(window, "DOSE", 500, 250)
        self.frame.frame.grid(column=0, columnspan=2, row=2, padx=5, pady=5)

        self.contents = Frame(self.frame.frame, bg=BG_PRIMARY)
        self.contents.grid(sticky='NSWE', column=0, row=1)
        self.contents.grid_columnconfigure(4, weight=1)

        self.filename_lbl = XcistLabel(self.contents, text="Filename")
        self.filename = XcistNameLabel(self.contents, "\"example\"")
        self.filename_lbl.grid(sticky='W', column=0, row=0, pady=4)
        self.filename.grid(sticky='W', column=1, columnspan=3, row=0, padx=10, pady=4)

        self.fov_lbl = XcistLabel(self.contents, text="FOV (mm)")
        self.fov_field = Entry(self.contents, width=8)
        self.fov_lbl.grid(sticky='W', column=0, row=1, pady=4)
        self.fov_field.grid(column=1, row=1, padx=(10,5), pady=4)

        self.size_lbl = XcistLabel(self.contents, text="Image size")
        self.size_field1 = Entry(self.contents, width=8)
        self.size_field2 = Entry(self.contents, width=8)
        self.size_field3 = Entry(self.contents, width=8)
        self.size_lbl.grid(sticky='W', column=0, row=2, pady=4)
        self.size_field1.grid(column=1, row=2, padx=(10,5), pady=4)
        self.size_field2.grid(column=2, row=2, padx=5, pady=4)
        self.size_field3.grid(column=3, row=2, padx=(5,0), pady=4)

        self.center_lbl = XcistLabel(self.contents, text="Center (mm)")
        self.center_field1 = Entry(self.contents, width=8)
        self.center_field2 = Entry(self.contents, width=8)
        self.center_field3 = Entry(self.contents, width=8)
        self.center_lbl.grid(sticky='W', column=0, row=3, pady=4)
        self.center_field1.grid(column=1, row=3, padx=(10,5), pady=4)
        self.center_field2.grid(column=2, row=3, padx=5, pady=4)
        self.center_field3.grid(column=3, row=3, padx=(5,0), pady=4)

        self.btn = Button(self.contents, width=10, text="DOSE", bg='red', fg=FG_PRIMARY, font=('Arial 12 bold'))
        self.btn.grid(sticky='SE', column=2, columnspan=3, row=4, padx=(20,0), pady=(20,0))

        self.image = XcistImage(self.frame.frame, "DOSE MAP")
        self.image.frame.grid(column=1, row=0, rowspan=2, padx=(20,0))
