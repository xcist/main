from tkinter import *
from tkinter import ttk
from colors import *


def XcistLabel(parent, text):
    return Label(parent, text=text, bg=BG_PRIMARY, fg=FG_PRIMARY, font=('Arial 10 bold'))

def XcistNameLabel(parent, text):
    return Label(parent, text=text, bg=BG_PRIMARY, fg=FG_SECONDARY, font=('Arial 10'))

class XcistFrame:
    def __init__(self, window, title, width, height):
        self.frame = Frame(window, width=width, height=height, padx=20, pady=20, bg=BG_PRIMARY)
        self.frame.grid_columnconfigure(0, weight=1)
        self.frame.grid_propagate(0)
        self.header = Frame(self.frame, bg=BG_PRIMARY)
        self.header.grid(sticky='NWE', column=0, row=0, pady=(0,20))
        self.header.grid_columnconfigure(0, weight=1)
        self.lbl = Label(self.header, text=title, bg=BG_PRIMARY, fg=FG_PRIMARY, font=('Arial 12 bold'))
        self.lbl.grid(sticky='W', column=0, row=0)
        self.load_btn = Button(self.header, text="LOAD")
        self.load_btn.grid(column=1, row=0, padx=5)
        self.save_btn = Button(self.header, text="SAVE")
        self.save_btn.grid(column=2, row=0)

class XcistImage:
    def __init__(self, parent, title):
        self.frame = Frame(parent, bg=BG_PRIMARY)
        self.title = Button(self.frame, text=title)
        self.title.grid(sticky='WE', column=0, row=0)
        self.image_frame = Frame(self.frame, bg=BG_PRIMARY, width=120, height=160)
        self.image_frame.grid(column=0, row=1)
        self.image_frame.pack_propagate(0)
        self.image = Label(self.image_frame, bg='gray')
        self.image.pack(fill=BOTH, expand=True)
        self.progress_bar = ttk.Progressbar(self.frame, length=100, mode='determinate')
        self.progress_bar.grid(column=0, row=2, pady=(10,0))
        self.progress_bar.start()
