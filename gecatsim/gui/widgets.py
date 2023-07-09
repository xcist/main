from tkinter import *
from tkinter import ttk
from colors import *
from PIL import Image, ImageTk
from skimage.transform import resize
import matplotlib.pyplot as plt


def XcistLabel(parent, text):
    return Label(parent, text=text, bg=BG_PRIMARY, fg=FG_PRIMARY, font=('Arial 10 bold'))

def XcistNameLabel(parent, text):
    return Label(parent, anchor='w', text=text, bg=BG_PRIMARY, fg=FG_SECONDARY, font=('Arial 10'))

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
        self.title = Button(self.frame, text=title, command=self.expand_image)
        self.title.grid(sticky='WE', column=0, row=0)
        self.canvas = Canvas(self.frame, width=120, height=160, bd=0, highlightthickness=0)
        self.canvas.grid(column=0, row=1)
        self.image_container = self.canvas.create_image(60, 0, anchor="n",image=[])
        self.progress_bar = ttk.Progressbar(self.frame, length=100, mode='indeterminate')
        self.progress_bar.grid(column=0, row=2, pady=(10,0))
    
    def set_image(self, image):
        self.image = image
        image = resize(image, (160, 160*image.shape[1]//image.shape[0]))
        image = image - image.min()
        image = image * 255 / image.max()
        photo_image = ImageTk.PhotoImage(image=Image.fromarray(image))
        self.photo_image = photo_image
        self.canvas.itemconfig(self.image_container, image=self.photo_image)

    def expand_image(self):
        plt.imshow(self.image, cmap='gray')
        plt.show()
