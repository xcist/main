from tkinter import *
from widgets import *
from colors import *

window = Tk()
window.title("XCIST 2.0")
window.configure(bg=BG_SECONDARY, padx=20)

title_lbl = Label(window, text="XCIST 2.0", bg=BG_SECONDARY, fg=TITLE_COLOR, font=('Arial 40 bold italic'))
title_lbl.grid(column=0, columnspan=3, row=0, pady=20)
credit_lbl = Label(window, text="The development of XCIST was supported by ITCR and the National Cancer Institute of the National Institutes of Health "
                   + "under Award Number U01CA231860 and Award Number U24CAXXXXXX.", bg=BG_SECONDARY, fg=TITLE_COLOR, font=('Arial 10'), wraplength=700)
credit_lbl.grid(column=0, columnspan=3, row=3, pady=20)

# DEFINITION OF 4 MAIN FRAMES
experiment_frame = XcistFrame(window, "EXPERIMENT", 300, 250)
experiment_frame.frame.grid(column=0, row=1, padx=5, pady=5)

protocol_frame = XcistFrame(window, "SCAN PROTOCOL", 700, 250)
protocol_frame.frame.grid(column=1, columnspan=2, row=1, padx=5, pady=5)

dose_frame = XcistFrame(window, "DOSE", 500, 250)
dose_frame.frame.grid(column=0, columnspan=2, row=2, padx=5, pady=5)

recon_frame = XcistFrame(window, "RECONSTRUCTION", 500, 250)
recon_frame.frame.grid(column=2, row=2, padx=5, pady=5)

# EXPERIMENT FRAME CONTENTS
experiment_contents = Frame(experiment_frame.frame, bg=BG_PRIMARY)
experiment_contents.grid(sticky='WE', column=0, row=1)
experiment_contents.grid_columnconfigure(1, weight=1)

name_lbl = XcistLabel(experiment_contents, "Name")
name_field = Entry(experiment_contents)
name_lbl.grid(sticky='W', column=0, row=0, pady=4)
name_field.grid(sticky='WE', column=1, row=0, padx=8, pady=4)

phantom_lbl = XcistLabel(experiment_contents, "Phantom")
phantom_name = XcistNameLabel(experiment_contents, "\"example\"")
phantom_btn = Button(experiment_contents, text="EDIT")
phantom_lbl.grid(sticky='W', column=0, row=1, pady=4)
phantom_name.grid(sticky='W', column=1, row=1, padx=4, pady=4)
phantom_btn.grid(sticky='E', column=2, row=1, pady=4)

scanner_lbl = XcistLabel(experiment_contents, "Scanner")
scanner_name = XcistNameLabel(experiment_contents, "\"example\"")
scanner_btn = Button(experiment_contents, text="EDIT")
scanner_lbl.grid(sticky='W', column=0, row=2, pady=4)
scanner_name.grid(sticky='W', column=1, row=2, padx=4, pady=4)
scanner_btn.grid(sticky='E', column=2, row=2, pady=4)

physics_lbl = XcistLabel(experiment_contents, "Physics")
physics_name = XcistNameLabel(experiment_contents, "\"example\"")
physics_btn = Button(experiment_contents, text="EDIT")
physics_lbl.grid(sticky='W', column=0, row=3, pady=4)
physics_name.grid(sticky='W', column=1, row=3, padx=4, pady=4)
physics_btn.grid(sticky='E', column=2, row=3, pady=4)

# SCAN PROTOCOL FRAME CONTENTS
protocol_contents = Frame(protocol_frame.frame, bg=BG_PRIMARY)
protocol_contents.grid(sticky='NSWE', column=0, row=1)
protocol_contents.grid_columnconfigure(0, weight=1)
protocol_contents.grid_columnconfigure(2, weight=1)
protocol_contents.grid_columnconfigure(4, weight=1)
protocol_contents.grid_rowconfigure(5, weight=1)

protocol_filename_lbl = XcistLabel(protocol_contents, text="Filename")
protocol_filename = XcistNameLabel(protocol_contents, "\"example\"")
protocol_filename_lbl.grid(sticky='W', column=0, columnspan=2, row=0, pady=4)
protocol_filename.grid(sticky='W', column=2, columnspan=4, row=0, pady=4)

kvp_label = XcistLabel(protocol_contents, "kVp")
kvp_field = Entry(protocol_contents, width=10)
ma_label = XcistLabel(protocol_contents, "mA")
ma_field = Entry(protocol_contents, width=10)
cu_label = XcistLabel(protocol_contents, "Cu (mm)")
cu_field = Entry(protocol_contents, width=10)
kvp_label.grid(sticky='W', column=0, row=1, pady=4)
kvp_field.grid(column=1, row=1, padx=10, pady=4)
ma_label.grid(sticky='W', column=2, row=1, pady=4)
ma_field.grid(column=3, row=1, padx=10, pady=4)
cu_label.grid(sticky='W', column=4, row=1, pady=4)
cu_field.grid(column=5, row=1, padx=(10,0), pady=4)

rt_label = XcistLabel(protocol_contents, "Rotation period (s)")
rt_field = Entry(protocol_contents, width=10)
bowtie_label = XcistLabel(protocol_contents, "Bowtie")
bowtie_field = Entry(protocol_contents, width=10)
rt_label.grid(sticky='W', column=0, columnspan=3, row=2, pady=4)
rt_field.grid(column=3, row=2, padx=10, pady=4)
bowtie_label.grid(sticky='W', column=4, row=2, pady=4)
bowtie_field.grid(column=5, row=2, padx=(10,0), pady=4)

ah_label = XcistLabel(protocol_contents, "Axial or helical")
ah_field = Entry(protocol_contents, width=10)
ma_mod_label = XcistLabel(protocol_contents, "mA-mod")
ma_mod_field = Entry(protocol_contents, width=10)
ah_label.grid(sticky='W', column=0, columnspan=3, row=3, pady=4)
ah_field.grid(column=3, row=3, padx=10, pady=4)
ma_mod_label.grid(sticky='W', column=4, row=3, pady=4)
ma_mod_field.grid(column=5, row=3, padx=(10,0), pady=4)

noise_lbl = XcistLabel(protocol_contents, "Noise")
noise_box = Checkbutton(protocol_contents, bg=BG_PRIMARY, activebackground=BG_PRIMARY)
mono_lbl = XcistLabel(protocol_contents, "Mono")
mono_box = Checkbutton(protocol_contents, bg=BG_PRIMARY, activebackground=BG_PRIMARY)
noise_lbl.grid(sticky='W', column=0, row=4, pady=4)
noise_box.grid(sticky='W', column=1, row=4, padx=(10,0), pady=4)
mono_lbl.grid(sticky='W', column=2, row=4, pady=4)
mono_box.grid(sticky='W', column=3, row=4, padx=(10,0), pady=4)

scan_btn = Button(protocol_contents, width=10, text="SCAN", bg='red', fg=FG_PRIMARY, font=('Arial 12 bold'))
scan_btn.grid(sticky='SE', column=6, row=4, rowspan=2, padx=(20,0), pady=(20,0))

protocol_image = XcistImage(protocol_frame.frame, "SINOGRAM")
protocol_image.frame.grid(column=1, row=0, rowspan=2, padx=(20,0))

# DOSE FRAME CONTENTS
dose_contents = Frame(dose_frame.frame, bg=BG_PRIMARY)
dose_contents.grid(sticky='NSWE', column=0, row=1)
dose_contents.grid_columnconfigure(4, weight=1)

dose_filename_lbl = XcistLabel(dose_contents, text="Filename")
dose_filename = XcistNameLabel(dose_contents, "\"example\"")
dose_filename_lbl.grid(sticky='W', column=0, row=0, pady=4)
dose_filename.grid(sticky='W', column=1, columnspan=3, row=0, padx=10, pady=4)

dose_fov_lbl = XcistLabel(dose_contents, text="FOV (mm)")
dose_fov_field = Entry(dose_contents, width=8)
dose_fov_lbl.grid(sticky='W', column=0, row=1, pady=4)
dose_fov_field.grid(column=1, row=1, padx=(10,5), pady=4)

dose_size_lbl = XcistLabel(dose_contents, text="Image size")
dose_size_field1 = Entry(dose_contents, width=8)
dose_size_field2 = Entry(dose_contents, width=8)
dose_size_field3 = Entry(dose_contents, width=8)
dose_size_lbl.grid(sticky='W', column=0, row=2, pady=4)
dose_size_field1.grid(column=1, row=2, padx=(10,5), pady=4)
dose_size_field2.grid(column=2, row=2, padx=5, pady=4)
dose_size_field3.grid(column=3, row=2, padx=(5,0), pady=4)

dose_center_lbl = XcistLabel(dose_contents, text="Center (mm)")
dose_center_field1 = Entry(dose_contents, width=8)
dose_center_field2 = Entry(dose_contents, width=8)
dose_center_field3 = Entry(dose_contents, width=8)
dose_center_lbl.grid(sticky='W', column=0, row=3, pady=4)
dose_center_field1.grid(column=1, row=3, padx=(10,5), pady=4)
dose_center_field2.grid(column=2, row=3, padx=5, pady=4)
dose_center_field3.grid(column=3, row=3, padx=(5,0), pady=4)

dose_btn = Button(dose_contents, width=10, text="DOSE", bg='red', fg=FG_PRIMARY, font=('Arial 12 bold'))
dose_btn.grid(sticky='SE', column=2, columnspan=3, row=4, padx=(20,0), pady=(20,0))

dose_image = XcistImage(dose_frame.frame, "DOSE MAP")
dose_image.frame.grid(column=1, row=0, rowspan=2, padx=(20,0))

# RECONSTRUCTION FRAME CONTENTS
recon_contents = Frame(recon_frame.frame, bg=BG_PRIMARY)
recon_contents.grid(sticky='NSWE', column=0, row=1)
recon_contents.grid_columnconfigure(4, weight=1)
recon_contents.grid_rowconfigure(5, weight=1)

recon_filename_lbl = XcistLabel(recon_contents, text="Filename")
recon_filename = XcistNameLabel(recon_contents, "\"example\"")
recon_filename_lbl.grid(sticky='W', column=0, row=0, pady=4)
recon_filename.grid(sticky='W', column=1, columnspan=3, row=0, padx=10, pady=4)

recon_fov_lbl = XcistLabel(recon_contents, text="FOV (mm)")
recon_fov_field = Entry(recon_contents, width=8)
recon_fov_lbl.grid(sticky='W', column=0, row=1, pady=4)
recon_fov_field.grid(column=1, row=1, padx=(10,5), pady=4)

recon_size_lbl = XcistLabel(recon_contents, text="Image size")
recon_size_field1 = Entry(recon_contents, width=8)
recon_size_field2 = Entry(recon_contents, width=8)
recon_size_field3 = Entry(recon_contents, width=8)
recon_size_lbl.grid(sticky='W', column=0, row=2, pady=4)
recon_size_field1.grid(column=1, row=2, padx=(10,5), pady=4)
recon_size_field2.grid(column=2, row=2, padx=5, pady=4)
recon_size_field3.grid(column=3, row=2, padx=(5,0), pady=4)

recon_center_lbl = XcistLabel(recon_contents, text="Center (mm)")
recon_center_field1 = Entry(recon_contents, width=8)
recon_center_field2 = Entry(recon_contents, width=8)
recon_center_field3 = Entry(recon_contents, width=8)
recon_center_lbl.grid(sticky='W', column=0, row=3, pady=4)
recon_center_field1.grid(column=1, row=3, padx=(10,5), pady=4)
recon_center_field2.grid(column=2, row=3, padx=5, pady=4)
recon_center_field3.grid(column=3, row=3, padx=(5,0), pady=4)

kernel_lbl = XcistLabel(recon_contents, text="Kernel")
kernel_field = Entry(recon_contents, width=12)
kernel_lbl.grid(sticky='W', column=0, row=4, pady=4)
kernel_field.grid(sticky='W', column=1, columnspan=2, row=4, padx=(10,5))

recon_btn = Button(recon_contents, width=10, text="RECON", bg='red', fg=FG_PRIMARY, font=('Arial 12 bold'))
recon_btn.grid(sticky='SE', column=2, columnspan=3, row=4, rowspan=2, padx=(20,0), pady=(20,0))

recon_image = XcistImage(recon_frame.frame, "IMAGE")
recon_image.frame.grid(column=1, row=0, rowspan=2, padx=(20,0))


window.mainloop()
