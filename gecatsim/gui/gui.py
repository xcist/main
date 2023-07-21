from tkinter import *
from colors import *
from experiment import *
from scan import *
from dose import *
from recon import *
from threading import Thread
import gecatsim as xc

window = Tk()
window.title("XCIST 2.0")
window.configure(bg=BG_SECONDARY, padx=20)

title_lbl = Label(window, text="XCIST 2.0", bg=BG_SECONDARY, fg=TITLE_COLOR, font=('Arial 40 bold italic'))
title_lbl.grid(column=0, columnspan=3, row=0, pady=20)
credit_lbl = Label(window, text="The development of XCIST was supported by ITCR and the National Cancer Institute of the National Institutes of Health "
                   + "under Award Number U01CA231860 and Award Number U24CAXXXXXX.", bg=BG_SECONDARY, fg=TITLE_COLOR, font=('Arial 10'), wraplength=700)
credit_lbl.grid(column=0, columnspan=3, row=3, pady=20)

# DEFINITION OF 4 MAIN GUI FRAMES
experiment_frame = XcistExperimentFrame(window)
scan_frame = XcistScanFrame(window)
dose_frame = XcistDoseFrame(window)
recon_frame = XcistReconFrame(window)

# MAIN FUNCTIONS

def scan():
    try:
        ct = xc.CatSim(experiment_frame.phantom_file, experiment_frame.scanner_file, scan_frame.protocol_file)
        ct.resultsName = "test"
        ct.protocol.viewsPerRotation = 50
        ct.protocol.viewCount = ct.protocol.viewsPerRotation
        ct.protocol.stopViewId = ct.protocol.viewCount-1
        ct.run_all()
        prep = xc.rawread(ct.resultsName+'.prep', [ct.protocol.viewCount, ct.scanner.detectorRowCount, ct.scanner.detectorColCount], 'float')
        scan_frame.image.set_image(prep[:,7,0::10])
    finally:
        scan_frame.image.progress_bar.stop()
        scan_frame.scan_btn.configure(state='normal')

def async_scan():
    scan_frame.scan_btn.configure(state='disabled')
    scan_frame.image.progress_bar.start(5)
    t = Thread(target=scan)
    t.daemon = True
    t.start()
    
scan_frame.scan_btn.configure(command=async_scan)

window.mainloop()
