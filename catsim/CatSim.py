# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

import numpy as np
from catsim.CommonTools import *
from catsim.OneScan import one_scan
import catsim.PrepView as PrepView

class CatSim:
    '''
    Run airscan, offset scan, phantom scan --> save rawdata (.air .offset .scan)
    Prep view and save .prep data
    Mingye Wu, GE Research
    
    '''
    def __init__(self, *para):
        cfg = CFG(*para)
        self.attrList = ['protocol', 'scanner', 'phantom', 'physics', 'recon', 'resultsName']
        self.cfg_to_self(cfg)
        
    def run_all(self, *para):
        cfg = self.get_current_cfg(*para)
        scanTypes = cfg.protocol.scanTypes
        if scanTypes[0]:
            cfg = self.air_scan(cfg)
        if scanTypes[1]:
            cfg = self.offset_scan(cfg)
        if scanTypes[2]:
            cfg = self.phantom_scan(cfg)
        if all(scanTypes):
            cfg = self.prep_view(cfg)
        return cfg

    def air_scan(self, *para):
        cfg = self.get_current_cfg(*para)
        print('Airscan')
        cfg.sim.thisScanType = [1, 0, 0]
        cfg = one_scan(cfg)
        return cfg
    
    def offset_scan(self, *para):
        cfg = self.get_current_cfg(*para)
        print('Offset scan')
        cfg.sim.thisScanType = [0, 1, 0]
        cfg = one_scan(cfg)
        return cfg

    def phantom_scan(self, *para):
        cfg = self.get_current_cfg(*para)
        print('Phantom scan')
        cfg.sim.thisScanType = [0, 0, 1]
        cfg = one_scan(cfg)
        return cfg
    
    def prep_view(self, *para):
        cfg = self.get_current_cfg(*para)
        print('Prep view')
        cfg = PrepView.prep_view(cfg)
        return cfg
    
    def get_current_cfg(self, *para):
        if len(para)>0:
            cfg = para[0]
        else:
            cfg = self.self_to_cfg()
        return cfg
        
    def load_cfg(self, *para):
        cfg = self.self_to_cfg()
        for cfgFile in para:
            cfg = source_cfg(cfgFile, cfg)
        self.cfg_to_self(cfg)
        return cfg
    
    def self_to_cfg(self):
        cfg = self.cfg
        for attr in self.attrList:
            setattr(cfg, attr, getattr(self, attr))
        return cfg
        
    def cfg_to_self(self, cfg):
        self.cfg = cfg
        for attr in self.attrList:
            setattr(self, attr, getattr(cfg, attr))
            