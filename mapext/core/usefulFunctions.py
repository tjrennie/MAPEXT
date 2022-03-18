import h5py
import os
import numpy as np

def h5PullDict(obj):
    out = {}
    if isinstance(obj, h5py.Group):
        for key, value in obj.items():
            if isinstance(value, h5py.Group):
                out[key] = h5PullDict(value)
            else:
                out[key] = np.asarray(value)
    return out

def h5PushDict(obj, dic, ow=True):
    if isinstance(obj, h5py.Group):
        for key, item in dic.items():
            if isinstance(item, dict):
                if key not in obj:
                    subgrp = obj.create_group(key)
                h5PushDict(obj[key], item)
            else:
                if (key in obj) and (ow):
                    del obj[key]
                try:
                    if item != None:
                        obj[key] = item
                except:
                    obj[key] = item

def ensure_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def set_lims(data,lower=0.05,sym=True,upper=0.95):
    if sym:
        upper = 1-lower
    mask = np.all([np.isfinite(data),data!=0],axis=0)
    hist,edge = np.histogram(data[mask], bins=np.linspace(np.nanmin(data),
                                                       np.nanmax(data),
                                                       1001))
    cen = 0.5*(edge[1:]+edge[:-1])
    hist = np.cumsum(hist)
    hist = hist/np.max(hist)
    lims = np.interp([lower,upper],hist,cen)
    return lims
