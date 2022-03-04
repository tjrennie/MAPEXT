import h5py
import os

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
