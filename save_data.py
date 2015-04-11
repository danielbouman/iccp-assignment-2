"""
Save determined physical quantities or other data to an external file
data_array  : data to export
name        : filename
header      : optional header string (default: blank line)
write_mode  : specify write modus e.g. w = overwrite, a = append (default: a)
"""
## Import libraries
import re   # string editing tools
import six
import numpy as np
def save(data,name,header="",write_mode="a",optional_data=""):
    if isinstance(data, six.string_types):
        write_data = str(optional_data)+" "+data+"\n"
    # elif type(data_array).__module__ == np.__name__:
    #     write_data = "\n".join(str(x) for x in data_array)
    else:
        write_data = str(optional_data)+" "+str(data)
    with open(name+".dat", write_mode) as file:                         # write to file
        file.write(header+"\n")
        file.write(write_data)
    file.close()
    return