import numpy as np

def get_character_arr(getter, arr_len, char_len):
    char_arr1 = np.empty((char_len,arr_len),dtype='c')
    char_arr1 = getter(char_arr1,char_arr1.shape[1]).T
    names = []
    for i in range(arr_len):
        name = ''.join([char_arr1[i,j].decode("utf-8") for j in range(char_len)])
        names.append(name.strip())
    return names