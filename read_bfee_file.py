import numpy as np
from Bfee import Bfee
from get_scale_csi import get_scale_csi

if __name__ == '__main__':
    bfee = Bfee.from_file("csi.dat", model_name_encode="gb2312")
    for i in range(len(bfee.all_csi)):
        csi = get_scale_csi(bfee.dicts[i])
        print(csi[:,:,i])

