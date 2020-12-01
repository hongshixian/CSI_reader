import matplotlib.pyplot as plt
from wifilib import *

path = r"./run_lh_1.dat"

bf = read_bf_file(path)
csi_list = list(map(get_scale_csi,bf))
csi_np = (np.array(csi_list))
csi_amp = np.abs(csi_np)

print("csi shape: ",csi_np.shape)
fig = plt.figure()
plt.plot(csi_amp[:,0,0,3])
plt.show()