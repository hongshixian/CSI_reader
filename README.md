
# 前言


数据采集工具csi_tool采集数据并保存为后缀.dat的数据文件，在csi_tool中提供一个c语言函数解析此文件。阅读了c语言的解析代码后发现，数据文件的组织方法与计网中数据十分相似，但略有不同。  

# 数据格式

总体上，整个文件仅由n个bfee组成，巧了，数据文件中应当包含有n个采样信息，这个bfee的意义不言而喻，就是和采样一一对应。  

bfee：
![在这里插入图片描述](https://img-blog.csdnimg.cn/20190821092946768.png)

bfee的数据结构如上图所示。  
前两字节是field_len，之后一字节是code，再之后便是可变长度的field。field_len等于code+field的字长。  
当code为187时，表示field中是信道信息；不是187时，表示field中是其他信息。  
我们关心的是信道信息，其他信息不解析，跳过该bfee即可。

field：
![在这里插入图片描述](https://img-blog.csdnimg.cn/20190821093009966.png)

若code等于187，field有如上图数据格式。  
到这里你一定感觉很熟悉了。
field分为头部和有效载荷(payload)两部分。头部有20字节的固定长度，有效载荷是个可变长度，字长为len。  
头部各字段的数据类型和意义如下表：  

![在这里插入图片描述](https://img-blog.csdnimg.cn/20190821093027667.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L1l1YW5EaWFuTw==,size_16,color_FFFFFF,t_70)

可以见得，头部中包含了主要的信道信息。  
而其中最重要的csi矩阵，分为30个subc，保存在有效载荷中。  
分别对应30个子载波。  

subc的结构如下表所示：

![在这里插入图片描述](https://img-blog.csdnimg.cn/20190821093057560.png)

复数的结构：

![在这里插入图片描述](https://img-blog.csdnimg.cn/20190821093114719.png)

每个subc的开始会有3位的非数据部分，因此subc的长度不是字节(8位)的整数倍，这将导致subc这部分的解析需要按比特操作，增加我解析工作的复杂度。

到这里，整个文件的数据结构都清楚了，开始试着用python来解析run-lxx.dat这个文件。
~~(真想交给王福超来写啊zzz)~~ 

# 文件解析

  我依旧会使用面向对象的方式构建类，不过构造方法无力，属性太多，我选择用静态方法添加属性的方式构建对象。


```python
import numpy as np
```


```python
class Bfee:

    def __init__(self):
        pass

    @staticmethod
    def from_file(filename, model_name_encode="shift-JIS"):

        with open(filename, "rb") as f:
            from functools import reduce
            array = bytes(reduce(lambda x, y: x+y, list(f))) # reduce(函数，list)，将list中元素依次累加

        bfee = Bfee()
        
#         vmd.current_index = 0
        bfee.file_len = len(array)
        bfee.dicts = []
        bfee.all_csi = []

#         vmd.timestamp_low0 = int.from_bytes(array[3:7], byteorder='little', signed=False)
        
#         array = array[3:]   
        
        #%% Initialize variables
        #ret = cell(ceil(len/95),1);    # % Holds the return values - 1x1 CSI is 95 bytes big, so this should be upper bound
        cur = 0                       # % Current offset into file
        count = 0                    # % Number of records output
        broken_perm = 0              # % Flag marking whether we've encountered a broken CSI yet
        triangle = [0, 1, 3]           # % What perm should sum to for 1,2,3 antennas
                    
        while cur < (bfee.file_len - 3):
            #% Read size and code
            #% 将文件数据读取到维度为 sizeA 的数组 A 中，并将文件指针定位到最后读取的值之后。fread 按列顺序填充 A。
            bfee.field_len = int.from_bytes(array[cur:cur+2], byteorder='big', signed=False)
            bfee.code = array[cur+2]
            cur = cur+3
            
            # there is CSI in field if code == 187，If unhandled code skip (seek over) the record and continue
            if bfee.code == 187:
                pass
            else:
                #% skip all other info
                cur = cur + bfee.field_len - 1
                continue
            
            # get beamforming or phy data
            if bfee.code == 187:
                count =count + 1
                
                bfee.timestamp_low = int.from_bytes(array[cur:cur+4], byteorder='little', signed=False)
                bfee.bfee_count = int.from_bytes(array[cur+4:cur+6], byteorder='little', signed=False)
                bfee.Nrx = array[cur+8]
                bfee.Ntx = array[cur+9]
                bfee.rssi_a = array[cur+10]
                bfee.rssi_b = array[cur+11]
                bfee.rssi_c = array[cur+12]
                bfee.noise = array[cur+13] - 256
                bfee.agc = array[cur+14]
                bfee.antenna_sel = array[cur+15]
                bfee.len = int.from_bytes(array[cur+16:cur+18], byteorder='little', signed=False)
                bfee.fake_rate_n_flags = int.from_bytes(array[cur+18:cur+20], byteorder='little', signed=False)
                bfee.calc_len = (30 * (bfee.Nrx * bfee.Ntx * 8 * 2 + 3) + 6) / 8
                bfee.csi = np.zeros(shape=(30,bfee.Nrx ,bfee.Ntx), dtype=np.dtype(np.complex))
                bfee.perm = [1,2,3]
                bfee.perm[0] = ((bfee.antenna_sel) & 0x3)
                bfee.perm[1] = ((bfee.antenna_sel >> 2) & 0x3)
                bfee.perm[2] = ((bfee.antenna_sel >> 4) & 0x3)
        

                cur = cur + 20
                
                #get payload
                payload = array[cur:cur+bfee.len]
                cur = cur + bfee.len
                
                index = 0

                #Check that length matches what it should
                if (bfee.len != bfee.calc_len):
                    print("MIMOToolbox:read_bfee_new:size","Wrong beamforming matrix size.")

                #Compute CSI from all this crap :
                #import struct
                for i in range(30):
                    index += 3
                    remainder = index % 8
                    for j in range(bfee.Nrx):
                        for k in range(bfee.Ntx):
                            real_bin = bytes([(payload[int(index / 8)] >> remainder) | (payload[int(index/8+1)] << (8-remainder)) & 0b11111111])
                            real = int.from_bytes(real_bin, byteorder='little', signed=True)
                            imag_bin = bytes([(payload[int(index / 8+1)] >> remainder) | (payload[int(index/8+2)] << (8-remainder)) & 0b11111111])
                            imag = int.from_bytes(imag_bin, byteorder='little', signed=True)
                            tmp = np.complex(float(real), float(imag))
                            bfee.csi[i, j, k] = tmp
                            index += 16
                
                # % matrix does not contain default values
                if sum(bfee.perm) != triangle[bfee.Nrx-1]:
                    print('WARN ONCE: Found CSI (',filename,') with Nrx=', bfee.Nrx,' and invalid perm=[',bfee.perm,']\n' )
                else:
                    temp_csi = np.zeros(bfee.csi.shape, dtype=np.dtype(np.complex))
                    # bfee.csi[:,bfee.perm[0:bfee.Nrx],:] = bfee.csi[:,0:bfee.Nrx,:]
                    for r in range(bfee.Nrx):
                        temp_csi[:,bfee.perm[r],:] = bfee.csi[:,r,:]
                    bfee.csi = temp_csi
                
                # print phy data
#                 print(vmd.file_len,
#                       vmd.field_len,
#                       vmd.code,
#                       vmd.timestamp_low,
#                       vmd.bfee_count,
#                       vmd.Nrx,
#                       vmd.Ntx,
#                       vmd.rssi_a,
#                       vmd.rssi_b,
#                       vmd.rssi_c,
#                       vmd.noise,
#                       vmd.agc,
#                       vmd.antenna_sel,
#                       vmd.len,
#                       vmd.fake_rate_n_flags,
#                       vmd.calc_len,
#                       vmd.perm,
#                       vmd.csi.shape
#                      )

                # 将类属性导出为dict，并返回
                bfee_dict = {}
                bfee_dict['timestamp_low'] = bfee.timestamp_low
                bfee_dict['bfee_count'] = bfee.bfee_count
                bfee_dict['Nrx'] = bfee.Nrx
                bfee_dict['Ntx'] = bfee.Ntx
                bfee_dict['rssi_a'] = bfee.rssi_a
                bfee_dict['rssi_b'] = bfee.rssi_b
                bfee_dict['rssi_c'] = bfee.rssi_c
                bfee_dict['noise'] = bfee.noise
                bfee_dict['agc'] = bfee.agc
                bfee_dict['antenna_sel'] = bfee.antenna_sel
                bfee_dict['perm'] = bfee.perm
                bfee_dict['len'] = bfee.len
                bfee_dict['fake_rate_n_flags'] = bfee.fake_rate_n_flags
                bfee_dict['calc_len'] = bfee.calc_len
                bfee_dict['csi'] = bfee.csi

                bfee.dicts.append(bfee_dict)
                bfee.all_csi.append(bfee.csi)
        
        return bfee

```


```python
if __name__ == '__main__':
    bfee = Bfee.from_file("run-lxx.dat", model_name_encode="gb2312")
    from pprint import pprint
    pprint(len(bfee.dicts))
    pprint(len(bfee.all_csi))
```

    4993
    4993


# 其他和总结

方法的返回两种结果：  
bfee.dicts字段等同于read_bfee_file() 函数的返回的结果，适用于原来的处理步骤。  
bfee.all_csi字段是所有csi矩阵的列表，可以直接转化成numpy数组，用来弥补字典性能低下的问题。  
两个长度一样。


```python
temp = np.array(vmd.all_csi) 
np.savez('run-lxx.npz', temp)
temp.shape
```




    (4993, 30, 3, 2)



保存为npz格式，  
run-lxx.dat大小1.9Mb，run-lxx.npz变成了14.4Mb  
两种文件的数据是一样多的，dat文件中复数的实部虚部用8位的sign int表示，npz文件中用64位的double表示，数据长度是原来的8倍，文件大小也变8倍。  
可见.dat文件占用比较小

正确的matlab解析步骤应该是：  
1.从文件将头部信息和csi矩阵读取到字典,即read_bfee_file()  
2.依次从字典中取出标准化CSI，即get_scale_csi()  
3.将所有csi整合到一起，保存为csv
