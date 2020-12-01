import numpy as np
import math


def read_bf_file(filename, decoder="python"):

    with open(filename, "rb") as f:
        bfee_list = []
        field_len = int.from_bytes(f.read(2), byteorder='big', signed=False)
        while field_len != 0:
            bfee_list.append(f.read(field_len))
            field_len = int.from_bytes(f.read(2), byteorder='big', signed=False)

    dicts = []


    count = 0                    # % Number of records output
    broken_perm = 0              # % Flag marking whether we've encountered a broken CSI yet
    triangle = [0, 1, 3]           # % What perm should sum to for 1,2,3 antennas

    for array in bfee_list:
        #% Read size and code
        code = array[0]

        # there is CSI in field if code == 187，If unhandled code skip (seek over) the record and continue
        if code != 187:
            #% skip all other info
            continue
        else:
            # get beamforming or phy data
            count =count + 1

            timestamp_low = int.from_bytes(array[1:5], byteorder='little', signed=False)
            bfee_count = int.from_bytes(array[5:7], byteorder='little', signed=False)
            Nrx = array[9]
            Ntx = array[10]
            rssi_a = array[11]
            rssi_b = array[12]
            rssi_c = array[13]
            noise = array[14] - 256
            agc = array[15]
            antenna_sel = array[16]
            b_len = int.from_bytes(array[17:19], byteorder='little', signed=False)
            fake_rate_n_flags = int.from_bytes(array[19:21], byteorder='little', signed=False)
            payload = array[21:]  #get payload
            
            calc_len = (30 * (Nrx * Ntx * 8 * 2 + 3) + 6) / 8
            perm = [1,2,3]
            perm[0] = ((antenna_sel) & 0x3)
            perm[1] = ((antenna_sel >> 2) & 0x3)
            perm[2] = ((antenna_sel >> 4) & 0x3)

            
            

            
            #Check that length matches what it should
            if (b_len != calc_len):
                print("MIMOToolbox:read_bfee_new:size","Wrong beamforming matrix size.")

            #Compute CSI from all this crap :
            if decoder=="python":
                csi = parse_csi(payload,Ntx,Nrx)
            else:
                csi = None
                print("decoder name error! Wrong encoder name:",decoder)
                return

            # % matrix does not contain default values
            if sum(perm) != triangle[Nrx-1]:
                print('WARN ONCE: Found CSI (',filename,') with Nrx=', Nrx,' and invalid perm=[',perm,']\n' )
            else:
                csi[:,perm,:] = csi[:,[0,1,2],:]


            # dict,and return
            bfee_dict = {
                'timestamp_low': timestamp_low,
                'bfee_count': bfee_count,
                'Nrx': Nrx,
                'Ntx': Ntx,
                'rssi_a': rssi_a,
                'rssi_b': rssi_b,
                'rssi_c': rssi_c,
                'noise': noise,
                'agc': agc,
                'antenna_sel': antenna_sel,
                'perm': perm,
                'len': b_len,
                'fake_rate_n_flags': fake_rate_n_flags,
                'calc_len': calc_len,
                'csi': csi}


            dicts.append(bfee_dict)

    return dicts


def parse_csi_new(payload,Ntx,Nrx):
    #Compute CSI from all this crap
    csi = np.zeros(shape=(30,Nrx ,Ntx), dtype=np.dtype(np.complex))
    index = 0

    for i in range(30):
        index += 3
        remainder = index % 8
        for j in range(Nrx):
            for k in range(Ntx):
                real_bin = (int.from_bytes(payload[int(index / 8):int(index/8+2)], byteorder='big', signed=True) >> remainder) & 0b11111111
                real = real_bin
                imag_bin = bytes([(payload[int(index / 8+1)] >> remainder) | (payload[int(index/8+2)] << (8-remainder)) & 0b11111111])
                imag = int.from_bytes(imag_bin, byteorder='little', signed=True)
                tmp = np.complex(float(real), float(imag))
                csi[i, j, k] = tmp
                index += 16
    return csi


def parse_csi(payload,Ntx,Nrx):
    #Compute CSI from all this crap
    csi = np.zeros(shape=(Ntx,Nrx,30), dtype=np.dtype(np.complex))
    index = 0

    for i in range(30):
        index += 3
        remainder = index % 8
        for j in range(Nrx):
            for k in range(Ntx):
                start = index // 8
                real_bin = bytes([(payload[start] >> remainder) | (payload[start+1] << (8-remainder)) & 0b11111111])
                real = int.from_bytes(real_bin, byteorder='little', signed=True)
                imag_bin = bytes([(payload[start+1] >> remainder) | (payload[start+2] << (8-remainder)) & 0b11111111])
                imag = int.from_bytes(imag_bin, byteorder='little', signed=True)
                tmp = np.complex(float(real), float(imag))
                csi[k, j, i] = tmp
                index += 16
    return csi


def db(X, U):
    R = 1
    if 'power'.startswith(U):
        assert X >= 0
    else:
        X = math.pow(abs(X), 2) / R

    return (10 * math.log10(X) + 300) - 300


def dbinv(x):
    return math.pow(10, x / 10)


def get_total_rss(csi_st):
    # Careful here: rssis could be zero
    rssi_mag = 0
    if csi_st['rssi_a'] != 0:
        rssi_mag = rssi_mag + dbinv(csi_st['rssi_a'])
    if csi_st['rssi_b'] != 0:
        rssi_mag = rssi_mag + dbinv(csi_st['rssi_b'])
    if csi_st['rssi_c'] != 0:
        rssi_mag = rssi_mag + dbinv(csi_st['rssi_c'])
    return db(rssi_mag, 'power') - 44 - csi_st['agc']

    
def get_scale_csi(csi_st):
    #Pull out csi
    csi = csi_st['csi']
    # print(csi.shape)
    # print(csi)
    #Calculate the scale factor between normalized CSI and RSSI (mW)
    csi_sq = np.multiply(csi, np.conj(csi)).real
    csi_pwr = np.sum(csi_sq, axis=0)
    csi_pwr = csi_pwr.reshape(1, csi_pwr.shape[0], -1)
    rssi_pwr = dbinv(get_total_rss(csi_st))

    scale = rssi_pwr / (csi_pwr / 30)

    if csi_st['noise'] == -127:
        noise_db = -92
    else:
        noise_db = csi_st['noise']
    thermal_noise_pwr = dbinv(noise_db)

    quant_error_pwr = scale * (csi_st['Nrx'] * csi_st['Ntx'])

    total_noise_pwr = thermal_noise_pwr + quant_error_pwr
    ret = csi * np.sqrt(scale / total_noise_pwr)
    if csi_st['Ntx'] == 2:
        ret = ret * math.sqrt(2)
    elif csi_st['Ntx'] == 3:
        ret = ret * math.sqrt(dbinv(4.5))
    return ret
