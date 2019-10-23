#################################
#  filter + plot 4 type :
#  0. txt  to report
#  1. RBER vs FER + eqBCH
#  2. EbN0 vs FER : EsN0
#  3. err  vs exact_FER + disable_cnt
#  4. err  vs throughput : RBER or SB
#################################
# python -m pip install --upgrade pip
# python -m pip install xxxx
import glob, re, os
import numpy as np
import matplotlib.pyplot as plt

from ip1kb_qc64  import ip1kb_qc64
from ip2kb_qc64  import ip2kb_qc64
from ip4kb_qc128 import ip4kb_qc128
from WSH         import wsh

type_func = 1   # 0, 1, 2, 3, 4
 
if type_func==1:
    en_bch = 0
    code_sel = 'IP.*8256_QC64'
    #code_sel = 'IP.*16640_QC64'
    #code_sel = 'IP.*QC128'
elif type_func==2:
    esn0 = 0
    #code_sel = 'N.*QC128'
    #code_sel = 'N.*128_QC128'
    #code_sel = 'N.*256_QC128'
    #code_sel = 'N.*512_QC128'
    #code_sel = 'N.*1024_QC128'
    #code_sel = 'N.*2048_QC128'
    #code_sel = 'N.*4096_QC128'
    code_sel = 'N.*8320_QC128'
elif type_func==3:
    xlmt     = 250
    code_sel = 'IP9280_8256_QC64_DEG4_5363'
    #code_sel = 'IP18560_16640_QC64'
    #xlmt     = 500
    #code_sel = 'IP37760_32896_QC128'
elif type_func==4:
    xlmt     = 150
    gap      = 10
    clkf     = 220
    max_th   = 400
    code_sel = 'IP.*8256_QC64'
    #xlmt     = 500
    #gap      = 20
    #clkf     = 480
    #max_th   = 3000
    #code_sel = 'IP.*QC128'


try:
    code_sel
except NameError:
    code_sel = ''
    
style = ['bo-', 'rs-', 'y>-', 'mx-', 'cp-', 'gD-', 'b^-', 'r*-', 'yH-', 'm+-', 'c<-']
array_rpt = {}
cycle_rpt = {}


def rpt_txt():
    """
    translate c_report to plot_format
    """
    path_dir = "./"
    fout = open("./plot_report.txt", "w")
    fout.write("\tarray_rpt = {}\n")
    fout.write("\n")
    
    flist = []
    for file in glob.glob(path_dir+"**/report*txt", recursive=True):
        flist.append(file)
    for file in glob.glob(path_dir+"**/fpga*txt", recursive=True):
        flist.append(file)
    
    for file in flist:
        fout.write("\t# %s \n\t" % file)
        mth = re.search(r'ip\dkb', file)
        if mth:
            name = "IP"
        else:
            name = "N"
        mth = re.search(r'[report|fpga](\d+_\d+_QC\d+_DEG\d_\d+.*)\.txt', file)
        label = "array_rpt[\'%s%s\'] = {\n\t" %(name, mth.group(1))
        fout.write(label)
        
        mth = re.search(r'[report|fpga](\d+)_(\d+)_QC(\d+)_DEG(\d)_(\d+).*\.txt', file)
        fout.write('\'len_cw\':%s,\n\t'      % mth.group(1))
        fout.write('\'len_info\':%s,\n\t'    % mth.group(2))
        fout.write('\'qc_size\':%s,\n\t'     % mth.group(3))
        fout.write('\'num_deg\':%s,\n\t'     % mth.group(4))
        fout.write('\'version\':\'%s\',\n\t' % mth.group(5))
        
        ii = 0
        while (1<<ii) < int(mth.group(1)):
            ii += 1
        bch = int((int(mth.group(1))-int(mth.group(2)))*1.0/ii)+1
        
        ### manual bch ###
        if re.search(r'16640_QC64_DEG4', file):
            bch = bch+8
        if re.search(r'report9280_8256_QC64_DEG4_5363', file):
            bch = 79
        ##################
               
        with open(file, "r") as fin:
            snr = []
            average_error = []
            error_block = []
            total_block = []
            bch_err = 0
            xbch = []
            ocorr = [0]*1000
            xcorr = [0]*1000
            for line in fin:
                text = line.rstrip('\n')
                mth = re.search(r'----------------------------------------------', text)
                if mth:
                    xbch.append(bch_err)
                    bch_err = 0
                mth = re.search(r'SNR\s+=\s+(-?\d+\.?\d?)', text)
                if mth:
                    snr.append(mth.group(1))
                mth = re.search(r'Average_Error\s+=\s+(-?\d+\.\d+)', text)
                if mth:
                    average_error.append(mth.group(1))
                mth = re.search(r'Error_Block\s+=\s+(\d+)', text)
                if mth:
                    error_block.append(mth.group(1))
                mth = re.search(r'Total_Block\s+=\s+(\d+)', text)
                if mth:
                    total_block.append(mth.group(1))
                mth = re.search(r'Correctable_Error_(\d+)\s+=\s+(\d+)', text)
                if mth:
                    ocorr[int(mth.group(1))] += int(mth.group(2))
                    if(int(mth.group(1))>bch):
                        bch_err += int(mth.group(2))
                mth = re.search(r'Uncorrectable_Error_(\d+)\s+=\s+(\d+)', text)
                if mth:
                    xcorr[int(mth.group(1))] += int(mth.group(2))
                    if(int(mth.group(1))>bch):
                        bch_err += int(mth.group(2))
                    #print '%s, %s' % (text, mth.group(1))
        
        fout.write('\'snr\':[',)
        for ii in snr:
            fout.write(" %s," % ii)
        fout.write(' ],\n\t')
        fout.write('\'average_error\':[')
        for ii in average_error:
            fout.write(" %s," % ii)
        fout.write(' ],\n\t')
        fout.write('\'error_block\':[')
        for ii in error_block:
            fout.write(" %s," % ii,)
        fout.write(' ],\n\t')
        fout.write('\'total_block\':[')
        for ii in total_block:
            fout.write(" %s," % ii)
        fout.write(' ],\n\t')
        fout.write('\'correctable\':[')
        for ii in ocorr:
            fout.write( " %s," % ii)
        fout.write(' ],\n\t')
        fout.write('\'uncorrectable\':[')
        for ii in xcorr:
            fout.write(" %s," % ii)
        fout.write(' ],\n\t')
        fout.write('\'bch_t\':\'%d\',\n\t' % bch)
        fout.write('\'bch_err\':[')
        for ii in xbch:
            fout.write(" %s," % ii)
        fout.write(' ],\n\t')        
        fout.write('}\t\n')
        fout.write('\n')
    fout.close()
    os.system("\"C:\\Program Files (x86)\\Notepad++\\notepad++.exe\" plot_report.txt")

def keycwlen(elem):
    mth = re.search(r'(\d+)_(\d+)_QC(\d+)_DEG(\d)_(\d+)', elem)
    if mth:
        #return int(mth.group(1))
        return int(mth.group(2))*10000+int(mth.group(1))
    return 0

def plot_rber(array_rpt):
    """
    plot multiple code in rber vs fer
    """
    fig, ax1 = plt.subplots()
    ss = 0   
    #
    kkey = []
    for mode in array_rpt.keys():
        mth = re.search(code_sel, mode)
        if mth:
            kkey.append(mode)
    kkey.sort(key=keycwlen)
    
    for mode in kkey:
        xdata = []
        for ii in array_rpt[mode]['average_error']:
            xdata.append(ii/array_rpt[mode]['len_cw'])
        ydata = []
        for ii in range(len(array_rpt[mode]['error_block'])):
            xx = array_rpt[mode]['error_block'][ii];
            yy = array_rpt[mode]['total_block'][ii];
            ydata.append(xx*1.0/yy)
        #legend = mode
        mth = re.search(r'(\d+)_(\d+)_QC(\d+)_DEG(\d)_(\d+)', mode)
        legend = 'N%s_K%s' % (mth.group(1), mth.group(2))
        ax1.loglog(xdata, ydata, style[ss], label=legend)
        ss = (ss+1)%(len(style)) 
    # BCH
    if en_bch==1:
        for mode in kkey:
            xdata = []
            for ii in array_rpt[mode]['average_error']:
                xdata.append(ii/array_rpt[mode]['len_cw'])
            ydata = []
            for ii in range(len(array_rpt[mode]['error_block'])):
                xx = array_rpt[mode]['bch_err'][ii];
                yy = array_rpt[mode]['total_block'][ii];
                ydata.append(xx*1.0/yy)
            legend = 'bch%s' % array_rpt[mode]['bch_t']
            ax1.loglog(xdata, ydata, 'k--', label=legend)
    
    #color = 'tab:blue'
    color = 'black'
    ax1.set_xlabel('Raw Bit Error Rate', fontsize=20)
    ax1.set_ylabel('Frame Error Rate', color=color, fontsize=20)
    ax1.xaxis.set_tick_params(labelsize=20)
    ax1.yaxis.set_tick_params(labelsize=20)
    #ax1.tick_params(axis='y', labelcolor=color)
    ax1.tick_params(which='major', length=7)
    ax1.tick_params(which='minor', length=5)
    plt.xlim(1e-3, 0.3)
    plt.ylim(2e-8, 1.0)
    plt.title("Dv4 QC-LDPC code", fontsize=20)
    #plt.legend(loc='lower right')
    plt.legend(loc='upper left')
    ax1.grid(True, which="both")

    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    fig.set_size_inches(12,8)
    plt.savefig('rber_fer.png', dpi=200)
    plt.show()    


def plot_ebn0(array_rpt):
    """
    plot ebn0 vs FER : or esn0
    """
    fig, ax1 = plt.subplots()
    ss = 0
    #if esn0==1 :
    #    code_sel = ''
    kkey = []
    for mode in array_rpt.keys():
        mth = re.search(code_sel, mode)
        if not code_sel:
            kkey.append(mode)
        elif mth:
            kkey.append(mode)
    kkey.sort(key=keycwlen)
    #print kkey
    
    for mode in kkey:
        xdata = []
        for ii in array_rpt[mode]['snr']:
            if esn0==1:
                ii += (np.log10(array_rpt[mode]['len_info']*2.0/array_rpt[mode]['len_cw'])*10)
            xdata.append(ii)
        ydata = []
        for ii in range(len(array_rpt[mode]['error_block'])):
            xx = array_rpt[mode]['error_block'][ii];
            yy = array_rpt[mode]['total_block'][ii];
            ydata.append(xx*1.0/yy)    
        mth = re.search(r'(\d+)_(\d+)_QC(\d+)_DEG(\d)_(\d+)', mode)
        legend = 'N%s_K%s' % (mth.group(1), mth.group(2))
        ax1.semilogy(xdata, ydata, style[ss], label=legend)
        ss = (ss+1)%(len(style))    

    if esn0==1 :
        for mode in kkey:
            xdata = []
            for ii in array_rpt[mode]['snr']:
                ii += (np.log10(array_rpt[mode]['len_info']*2.0/array_rpt[mode]['len_cw'])*10)
                xdata.append(ii)
            ydata = []
            for ii in array_rpt[mode]['average_error']:
                ydata.append(ii*1.0/array_rpt[mode]['len_cw']) 
            #legend = 'raw ber'
            #ax1.semilogy(xdata, ydata, 'k--', label=legend)  
            ax1.semilogy(xdata, ydata, 'k--')  
        xname = "Signal Noise Ratio EsN0"        
    else :
        xname = "Signal Noise Ratio EbN0"
    #color = 'tab:blue'
    color = 'black'
    ax1.set_xlabel(xname, fontsize=20)
    ax1.set_ylabel('Frame Error Rate', color=color, fontsize=20)
    ax1.xaxis.set_tick_params(labelsize=16)
    ax1.yaxis.set_tick_params(labelsize=16)
    #ax1.tick_params(axis='y', labelcolor=color)
    ax1.tick_params(which='major', length=7)
    ax1.tick_params(which='minor', length=5)
    plt.xlim(-2.0-esn0*6, 10.0)
    plt.ylim(7e-8, 1.0)
    plt.title("Dv4 QC-LDPC code", fontsize=20)
    #plt.legend(loc='lower right')
    plt.legend(loc='lower left')
    ax1.grid(True, which="both")

    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    fig.set_size_inches(8,6.4)
    plt.savefig('snr_fer.png', dpi=200)
    plt.show()                


def rber2esn0(ber):
    """
    input rber, return EsN0
    """
    rber2snr = {0.308086914:-6.021,0.304861328:-5.821,0.299851563:-5.621,0.298564236:-5.532,0.296335937:-5.441,0.296270508:-5.421,0.294378472:-5.332,0.292434152:-5.241,0.292004883:-5.221,0.289886285:-5.132,0.287700893:-5.041,0.287335938:-5.021,0.285297743:-4.932,0.283637277:-4.841,0.282366211:-4.821,0.281825521:-4.771,0.281181424:-4.732,0.279052455:-4.641,0.278442383:-4.621,0.277128906:-4.571,0.276778646:-4.532,0.274003348:-4.441,0.273930664:-4.421,0.272697917:-4.371,0.271605903:-4.332,0.269302455:-4.241,0.268720703:-4.221,0.268113281:-4.171,0.267409722:-4.132,0.264536830:-4.041,0.263815430:-4.021,0.263605469:-3.971,0.262631944:-3.932,0.259833705:-3.841,0.259461914:-3.821,0.258602865:-3.771,0.257558160:-3.732,0.255454241:-3.641,0.254854492:-3.621,0.254050000:-3.579,0.253369792:-3.571,0.252960069:-3.532,0.250676339:-3.441,0.249507813:-3.421,0.248870312:-3.379,0.248601562:-3.371,0.248004340:-3.332,0.244960938:-3.241,0.244819336:-3.221,0.244139844:-3.179,0.243743490:-3.171,0.242670139:-3.132,0.242864583:-3.122,0.240407366:-3.041,0.239910156:-3.021,0.239231250:-2.979,0.238584635:-2.971,0.237766493:-2.932,0.237539062:-2.922,0.235281250:-2.841,0.234460937:-2.821,0.234000781:-2.779,0.233398438:-2.771,0.232892361:-2.732,0.231953993:-2.722,0.229763393:-2.641,0.229877930:-2.621,0.229240234:-2.610,0.228748438:-2.579,0.228003906:-2.571,0.227569444:-2.532,0.227393229:-2.522,0.224716518:-2.441,0.224562500:-2.421,0.224216797:-2.410,0.223450000:-2.379,0.222546875:-2.371,0.222111111:-2.332,0.221883681:-2.322,0.219918527:-2.241,0.219270508:-2.221,0.219166992:-2.210,0.218378125:-2.179,0.218061198:-2.171,0.217103299:-2.132,0.216800347:-2.122,0.214501116:-2.041,0.214382812:-2.030,0.213983398:-2.021,0.213672852:-2.010,0.213100000:-1.979,0.213076823:-1.971,0.211717882:-1.932,0.211359375:-1.922,0.209195312:-1.841,0.209180804:-1.830,0.208712891:-1.821,0.208370117:-1.810,0.207593750:-1.779,0.207160156:-1.771,0.207117839:-1.761,0.206312500:-1.732,0.206213542:-1.722,0.204262370:-1.661,0.203863839:-1.641,0.203643973:-1.630,0.203360352:-1.621,0.203304687:-1.610,0.202329688:-1.579,0.201811198:-1.571,0.201765625:-1.561,0.200919271:-1.532,0.200572049:-1.522,0.199032552:-1.461,0.198435268:-1.441,0.198246652:-1.430,0.197875000:-1.421,0.197375977:-1.410,0.196872869:-1.383,0.196829688:-1.379,0.196623698:-1.371,0.196136068:-1.361,0.195495660:-1.332,0.195255208:-1.322,0.194252841:-1.283,0.193742188:-1.261,0.192960938:-1.241,0.192717634:-1.230,0.192463867:-1.221,0.192161133:-1.210,0.191713778:-1.183,0.191358594:-1.179,0.190983073:-1.171,0.190822266:-1.161,0.190026910:-1.132,0.189718750:-1.122,0.188313210:-1.083,0.188099609:-1.061,0.187532366:-1.041,0.186848214:-1.030,0.186963867:-1.021,0.186549805:-1.010,0.185663352:-0.983,0.185835156:-0.979,0.185630208:-0.971,0.185252344:-0.969,0.185319010:-0.961,0.184361111:-0.922,0.182711648:-0.883,0.182718750:-0.869,0.182610677:-0.861,0.182011161:-0.841,0.181929687:-0.830,0.181431641:-0.821,0.181128906:-0.810,0.180133523:-0.783,0.180304688:-0.779,0.180093750:-0.771,0.179931250:-0.769,0.179841797:-0.761,0.178654514:-0.722,0.177706676:-0.683,0.177189062:-0.669,0.177048828:-0.661,0.176484375:-0.641,0.175953125:-0.630,0.175919922:-0.621,0.175632813:-0.610,0.174745028:-0.583,0.174773438:-0.579,0.174555990:-0.571,0.174038281:-0.569,0.174081380:-0.561,0.173184028:-0.522,0.172543403:-0.512,0.172486506:-0.483,0.171302344:-0.469,0.171342448:-0.461,0.170925223:-0.441,0.170656250:-0.430,0.170363281:-0.421,0.170026910:-0.412,0.169995117:-0.410,0.168945312:-0.383,0.169224219:-0.379,0.168988281:-0.371,0.168554687:-0.369,0.168723958:-0.361,0.167616319:-0.322,0.167317708:-0.312,0.166250710:-0.283,0.166017969:-0.269,0.165861328:-0.261,0.165363839:-0.241,0.165229911:-0.230,0.164463542:-0.212,0.164462891:-0.210,0.163852273:-0.183,0.163651562:-0.179,0.163417969:-0.171,0.163419531:-0.169,0.162870443:-0.161,0.161679688:-0.112,0.160710938:-0.083,0.160347656:-0.069,0.160212240:-0.061,0.159791295:-0.041,0.159330357:-0.030,0.158475694:-0.012,0.158916992:-0.010,0.158005682:0.017,0.157848958:0.029,0.157892188:0.031,0.157509766:0.039,0.156395833:0.088,0.155105824:0.117,0.155527344:0.131,0.154789714:0.139,0.154223214:0.159,0.153915179:0.170,0.153346354:0.188,0.153366211:0.190,0.152102273:0.217,0.152287760:0.229,0.151727344:0.231,0.151917969:0.239,0.150884375:0.280,0.150578125:0.288,0.149953835:0.317,0.149088281:0.331,0.149171875:0.339,0.148383929:0.370,0.147591146:0.388,0.147810547:0.390,0.146970170:0.417,0.146692969:0.431,0.146383464:0.439,0.145397917:0.480,0.145003472:0.488,0.144230114:0.517,0.144052344:0.531,0.143655599:0.539,0.142802455:0.570,0.142409040:0.580,0.142074653:0.588,0.142240234:0.590,0.141441051:0.617,0.141555469:0.631,0.140880859:0.639,0.139861458:0.680,0.139703125:0.688,0.138649858:0.717,0.138512500:0.731,0.138114583:0.739,0.137273438:0.770,0.137098214:0.780,0.136770833:0.788,0.135943182:0.817,0.135590625:0.831,0.135343099:0.839,0.134158854:0.880,0.134000000:0.888,0.133684495:0.902,0.133183949:0.917,0.132639844:0.931,0.132596354:0.939,0.131767857:0.970,0.131734933:0.980,0.131194444:0.988,0.130711914:1.000,0.130471591:1.017,0.130122656:1.031,0.129860026:1.039,0.128843750:1.080,0.128573785:1.088,0.128201172:1.100,0.128157452:1.102,0.127723722:1.117,0.127331250:1.131,0.126040179:1.180,0.125955729:1.188,0.125439941:1.200,0.124993608:1.217,0.124711719:1.231,0.123320313:1.280,0.122785590:1.288,0.122725098:1.300,0.122685096:1.302,0.122265625:1.317,0.121878125:1.331,0.120592634:1.380,0.120415625:1.380,0.120225694:1.388,0.120146973:1.400,0.119558239:1.417,0.119189063:1.431,0.117773437:1.480,0.117528646:1.488,0.117298828:1.500,0.117450721:1.502,0.116461719:1.531,0.115046317:1.580,0.115177604:1.580,0.114988715:1.588,0.114517578:1.600,0.113792188:1.631,0.112551562:1.680,0.112235243:1.688,0.112014160:1.700,0.112052284:1.702,0.111142969:1.731,0.109934710:1.780,0.109686458:1.780,0.109640625:1.788,0.109215332:1.800,0.108476562:1.831,0.107076562:1.880,0.106937500:1.888,0.106675293:1.900,0.106506010:1.902,0.105840625:1.931,0.104503348:1.980,0.104488542:1.980,0.104300347:1.988,0.104069824:2.000,0.103947115:2.002,0.101852604:2.080,0.101717882:2.088,0.101464844:2.100,0.101373197:2.102,0.099173549:2.180,0.099410938:2.180,0.099134549:2.188,0.098879395:2.200,0.098986779:2.202,0.096682813:2.280,0.096563368:2.288,0.096234863:2.300,0.096358173:2.302,0.094154018:2.380,0.094236458:2.380,0.093994792:2.388,0.093641113:2.400,0.093592548:2.402,0.091671875:2.480,0.091220215:2.500,0.091322115:2.502,0.089186384:2.580,0.089166146:2.580,0.088719727:2.600,0.088712740:2.602,0.087424154:2.649,0.086679688:2.680,0.086201172:2.700,0.086167067:2.702,0.084921875:2.749,0.084226563:2.780,0.084208854:2.780,0.083738281:2.800,0.083752404:2.802,0.082848505:2.834,0.082530599:2.849,0.081783333:2.880,0.081309570:2.900,0.081247596:2.902,0.080379076:2.934,0.080002604:2.949,0.079347098:2.980,0.079365625:2.980,0.078897949:3.000,0.078800481:3.002,0.078219460:3.027,0.078093410:3.034,0.077820313:3.049,0.076567308:3.102,0.075756037:3.127,0.075544158:3.134,0.075471680:3.149,0.074627232:3.180,0.074145433:3.202,0.073606889:3.227,0.073417459:3.234,0.073108398:3.249,0.071856971:3.302,0.071311435:3.327,0.070901155:3.334,0.070663086:3.349,0.070010045:3.380,0.069512019:3.402,0.068886719:3.427,0.069086310:3.429,0.068831182:3.434,0.068565625:3.441,0.068410156:3.449,0.067219952:3.502,0.067023330:3.506,0.066560724:3.527,0.066532948:3.534,0.066311914:3.541,0.066159831:3.549,0.066187700:3.551,0.065873806:3.566,0.065030649:3.602,0.064895976:3.606,0.064587358:3.627,0.064315104:3.629,0.064243886:3.634,0.064073047:3.641,0.063946940:3.649,0.064004607:3.651,0.063565584:3.664,0.063628147:3.666,0.063117969:3.688,0.062834736:3.702,0.062730522:3.706,0.062273082:3.727,0.062247283:3.734,0.061986914:3.741,0.061820638:3.749,0.061730970:3.751,0.061456620:3.764,0.061364366:3.766,0.061149071:3.780,0.060986272:3.788,0.060673678:3.802,0.060474101:3.806,0.060223011:3.827,0.059850818:3.829,0.060080503:3.834,0.059794922:3.841,0.059658203:3.849,0.059560296:3.851,0.059329564:3.864,0.059244249:3.866,0.059002745:3.880,0.058827232:3.888,0.058549880:3.902,0.058317744:3.906,0.057821733:3.927,0.057880774:3.934,0.057667578:3.941,0.057553060:3.949,0.057514824:3.951,0.057241571:3.964,0.057172201:3.966,0.056895059:3.980,0.056697321:3.988,0.056456130:4.002,0.056345569:4.006,0.055829545:4.027,0.055903646:4.029,0.055831861:4.034,0.055631055:4.041,0.055468750:4.049,0.055472957:4.051,0.055227796:4.064,0.055031359:4.066,0.054921242:4.080,0.054705022:4.088,0.054231485:4.106,0.053928622:4.127,0.053723505:4.134,0.053609766:4.141,0.053437826:4.149,0.053415064:4.151,0.053219572:4.164,0.053033312:4.166,0.052833826:4.180,0.052574888:4.188,0.052237800:4.206,0.051860795:4.227,0.051951265:4.229,0.051736073:4.234,0.051569141:4.241,0.051404046:4.251,0.051112664:4.264,0.051180339:4.266,0.050816723:4.280,0.050674219:4.288,0.050366438:4.306,0.049850497:4.327,0.049772079:4.334,0.049625195:4.341,0.049471955:4.351,0.049207031:4.364,0.049162869:4.366,0.048926098:4.380,0.048732589:4.388,0.048356914:4.406,0.047924716:4.427,0.047937500:4.429,0.047839334:4.434,0.047659570:4.441,0.047534655:4.451,0.047239309:4.464,0.047226345:4.466,0.046951014:4.480,0.046820536:4.488,0.046410210:4.506,0.046075284:4.527,0.045847266:4.541,0.045643630:4.551,0.045342722:4.564,0.045370768:4.566,0.045064400:4.580,0.044965290:4.588,0.044693279:4.606,0.044219105:4.627,0.044176711:4.629,0.044025000:4.641,0.043802083:4.651,0.043647821:4.664,0.043477105:4.666,0.043249578:4.680,0.043135491:4.688,0.042784461:4.706,0.042412997:4.727,0.042126758:4.741,0.041922075:4.751,0.041677426:4.764,0.041657552:4.766,0.041472551:4.780,0.041338058:4.788,0.041061965:4.806,0.040646307:4.827,0.040595982:4.829,0.040419922:4.841,0.040304888:4.851,0.040003906:4.864,0.039956706:4.866,0.039739020:4.880,0.039600558:4.888,0.039274615:4.906,0.038625781:4.941,0.038484575:4.951,0.038282072:4.964,0.038255317:4.966,0.038009079:4.980,0.037917187:4.988,0.037572132:5.006,0.037183408:5.029,0.037012305:5.041,0.036812099:5.051,0.036612459:5.064,0.036586697:5.066,0.036343117:5.080,0.036213281:5.088,0.035853489:5.106,0.035354883:5.141,0.035217949:5.151,0.035036184:5.164,0.034953342:5.166,0.034778927:5.180,0.034571987:5.188,0.034259418:5.206,0.033934561:5.229,0.033745703:5.241,0.033541667:5.251,0.033342722:5.264,0.033362739:5.266,0.033128801:5.280,0.033020871:5.288,0.032745719:5.306,0.032194727:5.341,0.032010617:5.351,0.031851151:5.364,0.031728624:5.366,0.031588894:5.380,0.031528237:5.388,0.031205908:5.406,0.030679492:5.441,0.030539263:5.451,0.030376439:5.464,0.030314779:5.466,0.030176943:5.480,0.029969643:5.488,0.029674443:5.506,0.029210156:5.541,0.029065304:5.551,0.028867599:5.564,0.028824544:5.566,0.028679265:5.580,0.028480246:5.588,0.028319028:5.606,0.027636418:5.651,0.027467311:5.664,0.027403646:5.666,0.027216427:5.680,0.027134710:5.688,0.026889983:5.706,0.026256410:5.751,0.026113076:5.764,0.026070421:5.766,0.025900549:5.780,0.025767522:5.788,0.025451306:5.806,0.024751439:5.864,0.024712565:5.866,0.024601351:5.880,0.024426004:5.888,0.024161922:5.906,0.023461760:5.964,0.023453885:5.966,0.023266258:5.980,0.023148549:5.988,0.022920163:6.006,0.022216488:6.064,0.022174154:6.066,0.022006757:6.080,0.021879464:6.088,0.021731807:6.106,0.021426827:6.127,0.020964952:6.166,0.020815878:6.180,0.020723103:6.188,0.020538099:6.206,0.020254181:6.227,0.019787543:6.266,0.019668285:6.280,0.019593415:6.288,0.019340860:6.306,0.019090779:6.327,0.018695964:6.366,0.018561867:6.380,0.018461272:6.388,0.018262307:6.406,0.018041703:6.427,0.017615668:6.466,0.017400893:6.488,0.017217894:6.506,0.017029599:6.527,0.016607747:6.566,0.016396763:6.588,0.016200664:6.606,0.016002641:6.627,0.015605469:6.666,0.015411161:6.688,0.014998900:6.727,0.014770267:6.751,0.014658529:6.766,0.014448549:6.788,0.014094850:6.827,0.013885417:6.851,0.013744249:6.866,0.013553460:6.888,0.013218420:6.927,0.013010417:6.951,0.012875760:6.966,0.012670424:6.988,0.012364217:7.027,0.012157609:7.051,0.011856696:7.088,0.011553587:7.127,0.011359375:7.151,0.011072098:7.188,0.010782372:7.227,0.010599864:7.251,0.010321875:7.288,0.009880661:7.351,0.009613616:7.388,0.009148211:7.451,0.008525136:7.551,0.007928895:7.651,0.007308537:7.751,0.006779325:7.851,0.006227129:7.951,0.005770494:8.051,0.005296988:8.151,0.004845335:8.251,0.004448143:8.351,0.004077785:8.451,0.003724638:8.551,0.003396626:8.651,0.003073483:8.751,0.002798687:8.851,0.002538723:8.951,0.002287024:9.051,0.002067822:9.151,0.001858809:9.251,0.001670969:9.351,0.001495471:9.451,0.001336843:9.551}
    prv_ii = 0
    for ii in sorted(rber2snr.keys()):
        if(ber<=ii):
            break
        prv_ii = ii
            
    print("%f:%f, %f:%f"%(prv_ii, rber2snr[prv_ii], ii, rber2snr[ii]))
    if ber>0.308086914:
        print("Error: rber2esn0 input=%f > 0.308086914"%ber)
        ans = -6.020600
    elif ber < 0.001336843:
        print("Error: rber2esn0 input=%f < 0.001336843"%ber)
        ans = 9.551
    else:
        ans = rber2snr[prv_ii] - (rber2snr[prv_ii]-rber2snr[ii])*(ber-prv_ii)/(ii-prv_ii)
    return ans
    
def plot_logy0(ax1,xdata,ydata):
    jj = 0
    legend = 'exact FER'
    for ii in range(len(ydata)):
        if ydata[ii]==0:
            if jj==1:
                ax1.semilogy(xdata[ii-jj:ii], ydata[ii-jj:ii], 'ko', markersize=8, fillstyle='none')
                jj = 0
            elif jj>1:
                ax1.semilogy(xdata[ii-jj:ii], ydata[ii-jj:ii], 'ko-', linewidth=3, markersize=8, fillstyle='none')
                jj = 0
        else:
            jj +=1
    if jj>1:
        ax1.semilogy(xdata[ii-jj+1:ii], ydata[ii-jj+1:ii], 'ko-', linewidth=3, markersize=8, fillstyle='none')
    ax1.set_ylim(2e-8, 1e0)    

def plot_exact(array_rpt):
    """
    plot error vs exact_FER + diable_cnt
    """
    fig, ax1 = plt.subplots()
    kkey = []
    for mode in array_rpt.keys():
        mth = re.search(code_sel, mode)
        if not code_sel:
            kkey.append(mode)
        elif mth:
            kkey.append(mode)
    print(kkey)
    mode = kkey[0] 
    
    xdata = range(xlmt)
    ydata = []
    for ii in range(len(array_rpt[mode]['uncorrectable'])):
        xx = array_rpt[mode]['uncorrectable'][ii];
        yy = array_rpt[mode]['correctable'][ii];
        if ((yy+xx)!=0):
            ydata.append(xx*1.0/(yy+xx))
        else:
            ydata.append(0)
    mth = re.search(r'(\d+)_(\d+)_QC(\d+)_DEG(\d)_(\d+)', mode)
    plot_logy0(ax1, xdata, ydata[0:xlmt])

    # BCH
    xx = int(array_rpt[mode]['bch_t'])
    legend = 'Eq BCH=%s' % xx  
    ax1.semilogy([xx,xx], [1,1e-12], 'r--', label=legend, linewidth=4)
    
    # instantiate a second axes that shares the same x-axis
    xdata = range(xlmt)
    ydata = array_rpt[mode]['uncorrectable']
    
    ax2 = ax1.twinx()  
    color = 'tab:blue'
    ax2.set_ylabel('disable histogram', color=color, fontsize=20)  # we already handled the x-label with ax1
    legend = 'Uncorrectable'
    #ax2.plot(xdata, ydata[0:xlmt], 'bs', color=color, label=legend, markersize=8)
    ax2.semilogy(xdata, ydata[0:xlmt], 'bs-', color=color, label=legend, markersize=6)
    ax2.yaxis.set_tick_params(labelsize=16)    
    #ax2.tick_params(axis='y', labelcolor=color)    
    ax2.tick_params(which='major', labelcolor=color, length=5)
    ax2.tick_params(which='minor', labelcolor=color, length=3)
    ax2.set_ylim(1, 30000)
    
    color = 'black'
    ax1.set_xlabel('exact error bit', fontsize=20)
    ax1.set_ylabel('exact Frame Error Rate', color=color, fontsize=24)
    ax1.xaxis.set_tick_params(labelsize=16)
    ax1.yaxis.set_tick_params(labelsize=16)
    #ax1.tick_params(axis='y', labelcolor=color)
    ax1.tick_params(which='major', length=7)
    ax1.tick_params(which='minor', length=5)    
    plt.xlim(0, xlmt)

    xdata = 0
    for ii in array_rpt[mode]['total_block']:
        xdata += ii 
    xdata = xdata*1.0/1000000*int(int(mth.group(2))/8192)
    tbox = '\n'.join((r'N   =%s b'%(mth.group(1)), 
                      r'K   =%s b'%(mth.group(2)), 
                      r'P   =%s b'%(int(mth.group(1))-int(mth.group(2))),
                      r'ver=%s'%(mth.group(5)),
                      r'$%5.2f$ GB'%(xdata) ))
    plt.text(0.7*xlmt, 2, tbox, fontsize=18)
    
    plt.title("Dv4 QC-LDPC code", fontsize=20)
    #plt.legend(loc='lower right')
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    ax1.grid(True, which="both")
    
    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    fig.set_size_inches(10,8)
    plt.savefig('exact_err.png', dpi=100)
    plt.show()   
    

def plot_cycle(cycle_rpt):
    """
    plot single code in error vs throughput
    """
    fig, ax1 = plt.subplots()

    kkey = []
    for mode in cycle_rpt.keys():
        mth = re.search(code_sel, mode)
        if not code_sel:
            kkey.append(mode)
        elif mth:
            kkey.append(mode)
    ss   = 0
    dnum = int(xlmt/gap)+1
    for mode in kkey:
        ydata = [0]*dnum
        pltnn = [0]*dnum
        for ii in range(len(cycle_rpt[mode]['inerr'])):
            kk = round(cycle_rpt[mode]['inerr'][ii]/gap)
            if kk>=dnum:
                kk = dnum-1
            ydata[kk] += cycle_rpt[mode]['cycle'][ii]
            pltnn[kk] += 1
        for ii in range(len(ydata)):
            if pltnn[ii]==0:
                ydata[ii] = 0
            else:
                ydata[ii] = cycle_rpt[mode]['len_info']*clkf/8.0/(ydata[ii]*1.0/pltnn[ii])
        xdata = range(0, dnum*gap, gap)
        mth = re.search(r'(\d+)_(\d+)_QC(\d+)_DEG(\d)_(\d+)', mode)
        legend = 'HB N%s_K%s' % (mth.group(1), mth.group(2))    
        ax1.plot(xdata, ydata, style[ss][0], label=legend)
        ss = (ss+1)%(len(style))
    
    #xdata = []
    #for ii in cycle_rpt[code_sel]['inerr']:
    #    for jj in range(0,xlmt,gap):
    #        if jj > ii:
    #            ii = jj-int(gap/2.0)
    #            break
    #    xdata.append(ii)
    #ydata = []
    #for ii in cycle_rpt[code_sel]['cycle']:
    #    ii = cycle_rpt[code_sel]['len_info']*clkf/8.0/ii
    #    if ii > max_th:
    #        ii = max_th
    #    ydata.append(ii)
    #mth = re.search(r'N(\d+)_(\d+)_QC(\d+)_DEG(\d)_(\d+)', code_sel)
    #legend = 'HB N%s_K%s' % (mth.group(1), mth.group(2))    
    #ax1.plot(xdata, ydata, 'bo', label=legend, fillstyle='none', markersize=8)
    
    ax1.set_xlabel('exact error bit', fontsize=24)
    ax1.set_ylabel('throughput MBs', fontsize=24)
    ax1.xaxis.set_tick_params(labelsize=24)
    ax1.yaxis.set_tick_params(labelsize=24)
    #ax1.tick_params(axis='y', labelcolor=color)
    ax1.tick_params(which='major', length=7)
    ax1.tick_params(which='minor', length=5)    
    plt.xlim(0, xlmt)
    plt.ylim(0, max_th)
    plt.title("QC-LDPC codec ip clk=%d" % clkf, fontsize=20)
    #plt.legend(loc='lower right')
    plt.legend(loc='upper right', fontsize=24)
    ax1.grid(True, which="both")
    
    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    fig.set_size_inches(12,8)
    plt.savefig('throughput.png', dpi=200)
    plt.show()    
    

if __name__=="__main__":
    ip1kb_qc64.rpt(array_rpt)
    ip2kb_qc64.rpt(array_rpt)
    ip4kb_qc128.rpt(array_rpt)
    wsh.rpt(array_rpt)
    
    ip1kb_qc64.cycle(cycle_rpt)
    ip4kb_qc128.cycle(cycle_rpt)
    
    if type_func==0:
        rpt_txt()
    elif type_func==1:
        plot_rber(array_rpt)
    elif type_func==2:
        plot_ebn0(array_rpt)
    elif type_func==3:
        plot_exact(array_rpt)
    elif type_func==4:
        plot_cycle(cycle_rpt)

    