import torch
import tifffile
import glob
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
from I3Net import I3Net
def main():
    head_dir = r'I:\Lattice\20231224_SIM\128data\Training'
    Raw_path = head_dir + r'\ValRaw'
    model_dir = head_dir + '\\' + 'logfile'
    all_imgs_path = glob.glob(Raw_path + '\\'+'*.tif')
    num = len(all_imgs_path)
    Endstr ='.tif'
    oszie = 512
    channel = 7
    save_path = head_dir + '\\'+'output'
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    device = 'cuda'
    method = 'final_network.pth'

    model_label = r'I3NetHyb'
    model = I3Net()
    model.load_state_dict(torch.load(model_dir + '\\' + model_label + '\\'+ method))
    save_path = head_dir + '\\'+'output'+'\\'+model_label
    if not os.path.exists(save_path):
        os.mkdir(save_path)

    pbar = tqdm(total=num, desc="Processing images")
    for i in range(num):
        numstr = np.array2string(np.array(i+1))
        I = tifffile.imread(Raw_path+ '\\' +numstr+Endstr)
        I = np.array(I)/I.max()  # 将PIL Image转化为numpy数组
        I = I.reshape(1,channel,oszie,oszie)
        im = torch.tensor(I, dtype=torch.float32)
        imnose = np.transpose(im.data.numpy(),(0,2,3,1))
        imnose = imnose[0,...]
        model = model.to(device)
        with torch.no_grad():
            model.eval()
            output = model(im.to(device))
            output = output.to('cpu')
            output = output.detach()
            imde = np.transpose(output.data.numpy(),(0,2,3,1))
            imdeHybrid = imde[0,...]
            imde = imdeHybrid
            imde = (2**16-1)*(imde-imde.min())/(imde.max()-imde.min())
            imde = imde.astype(np.uint16)
            tifffile.imwrite(save_path+ '\\' +numstr+model_label+Endstr,imde)
        pbar.update(1)
    pbar.close()

if __name__ == "__main__":
    main()