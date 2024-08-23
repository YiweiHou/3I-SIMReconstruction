from RawLoad import *
import torch.utils.data as Data
from GeneratePSF import *
from loss import *
import glob
import numpy as np
import matplotlib.pyplot as plt
from ModelTrainFunc import *
from I3Net import I3Net

def main():
    # Loading 128 training and test data
    head_dir = r"I:\Lattice\20231224_SIM\128data\Training"
    Training_GT_path = head_dir + '\\' + r'GT\*.tif'
    Training_Raw_path = head_dir + '\\' + r'Raw\*.tif'
    Validation_GT_path = head_dir + '\\' + r'ValGT\*.tif'
    Validation_Raw_path = head_dir + '\\' + r'ValRaw\*.tif'
    image_size = 128
    y_train = loaddata(Training_GT_path, image_size*2, 2)
    X_train = loaddata(Training_Raw_path, image_size, 7)
    image_size = 512
    y_val = loaddata(Validation_GT_path, image_size*2, 2)
    X_val = loaddata(Validation_Raw_path, image_size, 7)
    train_data = Data.TensorDataset(X_train, y_train)
    val_data = Data.TensorDataset(X_val, y_val)
    del X_train, y_train

    # Define PSF parameter
    wavelength = 525
    pixelsize = 32.5
    NA = 1.49
    PSF = GPSF(wavelength, pixelsize, NA)
    # Define hyperparameter
    batchsize = 16
    savestep = 1000
    epoch_num = 20
    LR = 0.0001
    reg = 0.01
    SSIM_weight = 0.1
    save_path=head_dir + '\\' + r'logfile'
    if not os.path.exists(save_path):
        os.mkdir(save_path)

    model = I3Net()
    save_path =  head_dir + r'\logfile\I3NetHyb'
    loss_mode = 1 #Dual-GT hybrid loss
    model = modeltraining(model, LR, epoch_num, loss_mode, save_path, savestep, train_data, val_data, batchsize,  PSF, reg ,SSIM_weight)

    #model = I3Net()
    #save_path =  head_dir + r'\128data\Training\logfile\I3Net'
    #loss_mode = 0 #Conventional loss
    #model = modeltraining(model, LR, epoch_num, loss_mode, save_path, savestep, train_data, val_data, batchsize,  PSF, reg ,SSIM_weight)

if __name__ == "__main__":
    main()