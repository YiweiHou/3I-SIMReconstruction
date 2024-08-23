import numpy as np
from skimage import io, transform
import torch
import cv2
import glob
from tqdm import tqdm


def loaddata(imgs_path, aimingsize, channel):
    all_imgs_path = glob.glob(imgs_path)
    images = np.zeros((len(all_imgs_path), aimingsize, aimingsize, channel))

    # 使用tqdm来显示进度条
    pbar = tqdm(total=len(all_imgs_path), desc="Loading images")
    for i, img_path in enumerate(all_imgs_path):
        img = io.imread(img_path, cv2.IMREAD_GRAYSCALE)
        img = np.array(img) / img.max()
        images[i] = img.transpose([1, 2, 0])
        pbar.update()

    pbar.close()
    images = np.transpose(images, (0, 3, 2, 1))
    images = torch.tensor(images, dtype=torch.float32)
    return images