import cv2, os
import numpy as np


bmpdir = './color_bmp/'
spdir = './superpixel_10000/'

img_names = []
# get img names
for files in os.listdir(bmpdir):
    if files.endswith('.bmp'):
        img_names.append(files.split('.bmp')[0])

for img_name in img_names:
    # read img
    bmp_file = '%s%s.bmp' % (bmpdir, img_name)
    img = cv2.imread(bmp_file)
    img = np.array(img)
    rows, cols, ch = img.shape

    # read dat file
    dat_filename = '%s/%s.dat' % (spdir, img_name)
    dat = np.fromfile(dat_filename, dtype=np.int32)
    # actual superpixel number
    spnum = max(d for d in dat)
    print spnum

    img_i = np.reshape(img, (rows*cols, 3))
    img_o = np.zeros((rows*cols, 3))
    color_space = np.zeros((spnum, 3))

    for i in range(spnum):
        inds = np.where(dat == i)[0]
        mean = np.mean(img_i[inds], axis=0)
        color_space[i] = np.floor(mean)
        img_o[inds] = np.floor(mean)

    np.savetxt(img_name + '_color_space.txt', color_space, fmt='%d')
    img_o = np.reshape(img_o, (rows, cols, 3))
    cv2.imwrite(img_name + '.bmp', img_o)

