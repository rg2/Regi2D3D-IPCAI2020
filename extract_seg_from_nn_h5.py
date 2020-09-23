# MIT License
#
# Copyright (c) 2020 Robert Grupp
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# This program extracts a segmentation, for a single projection, inferred by a neural network.
# See https://github.com/rg2/DeepFluoroLabeling-IPCAI2020/tree/master/train_test_code
# This allows the segmentation to be used for the intraoperative registration strategy.
# Use a lossless output format, such as PNG.
# Make sure that the dimensions of this segmentation match the projection dimensions of
# the registration strategy where this segmentation is to be used. e.g. check the amount of
# cropping and downsampling applied.

import sys
import os.path

import h5py as h5

from PIL import Image

if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.stderr.write('USAGE: {} <NN H5 segmentation file> <Proj. Index> <Output Image Path>\n'.format(os.path.basename(sys.argv[0])))
        sys.stderr.flush()
        sys.exit(1)

    src_nn_h5_path = sys.argv[1]
    proj_idx       = int(sys.argv[2])
    dst_seg_path   = sys.argv[3]

    f = h5.File(src_nn_h5_path)
    
    seg = f['nn-segs'][proj_idx,:,:]
    
    Image.fromarray(seg).save(dst_seg_path)

