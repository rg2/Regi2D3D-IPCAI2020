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

# This program extracts anatomical landmark locations, for a single projection,
# that were inferred by a neural network and saved in CSV format.
# See https://github.com/rg2/DeepFluoroLabeling-IPCAI2020/tree/master/train_test_code
# The output file format is FCSV, which is a 3D Slicer and xReg compatible file format.
# This allows the estimated landmarks to be used for the intraoperative registration strategy.
# A downsampling factor needs to be passed which matches the amount of downsampling applied 
# to the full-resolution projections to create training data. (e.g. 0.125 for 8x downsampling).

import sys
import csv
import os.path

def write_lands_map_to_fcsv(lands, dst_fcsv_path, lps_to_ras=True):
    with open(dst_fcsv_path, 'w') as f:
        f.write('# Markups fiducial file version = 4.6\n')
        f.write('# CoordinateSystem = 0\n')
        f.write('# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n')

        for (land_name, land_pt) in lands.items():
            x = land_pt[0]
            y = land_pt[1]
            z = land_pt[2]

            if lps_to_ras:
                x *= -1
                y *= -1

            f.write(',{:.8f},{:.8f},{:.8f},0,0,0,1,1,1,0,{},,\n'.format(x,y,z,land_name))

        f.flush()

def extract_lands_map_for_proj_from_nn_csv(src_nn_csv_path, pat_idx, proj_idx, ds_factor):
    kLAND_NAMES = [ 'FH-l',   'FH-r',
                    'GSN-l',  'GSN-r',
                    'IOF-l',  'IOF-r',
                    'MOF-l',  'MOF-r',
                    'SPS-l',  'SPS-r',
                    'IPS-l',  'IPS-r',
                    'ASIS-l', 'ASIS-r']

    kCROP_WIDTH_FULL_RES_PIX = 50
    kORIG_PIX_SPACING = 0.194

    lands_map = { }

    with open(src_nn_csv_path) as f:
        csv_reader = csv.DictReader(f)
        
        for csv_row in csv_reader:
            if (int(csv_row['pat']) == pat_idx) and (int(csv_row['proj']) == proj_idx):
                r = int(csv_row['row'])
                c = int(csv_row['col'])

                if (r >= 0) and (c >= 0):
                    y = ((r / ds_factor) + kCROP_WIDTH_FULL_RES_PIX) * kORIG_PIX_SPACING
                    x = ((c / ds_factor) + kCROP_WIDTH_FULL_RES_PIX) * kORIG_PIX_SPACING

                    lands_map[kLAND_NAMES[int(csv_row['land'])]] = (x,y,0.0)

    return lands_map

if __name__ == '__main__':
    if len(sys.argv) < 6:
        sys.stderr.write('USAGE: {} <NN CSV Path> <Pat. Index> <Proj. Index> <DS. Factor> <FCSV Path>\n'.format(os.path.basename(sys.argv[0])))
        sys.stderr.flush()
        sys.exit(1)

    src_nn_csv_path = sys.argv[1]
    pat_idx         = int(sys.argv[2])
    proj_idx        = int(sys.argv[3])
    ds_factor       = float(sys.argv[4])
    dst_fcsv_path   = sys.argv[5]
    
    write_lands_map_to_fcsv(
            extract_lands_map_for_proj_from_nn_csv(src_nn_csv_path, pat_idx,
                                                   proj_idx, ds_factor),
            dst_fcsv_path)

