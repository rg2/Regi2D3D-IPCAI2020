/*
 * MIT License
 *
 * Copyright (c) 2020 Robert Grupp
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "xregProgOptUtils.h"
#include "xregHDF5.h"
#include "xregFCSVUtils.h"
#include "xregITKIOUtils.h"
#include "xregAnatCoordFrames.h"

#include "IPCAICommon.h"

int main(int argc, char* argv[])
{
  using namespace xreg;
  
  constexpr int kEXIT_VAL_SUCCESS   = 0;
  constexpr int kEXIT_VAL_BAD_USE   = 1;
  constexpr int kEXIT_VAL_BAD_INPUT = 2;

  // First, set up the program options

  ProgOpts po;

  xregPROG_OPTS_SET_COMPILE_DATE(po);

  po.set_help("Extract 3D data (intensity volume, volume segmentation, 3D landmarks) "
              "for a specific specimen from the full-resolution IPCAI 2020 HDF5 data file.");
  
  po.set_arg_usage("<Input H5 file> <Specimen ID> <Output Intensity Volume> "
                   "<Output Segmentation Volume> <Output Landmarks>");
  po.set_min_num_pos_args(5);
  po.set_help_epilogue(fmt::format("\nIPCAI version: {}", IPCAIVersionStr()));

  try
  {
    po.parse(argc, argv);
  }
  catch (const ProgOpts::Exception& e)
  {
    std::cerr << "Error parsing command line arguments: " << e.what() << std::endl;
    po.print_usage(std::cerr);
    return kEXIT_VAL_BAD_USE;
  }

  if (po.help_set())
  {
    po.print_usage(std::cout);
    po.print_help(std::cout);
    return kEXIT_VAL_SUCCESS;
  }

  std::ostream& vout = po.vout();
  
  const std::string src_h5_path     = po.pos_args()[0];
  const std::string spec_id_str     = po.pos_args()[1];
  const std::string dst_intens_path = po.pos_args()[2];
  const std::string dst_seg_path    = po.pos_args()[3];
  const std::string dst_lands_path  = po.pos_args()[4];

  vout << "opening source H5 for reading: " << src_h5_path << std::endl;

  H5::H5File h5(src_h5_path, H5F_ACC_RDONLY);
 
  if (ObjectInGroupH5(spec_id_str, h5))
  {
    H5::Group spec_g = h5.openGroup(spec_id_str);

    {
      vout << "reading intensity volume..." << std::endl;
      const auto vol = ReadITKImageH5Float3D(spec_g.openGroup("vol"));
      
      vout << "writing intensity volume to disk..." << std::endl;
      WriteITKImageToDisk(vol.GetPointer(), dst_intens_path);
    }

    {
      vout << "reading segmentation volume..." << std::endl;
      const auto seg = ReadITKImageH5UChar3D(spec_g.openGroup("vol-seg/image"));
      
      vout << "writing segmentation to disk..." << std::endl;
      WriteITKImageToDisk(seg.GetPointer(), dst_seg_path);
    }

    {
      vout << "reading 3D landmarks..." << std::endl;
      auto lands = ReadLandmarksMapH5Pt3(spec_g.openGroup("vol-landmarks"));

      vout << "converting landmarks to RAS for output to FCSV..." << std::endl;
      ConvertRASToLPS(&lands);

      vout << "writing landmarks to disk..." << std::endl;
      WriteFCSVFileFromNamePtMap(dst_lands_path, lands);
    }
  }
  else
  {
    std::cerr << "ERROR: specimen ID not in source file!" << std::endl;
    return kEXIT_VAL_BAD_INPUT;
  }

  vout << "exiting..." << std::endl;
  
  return kEXIT_VAL_SUCCESS;
}
