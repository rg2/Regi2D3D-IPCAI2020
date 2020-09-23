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
#include "xregH5ProjDataIO.h"
#include "xregStringUtils.h"

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

  po.set_help("Extract a single projection data file from the full-resolution "
              "IPCAI 2020 HDF5 data file.");
  
  po.set_arg_usage("<Input H5 file> <Specimen ID> <Projection Index> "
                   "<Output File>");
  po.set_min_num_pos_args(4);
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
  
  const std::string src_h5_path   = po.pos_args()[0];
  const std::string spec_id_str   = po.pos_args()[1];
  const size_type proj_idx        = StringCast<size_type>(po.pos_args()[2]);
  const std::string dst_proj_path = po.pos_args()[3];

  vout << "opening source H5 for reading: " << src_h5_path << std::endl;

  H5::H5File h5(src_h5_path, H5F_ACC_RDONLY);
 
  if (ObjectInGroupH5("proj-params", h5))
  {
    ProjDataF32 pd;
    { 
      vout << "setting up camera..." << std::endl;

      H5::Group proj_params_g = h5.openGroup("proj-params");

      pd.cam.setup(ReadMatrixH5CoordScalar("intrinsic", proj_params_g),
                   ReadMatrixH5CoordScalar("extrinsic", proj_params_g),
                   ReadSingleScalarH5ULong("num-rows", proj_params_g),
                   ReadSingleScalarH5ULong("num-cols", proj_params_g),
                   ReadSingleScalarH5CoordScalar("pixel-row-spacing", proj_params_g),
                   ReadSingleScalarH5CoordScalar("pixel-col-spacing", proj_params_g));
    }

    if (ObjectInGroupH5(spec_id_str, h5))
    {
      H5::Group projs_g = h5.openGroup(fmt::format("{}/projections", spec_id_str));

      const std::string proj_idx_str = fmt::format("{:03d}", proj_idx);

      if (ObjectInGroupH5(proj_idx_str, projs_g))
      {
        H5::Group proj_g = projs_g.openGroup(proj_idx_str);
        
        vout << "reading pixels..." << std::endl;

        pd.img = ReadITKImageH5Float2D(proj_g.openGroup("image"));

        vout << "setting rot-up field..." << std::endl;

        pd.rot_to_pat_up = ReadSingleScalarH5Bool("rot-180-for-up", proj_g) ?
                                  ProjDataRotToPatUp::kONE_EIGHTY : ProjDataRotToPatUp::kZERO;
        
        vout << "writing new projection file to disk..." << std::endl;
        WriteProjDataH5ToDisk(pd, dst_proj_path);
      }
      else
      {
        std::cerr << "ERROR: projection group not in source file!" << std::endl;
        return kEXIT_VAL_BAD_INPUT;
      }
    }
    else
    {
      std::cerr << "ERROR: specimen ID not in source file!" << std::endl;
      return kEXIT_VAL_BAD_INPUT;
    }
  }
  else
  {
    std::cerr << "ERROR: proj-params group not in source file!" << std::endl;
    return kEXIT_VAL_BAD_INPUT;
  }

  vout << "exiting..." << std::endl;
  
  return kEXIT_VAL_SUCCESS;
}
