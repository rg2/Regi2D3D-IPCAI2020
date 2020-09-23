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

#include <fmt/format.h>

#include "xregFCSVUtils.h"
#include "xregITKIOUtils.h"
#include "xregLandmarkMapUtils.h"
#include "xregAnatCoordFrames.h"
#include "xregH5ProjDataIO.h"
#include "xregOpenCVUtils.h"
#include "xregHDF5.h"
#include "xregProjPreProc.h"

#include "IPCAICommon.h"

int main(int argc, char* argv[])
{
  using namespace xreg;
  
  constexpr int kEXIT_VAL_SUCCESS = 0;
  constexpr int kEXIT_VAL_BAD_USE = 1;

  // First, set up the program options

  ProgOpts po;

  xregPROG_OPTS_SET_COMPILE_DATE(po);

  po.set_help("Projects 3D pelvis/vertebra/sacrum/femur labels and landmarks into 2D. "
              "A single-channel label map is written to disk. When the optional final "
              "positional argument is set, then a multi-channel label map is written in "
              "HDF5 format. This format accounts for the possible overlap due to the "
              "line integral nature of X-ray.");
  
  po.set_arg_usage("<Input vol. seg.> <3D Landmarks> <Proj. Data File> "
                   "<Pelvis Pose> <L. Femur Pose> <R. Femur Pose> "
                   "<Output Single-Chan. 2D Seg.> [<Output Multi-Chan. 2D Seg.>]");
  po.set_min_num_pos_args(7);
  po.set_help_epilogue(fmt::format("\nIPCAI version: {}", IPCAIVersionStr()));

  po.add("pat-up", ProgOpts::kNO_SHORT_FLAG, ProgOpts::kSTORE_TRUE, "pat-up",
         "If necessary, rotate the output projections to have the patient oriented \"up.\" "
         "This also removes the rot-pat-up field from the output projections.")
    << false;

  po.add("preproc", ProgOpts::kNO_SHORT_FLAG, ProgOpts::kSTORE_TRUE, "preproc",
         "Apply standard preprocessing such as cropping and, if the \"intensity\" flag is set, "
         "log intensity conversion.")
    << false;

  po.add("ds-factor", 'd', ProgOpts::kSTORE_DOUBLE, "ds-factor",
         "Downsample the camera (2D image) dimensions by this factor, this is done before "
         "projection. e.g. 0.5 --> shrink ny half, 1.0 --> no downsampling.")
    << 1.0;

  po.add("proj-idx", 'p', ProgOpts::kSTORE_UINT32, "proj-idx",
         "Index of the projection to project to (the projection data file may store several projections)")
    << ProgOpts::uint32(0);

  po.add("intensity", 'i', ProgOpts::kSTORE_STRING, "intensity",
         "Save a copy of the intensity projection data at this path, after cropping, downsampling, etc.")
    << "";

  po.add("no-ras2lps", ProgOpts::kNO_SHORT_FLAG, ProgOpts::kSTORE_TRUE, "no-ras2lps",
         "Do NOT convert RAS to LPS (or LPS to RAS) for the 3D landmarks; "
         "RAS to LPS conversion negates the first and second components.")
    << false;

  po.add_backend_flags();

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

  const bool verbose = po.get("verbose");
  std::ostream& vout = po.vout();
  
  const std::string seg_3d_path             = po.pos_args()[0];
  const std::string fcsv_3d_path            = po.pos_args()[1];
  const std::string proj_data_path          = po.pos_args()[2];
  const std::string pelvis_pose_path        = po.pos_args()[3];
  const std::string left_femur_pose_path    = po.pos_args()[4];
  const std::string right_femur_pose_path   = po.pos_args()[5];
  const std::string single_chan_2d_seg_path = po.pos_args()[6];

  const std::string multi_chan_2d_seg_path = (po.pos_args().size() > 7) ? po.pos_args()[7] : std::string();

  const double ds_factor = po.get("ds-factor");
  
  const bool do_ds = std::abs(ds_factor - 1.0) > 0.001;

  const bool ras2lps = !po.get("no-ras2lps").as_bool();

  const size_type proj_idx = po.get("proj-idx").as_uint32();

  const std::string dst_intens_proj_path = po.get("intensity");
  
  const bool save_intens_proj = !dst_intens_proj_path.empty();

  const bool apply_preproc = po.get("preproc");

  const bool rot_pat_up = po.get("pat-up");

  vout << "reading label map..." << std::endl;
  auto ct_labels = ReadITKImageFromDisk<itk::Image<unsigned char,3>>(seg_3d_path);

  //////////////////////////////////////////////////////////////////////////////
  // Get the landmarks

  vout << "reading 3D landmarks..." << std::endl;
  LandMap3 lands_3d = ReadFCSVFileNamePtMap(fcsv_3d_path);

  if (ras2lps)
  {
    vout << "converting landmarks from RAS -> LPS" << std::endl;
    ConvertRASToLPS(&lands_3d);
  }

  vout << "3D Landmarks:\n";
  PrintLandmarkMap(lands_3d, vout);

  vout << "reading projection data..." << std::endl;
  DeferredProjReader pd_reader(proj_data_path);
  
  xregASSERT(proj_idx < pd_reader.num_projs_on_disk());

  auto pd = pd_reader.proj_data_F32()[proj_idx];

  if (save_intens_proj)
  {
    pd.img = pd_reader.read_proj_F32(proj_idx);
  }

  if (apply_preproc)
  {
    vout << "projection preprocessing..." << std::endl;
    
    ProjPreProc preproc;
    preproc.set_debug_output_stream(vout, verbose);
    
    preproc.input_projs = { pd };

    preproc();

    pd = preproc.output_projs[0];
  }

  if (do_ds)
  {
    vout << "downsampling input projection data..." << std::endl;
    pd = DownsampleProjData(pd, ds_factor);
  }

  const auto& cam = pd.cam;

  vout << "reading pelvis pose..." << std::endl;
  const FrameTransform pelvis_pose = ReadITKAffineTransformFromFile(pelvis_pose_path);

  vout << "reading left femur pose..." << std::endl;
  const FrameTransform left_femur_pose = ReadITKAffineTransformFromFile(left_femur_pose_path);

  vout << "reading right femur pose..." << std::endl;
  const FrameTransform right_femur_pose = ReadITKAffineTransformFromFile(right_femur_pose_path);
 
  auto proj_labels = ProjectHipLabels(ct_labels, cam, pelvis_pose, left_femur_pose, right_femur_pose,
                                      ProgOptsDepthRayCasterFactory(po), verbose, vout);
 
  vout << "creating proj. data. struct for the single-channel segmentation..."  << std::endl;
  ProjDataU8 seg_pd;
  
  {
    const FrameTransform vol_to_cam = pelvis_pose.inverse();

    seg_pd.cam = cam;
    seg_pd.img = std::get<0>(proj_labels);

    for (const auto& lkv : lands_3d)
    {
      seg_pd.landmarks.emplace(lkv.first, cam.phys_pt_to_ind_pt(vol_to_cam * lkv.second).head(2));
    }
  }

  if (rot_pat_up && pd.rot_to_pat_up && (*pd.rot_to_pat_up != ProjDataRotToPatUp::kZERO))
  {
    xregASSERT(*pd.rot_to_pat_up == ProjDataRotToPatUp::kONE_EIGHTY);
    
    vout << "rotating projections for patient up..." << std::endl;

    ModifyForPatUp(seg_pd.img.GetPointer(), *pd.rot_to_pat_up);

    for (auto& lkv : seg_pd.landmarks)
    {
      lkv.second[0] = static_cast<CoordScalar>(cam.num_det_cols - 1) - lkv.second[0];
      lkv.second[1] = static_cast<CoordScalar>(cam.num_det_rows - 1) - lkv.second[1];
    }

    if (save_intens_proj)
    {
      ModifyForPatUp(pd.img.GetPointer(), *pd.rot_to_pat_up);
      pd.rot_to_pat_up.reset();
    }
  }
  else
  {
    seg_pd.rot_to_pat_up = pd.rot_to_pat_up;
  }

  if (save_intens_proj)
  {
    pd.landmarks = seg_pd.landmarks;
  
    vout << "writing processed intensity proj. data to disk..." << std::endl;
    WriteProjDataH5ToDisk(pd, dst_intens_proj_path);
  }

  vout << "writing single-channel segmentation to disk..." << std::endl;
  WriteProjDataH5ToDisk(seg_pd, single_chan_2d_seg_path);

  if (!multi_chan_2d_seg_path.empty())
  {
    vout << "writing multi-channel segmentation to disk..." << std::endl;

    H5::H5File h5(multi_chan_2d_seg_path, H5F_ACC_TRUNC);

    WriteImgsAsMultiChannelH5(std::get<1>(proj_labels), "seg", &h5);
    
    h5.flush(H5F_SCOPE_GLOBAL);
    h5.close();
  }

  vout << "exiting..." << std::endl;
  
  return kEXIT_VAL_SUCCESS;
}
