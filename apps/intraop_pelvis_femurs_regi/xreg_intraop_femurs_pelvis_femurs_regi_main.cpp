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

#include <opencv2/imgcodecs.hpp>

#include "xregStringUtils.h"
#include "xregFCSVUtils.h"
#include "xregITKIOUtils.h"
#include "xregITKLabelUtils.h"
#include "xregLandmarkMapUtils.h"
#include "xregAnatCoordFrames.h"
#include "xregH5ProjDataIO.h"
#include "xregHUToLinAtt.h"
#include "xregMultiObjMultiLevel2D3DRegiDebug.h"

#include "IPCAICommon.h"

int main(int argc, char* argv[])
{
  using namespace xreg;
  
  constexpr int kEXIT_VAL_SUCCESS = 0;
  constexpr int kEXIT_VAL_BAD_USE = 1;

  // First, set up the program options

  ProgOpts po;

  xregPROG_OPTS_SET_COMPILE_DATE(po);

  po.set_help("Register both femurs and pelvis model to a single fluoroscopic view using one "
              "of three intraoperative strategies. This program outputs the registered poses "
              "(e.g. mapping camera frame to volume) of the pelvis, left femur, and right femur "
              "as well as the relative poses of each femur with respect to the APP. The APP is "
              "defined to have an origin located at the femoral head center, and will therefore "
              "be different for each femur. The first argument to this program is a numeric ID "
              "of the registration method to use. This number corresponds to either method 1, "
              "2, or 3 as defined in the IPCAI/IJCARS paper. After the initial pelvis "
              "registration, each method registers each femur and finally performs a "
              "simultaneous registration of all objects. Methods 2 and 3 use the 2D segmentation "
              "to weight the patch computations at the 8x downsample level. "
              "Method 1: Unconstrained/unregularized pelvis registration using nominal AP view "
              "as initialization. "
              "Method 2: Regularized pelvis registration following a paired point PnP initialization. "
              "Method 3: Simultaneous intensity and landmark-based pelvis registration.");
  
  po.set_arg_usage("<Method ID> <Input CT vol.> <Input vol. seg.> <3D Landmarks> "
                   "<Proj. Data File> <2D Segmentation> "
                   "<Output Pelvis Regi. Pose> "
                   "<Output L. Femur Regi. Pose> <Output R. Femur Regi. Pose> "
                   "<Output L. Femur Rel. Pose> <Output R. Femur Rel. Pose>"
                   "[<output debug data file>]");
  po.set_min_num_pos_args(11);
  po.set_help_epilogue(fmt::format("\nIPCAI version: {}", IPCAIVersionStr()));

  po.add("proj-idx", 'p', ProgOpts::kSTORE_UINT32, "proj-idx",
         "Index of the projection to register (the projection data file may store several projections)")
    << ProgOpts::uint32(0);

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
 
  const int         regi_method_id                 = StringCast<int>(po.pos_args()[0]);
  const std::string ct_path                        = po.pos_args()[1];
  const std::string seg_path                       = po.pos_args()[2];
  const std::string fcsv_3d_path                   = po.pos_args()[3];
  const std::string proj_data_path                 = po.pos_args()[4];
  const std::string proj_seg_path                  = po.pos_args()[5];
  const std::string dst_pelvis_regi_pose_path      = po.pos_args()[6];
  const std::string dst_left_femur_regi_pose_path  = po.pos_args()[7];
  const std::string dst_right_femur_regi_pose_path = po.pos_args()[8];
  const std::string dst_left_femur_rel_pose_path   = po.pos_args()[9];
  const std::string dst_right_femur_rel_pose_path  = po.pos_args()[10];

  const bool save_debug = po.pos_args().size() > 11;

  const std::string dst_debug_path = save_debug ? po.pos_args()[11] : std::string();

  const bool ras2lps = !po.get("no-ras2lps").as_bool();

  const size_type proj_idx = po.get("proj-idx").as_uint32();

  if ((regi_method_id < 1) || (regi_method_id > 3))
  {
    std::cerr << "ERROR: invalid registration method ID: " << regi_method_id << std::endl;
    return kEXIT_VAL_BAD_USE;
  }

  const bool read_proj_seg = proj_seg_path != "-";
  
  if ((regi_method_id > 1) && !read_proj_seg)
  {
    std::cerr << "ERROR: must use 2D segmentation for methods 2 or 3!" << std::endl;
    return kEXIT_VAL_BAD_USE;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Read in input intensity volume
  
  vout << "reading in source intensity volume..." << std::endl;
  auto vol_hu = ReadITKImageFromDisk<RayCaster::Vol>(ct_path);
  
  vout << "converting HU --> Lin. Att." << std::endl;
  auto vol_att = HUToLinAtt(vol_hu.GetPointer());

  vout << "reading label map..." << std::endl;
  auto ct_labels = ReadITKImageFromDisk<itk::Image<unsigned char,3>>(seg_path);

  {
    vout << "remapping both hemi-pelves to have label 1, "
            "keeping femur labels the same, and and masking out the other labels..." << std::endl;      

    std::vector<unsigned char> lut(256, 0);
    lut[1] = 1;  // left hemi-pelvis
    lut[2] = 1;  // right hemi-pelvis
    lut[5] = 5;  // left femur
    lut[6] = 6;  // right femur
    
    ct_labels = RemapITKLabelMap<unsigned char>(ct_labels.GetPointer(), lut);
  }

  const std::vector<unsigned char> obj_labels = { 1, 5, 6 };

  vout << "cropping HU intensity volumes tightly around objects..." << std::endl;

  auto vols_hu = MakeVolListFromVolAndLabels(vol_hu.GetPointer(), ct_labels.GetPointer(),
                                             obj_labels, -1000);
  
  vout << "cropping linear att. intensity volumes tightly around objects..." << std::endl;

  auto vols_att = MakeVolListFromVolAndLabels(vol_att.GetPointer(), ct_labels.GetPointer(),
                                              obj_labels, 0);

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

  ProjDataF32 pd;

  {
    vout << "reading projection data..." << std::endl;
    DeferredProjReader proj_reader(proj_data_path);
    
    const auto proj_metas = proj_reader.proj_data_F32();

    xregASSERT(proj_idx < proj_metas.size());
    
    pd = proj_metas[proj_idx];
    
    pd.img = proj_reader.read_proj_F32(proj_idx);
    
    vout << "2D Landmarks:\n";
    PrintLandmarkMap(pd.landmarks, vout);
  }

  cv::Mat proj_seg;

  if (read_proj_seg)
  {
    vout << "reading 2D projection segmentation..." << std::endl;
    proj_seg = cv::imread(proj_seg_path, cv::IMREAD_GRAYSCALE);
  }

  MultiLevelMultiObjRegi ml_mo_regi;
  ml_mo_regi.set_debug_output_stream(vout, verbose);
  
  RunPelvisAndFemursIntraopRegi(&ml_mo_regi, regi_method_id, vols_hu, vols_att, lands_3d,
                                pd, read_proj_seg ? &proj_seg : nullptr, save_debug,
                                ProgOptsLineIntRayCasterFactory(po),
                                ProgOptsSimMetricFactory(po),
                                vout);

  vout << "writing pelvis regi pose to disk..." << std::endl;
  WriteITKAffineTransform(dst_pelvis_regi_pose_path, ml_mo_regi.cur_cam_to_vols[0]);
  
  vout << "writing left femur regi pose to disk..." << std::endl;
  WriteITKAffineTransform(dst_left_femur_regi_pose_path, ml_mo_regi.cur_cam_to_vols[1]);
  
  vout << "writing right femur regi pose to disk..." << std::endl;
  WriteITKAffineTransform(dst_right_femur_regi_pose_path, ml_mo_regi.cur_cam_to_vols[2]);
 
  {
    vout << "computing rel. pose of left femur and writing to disk..." << std::endl;
    
    const FrameTransform app_to_vol_fhl = AnteriorPelvicPlaneFromLandmarksMap(lands_3d, kAPP_ORIGIN_LEFT_FH);
    
    WriteITKAffineTransform(dst_left_femur_rel_pose_path,
                            app_to_vol_fhl.inverse() * ml_mo_regi.cur_cam_to_vols[0] *
                            ml_mo_regi.cur_cam_to_vols[1].inverse() * app_to_vol_fhl);
  }
  
  {
    vout << "computing rel. pose of right femur and writing to disk..." << std::endl;
    
    const FrameTransform app_to_vol_fhr = AnteriorPelvicPlaneFromLandmarksMap(lands_3d, kAPP_ORIGIN_RIGHT_FH);
    
    WriteITKAffineTransform(dst_right_femur_rel_pose_path,
                            app_to_vol_fhr.inverse() * ml_mo_regi.cur_cam_to_vols[0] *
                            ml_mo_regi.cur_cam_to_vols[2].inverse() * app_to_vol_fhr);
  }

  if (save_debug)
  {
    vout << "writing debug info to disk..." << std::endl;
    WriteMultiLevel2D3DRegiDebugToDisk(*ml_mo_regi.debug_info, dst_debug_path);
  }

  vout << "exiting..." << std::endl;
  
  return kEXIT_VAL_SUCCESS;
}
