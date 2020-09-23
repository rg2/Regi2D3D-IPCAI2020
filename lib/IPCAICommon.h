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

#ifndef IPCAI_COMMON_H_
#define IPCAI_COMMON_H_

#include "xregMultiObjMultiLevel2D3DRegi.h"
#include "xregProjData.h"
#include "xregProgOptUtils.h"

namespace xreg
{

std::string IPCAIVersionStr();

using RayCasterFactoryFn = std::function<std::shared_ptr<RayCaster>(void)>;

using SimMetricFactoryFn = std::function<std::shared_ptr<ImgSimMetric2D>(void)>;

RayCasterFactoryFn ProgOptsLineIntRayCasterFactory(ProgOpts& po);

RayCasterFactoryFn ProgOptsDepthRayCasterFactory(ProgOpts& po);

SimMetricFactoryFn ProgOptsSimMetricFactory(ProgOpts& po);

void RunGlobalPelvisRegi(MultiLevelMultiObjRegi* ml_mo_regi,
                         const int method_id,
                         itk::Image<float,3>::Pointer& vol_hu,
                         itk::Image<float,3>::Pointer& vol_att,
                         const LandMap3& lands,
                         const ProjDataF32& pd,
                         const bool save_debug,
                         RayCasterFactoryFn ray_cast_fact,
                         SimMetricFactoryFn sim_metric_fact,
                         const size_type de_pop_size,
                         const size_type grid_batch_size,
                         const size_type pso_num_parts,
                         std::ostream& vout);

void RunFemursAndPelvisRegiUsingPrevRegi(MultiLevelMultiObjRegi* ml_mo_regi,
                                         const FrameTransform& prev_pelvis_regi_cam_to_vol,
                                         std::vector<itk::Image<float,3>::Pointer>& vols_hu,
                                         std::vector<itk::Image<float,3>::Pointer>& vols_att,
                                         const LandMap3& lands,
                                         const ProjDataF32& pd,
                                         const bool save_debug,
                                         RayCasterFactoryFn ray_cast_fact,
                                         SimMetricFactoryFn sim_metric_fact,
                                         std::ostream& vout);

std::tuple<itk::Image<unsigned char,2>::Pointer, std::vector<cv::Mat>>
ProjectHipLabels(itk::Image<unsigned char,3>* ct_labels,
                 const CameraModel& cam,
                 const FrameTransform& pelvis_pose,
                 const FrameTransform& left_femur_pose,
                 const FrameTransform& right_femur_pose,
                 RayCasterFactoryFn ray_cast_fact,
                 const bool verbose,
                 std::ostream& vout);

void RunPelvisAndFemursIntraopRegi(MultiLevelMultiObjRegi* ml_mo_regi,
                                   const int method_id,
                                   std::vector<itk::Image<float,3>::Pointer>& vols_hu,
                                   std::vector<itk::Image<float,3>::Pointer>& vols_att,
                                   const LandMap3& lands_3d,
                                   const ProjDataF32& pd,
                                   const cv::Mat* pd_seg,
                                   const bool save_debug,
                                   RayCasterFactoryFn ray_cast_fact,
                                   SimMetricFactoryFn sim_metric_fact,
                                   std::ostream& vout);

}  // xreg

#endif


