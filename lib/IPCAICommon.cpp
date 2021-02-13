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

#include "IPCAICommon.h"

#include <fmt/format.h>

#include "xregProjPreProc.h"
#include "xregMultiObjMultiLevel2D3DRegiDebug.h"
#include "xregSE3OptVars.h"
#include "xregITKBasicImageUtils.h"
#include "xregAnatCoordFrames.h"
#include "xregRayCastProgOpts.h"
#include "xregRayCastInterface.h"
#include "xregImgSimMetric2DPatchCommon.h"
#include "xregImgSimMetric2DGradImgParamInterface.h"
#include "xregImgSimMetric2DProgOpts.h"
#include "xregIntensity2D3DRegiExhaustive.h"
#include "xregIntensity2D3DRegiPSO.h"
#include "xregIntensity2D3DRegiBOBYQA.h"
#include "xregIntensity2D3DRegiDiffEvo.h"
#include "xregRegi2D3DPenaltyFnGlobalPelvis.h"
#include "xregIntensity2D3DRegiCMAES.h"
#include "xregRegi2D3DPenaltyFnSE3EulerDecomp.h"
#include "xregNormDist.h"
#include "xregRegi2D3DPenaltyFnSE3Mag.h"
#include "xregFoldNormDist.h"
#include "xregNullDist.h"
#include "xregPnPUtils.h"
#include "xregRegi2D3DPenaltyFnLandReproj.h"
#include "xregOpenCVUtils.h"
#include "xregProj3DLabelsTo2D.h"
#include "xregITKOpenCVUtils.h"

extern const char* IPCAI_PROJ_VERSION;

namespace // un-named
{

using namespace xreg;
  
using UseCurEstForInit = MultiLevelMultiObjRegi::Level::SingleRegi::InitPosePrevPoseEst;

const size_type kPROJ_CROP_PIX = 50;

std::shared_ptr<ImgSimMetric2D> MakeSimMetric(SimMetricFactoryFn& sim_metric_fact, const CoordScalar ds_factor)
{
  auto sm = sim_metric_fact();
  {
    auto* grad_sm = dynamic_cast<ImgSimMetric2DGradImgParamInterface*>(sm.get());

    grad_sm->set_smooth_img_before_sobel_kernel_radius(5);
  }
 
  {
    auto* patch_sm = dynamic_cast<ImgSimMetric2DPatchCommon*>(sm.get());
    xregASSERT(patch_sm);

    patch_sm->set_patch_radius(std::lround(ds_factor * 41));
    patch_sm->set_patch_stride(1);
  }

  return sm;
}

void SetupGlobalPelvisRegiMethod1(MultiLevelMultiObjRegi* ml_mo_regi,
                                  const LandMap3& lands,
                                  RayCasterFactoryFn ray_cast_fact,
                                  SimMetricFactoryFn sim_metric_fact,
                                  const size_type de_pop_size,
                                  const bool pat_is_up,
                                  const size_type grid_batch_size,
                                  std::ostream& vout)
{
  vout << "setting up global method 2 (DE, Grid, CMAES, BOBYQA)..." << std::endl;
  
  if (ml_mo_regi->debug_info)
  {  
    ml_mo_regi->debug_info->regi_names = { { "DE", "GRID" }, { "CMAES" }, { "BOBYQA" } };
  }
  
  // se(3) lie algebra vector space for optimization
  auto se3_vars = std::make_shared<SE3OptVarsLieAlg>();

  // level 1: DE, GRID
  // level 2: CMA-ES
  // level 3: BOBYQA
  ml_mo_regi->levels.resize(3);
  
  // Level 1
  {
    vout << "  setting up pelvis regi level 1..." << std::endl;

    auto& lvl = ml_mo_regi->levels[0];

    lvl.fixed_imgs_to_use = { 0 };

    lvl.ds_factor = 0.03125;
    
    vout << "    setting up ray caster..." << std::endl;
    lvl.ray_caster = ray_cast_fact();
    
    vout << "    setting up sim metric..." << std::endl;
    lvl.sim_metrics = { MakeSimMetric(sim_metric_fact, lvl.ds_factor) };
    
    // Two regis, DE then Grid
    lvl.regis.resize(2);

    for (auto& regi : lvl.regis)
    {
      regi.mov_vols    = { 0 };  // pelvis vol pose is optimized over
      regi.ref_frames  = { 0 };  // use ref frame that is camera aligned
      regi.static_vols = { };    // No other objects exist

      auto init_guess_fn = std::make_shared<UseCurEstForInit>();
      init_guess_fn->vol_idx = 0;
      regi.init_mov_vol_poses = { init_guess_fn };
    }

    {
      auto de_regi = std::make_shared<Intensity2D3DRegiDiffEvo>();
      de_regi->set_opt_vars(se3_vars);

      de_regi->set_box_contraint({ 60 * kDEG2RAD, 40 * kDEG2RAD, 10 * kDEG2RAD, 200, 200, 250 });
      de_regi->set_max_num_iters(400);
      de_regi->set_pop_size(de_pop_size);

      // setup regularization
      auto reg = std::make_shared<Regi2D3DPenaltyFnGlobalPelvis>();
      reg->vol_idx            = 0;
      reg->pat_is_up          = pat_is_up;
      reg->left_fh_wrt_vol    = lands.find("FH-l")->second;
      reg->right_fh_wrt_vol   = lands.find("FH-r")->second;
      reg->left_asis_wrt_vol  = lands.find("ASIS-l")->second;
      reg->right_asis_wrt_vol = lands.find("ASIS-r")->second;
      reg->left_iof_wrt_vol   = lands.find("IOF-l")->second;
      reg->right_iof_wrt_vol  = lands.find("IOF-r")->second;

      de_regi->set_penalty_fn(reg);
      de_regi->set_img_sim_penalty_coefs(0.9, 0.1);

      lvl.regis[0].regi = de_regi;
    }

    {
      auto ex_regi = std::make_shared<Intensity2D3DRegiExhaustive>();
      ex_regi->set_opt_vars(se3_vars);

      ConstSpacedMeshGrid::Range1DList ranges;

      if (true)
      {
        // Set of ranges to use in production
        ranges = { { -5 * kDEG2RAD, 5 * kDEG2RAD, 1 * kDEG2RAD },
                   { -5 * kDEG2RAD, 5 * kDEG2RAD, 1 * kDEG2RAD },
                   { -1 * kDEG2RAD, 1 * kDEG2RAD, 1 * kDEG2RAD },
                   { -10,           10,           2            },
                   { -10,           10,           2            },
                   { -50,           50,           10           } };
      }
      else
      {
        // Set of ranges for debugging on limited GPU
        ranges = { {   0 * kDEG2RAD,  0 * kDEG2RAD,   0 * kDEG2RAD },
                   { -40 * kDEG2RAD, 40 * kDEG2RAD,   5 * kDEG2RAD },
                   {   0 * kDEG2RAD,  0 * kDEG2RAD,   0 * kDEG2RAD },
                   {   0,             0,              0            },
                   {   0,             0,              0            },
                   {   0,             0,              0            } };
      }
      
      ex_regi->set_cam_wrt_vols(ConstSpacedMeshGrid(ranges), grid_batch_size);
      
      lvl.regis[1].regi = ex_regi;
    }
  }  // Level 1
  
  // Level 2
  {
    vout << "  setting up pelvis regi level 2..." << std::endl;

    auto& lvl = ml_mo_regi->levels[1];

    lvl.fixed_imgs_to_use = { 0 };
    
    lvl.ds_factor = 0.125;

    vout << "    setting up ray caster..." << std::endl;
    lvl.ray_caster = ray_cast_fact();
    
    vout << "    setting up sim metric..." << std::endl;
    lvl.sim_metrics = { MakeSimMetric(sim_metric_fact, lvl.ds_factor) };
    
    // CMA-ES
    lvl.regis.resize(1);

    auto& regi = lvl.regis[0];

    regi.mov_vols    = { 0 };  // pelvis vol pose is optimized over
    regi.ref_frames  = { 0 };  // use ref frame that is camera aligned
    regi.static_vols = { };    // No other objects exist

    auto init_guess_fn = std::make_shared<UseCurEstForInit>();
    init_guess_fn->vol_idx = 0;
    regi.init_mov_vol_poses = { init_guess_fn };

    auto cmaes_regi = std::make_shared<Intensity2D3DRegiCMAES>();
    cmaes_regi->set_opt_vars(se3_vars);
    cmaes_regi->set_opt_x_tol(0.01);
    cmaes_regi->set_opt_obj_fn_tol(0.01);
    cmaes_regi->set_pop_size(100);
    cmaes_regi->set_sigma({ 15 * kDEG2RAD, 15 * kDEG2RAD, 30 * kDEG2RAD, 50, 50, 100 });

    auto pen_fn = std::make_shared<Regi2D3DPenaltyFnSE3EulerDecomp>();
    pen_fn->rot_x_pdf   = std::make_shared<NormalDist1D>(0, 10 * kDEG2RAD);
    pen_fn->rot_y_pdf   = std::make_shared<NormalDist1D>(0, 10 * kDEG2RAD);
    pen_fn->rot_z_pdf   = std::make_shared<NormalDist1D>(0, 10 * kDEG2RAD);
    pen_fn->trans_x_pdf = std::make_shared<NormalDist1D>(0, 20);
    pen_fn->trans_y_pdf = std::make_shared<NormalDist1D>(0, 20);
    pen_fn->trans_z_pdf = std::make_shared<NormalDist1D>(0, 100);

    cmaes_regi->set_penalty_fn(pen_fn);
    cmaes_regi->set_img_sim_penalty_coefs(0.9, 0.1);
    
    regi.regi = cmaes_regi;
  }  // Level 2
  
  // Level 3
  {
    vout << "  setting up pelvis regi level 3..." << std::endl;

    auto& lvl = ml_mo_regi->levels[2];

    lvl.fixed_imgs_to_use = { 0 };
    
    lvl.ds_factor = 0.25;

    vout << "    setting up ray caster..." << std::endl;
    lvl.ray_caster = ray_cast_fact();
    
    vout << "    setting up sim metric..." << std::endl;
    lvl.sim_metrics = { MakeSimMetric(sim_metric_fact, lvl.ds_factor) };
    
    // BOBYQA
    lvl.regis.resize(1);

    auto& regi = lvl.regis[0];

    regi.mov_vols    = { 0 };  // pelvis vol pose is optimized over
    regi.ref_frames  = { 0 };  // use ref frame that is camera aligned
    regi.static_vols = { };    // No other objects exist

    auto init_guess_fn = std::make_shared<UseCurEstForInit>();
    init_guess_fn->vol_idx = 0;
    regi.init_mov_vol_poses = { init_guess_fn };

    auto bobyqa_regi = std::make_shared<Intensity2D3DRegiBOBYQA>();
    bobyqa_regi->set_opt_vars(se3_vars);
    bobyqa_regi->set_opt_x_tol(0.0001);
    bobyqa_regi->set_opt_obj_fn_tol(0.0001);
    bobyqa_regi->set_bounds({ 2.5 * kDEG2RAD, 2.5 * kDEG2RAD, 2.5 * kDEG2RAD,
                              5, 5, 10 });
    
    regi.regi = bobyqa_regi;
  }  // Level 3
}

void SetupGlobalPelvisRegiMethod2(MultiLevelMultiObjRegi* ml_mo_regi,
                                  RayCasterFactoryFn ray_cast_fact,
                                  SimMetricFactoryFn sim_metric_fact,
                                  const size_type grid_batch_size,
                                  const size_type pso_num_parts,
                                  std::ostream& vout)
{
  vout << "setting up global method 2 (Grid, PSO, BOBYQA, BOBYQA)..." << std::endl;
  
  if (ml_mo_regi->debug_info)
  {  
    ml_mo_regi->debug_info->regi_names = { { "GRID", "PSO" }, { "BOBYQA1" }, { "BOBYQA2" } };
  }
  
  // se(3) lie algebra vector space for optimization
  auto se3_vars = std::make_shared<SE3OptVarsLieAlg>();
  
  // level 1: Grid and PSO
  // level 2: BOBYQA
  // level 3: BOBYQA
  
  ml_mo_regi->levels.resize(3);

  // Level 1
  {
    vout << "  setting up pelvis regi level 1..." << std::endl;

    auto& lvl = ml_mo_regi->levels[0];

    lvl.fixed_imgs_to_use = { 0 };

    lvl.ds_factor = 0.03125;
    
    vout << "    setting up ray caster..." << std::endl;
    lvl.ray_caster = ray_cast_fact();
    
    vout << "    setting up sim metric..." << std::endl;
    lvl.sim_metrics = { MakeSimMetric(sim_metric_fact, lvl.ds_factor) };
    
    // Two regis, Grid then PSO
    lvl.regis.resize(2);

    for (auto& regi : lvl.regis)
    {
      regi.mov_vols    = { 0 };  // pelvis vol pose is optimized over
      regi.ref_frames  = { 0 };  // use ref frame that is camera aligned
      regi.static_vols = { };    // No other objects exist

      auto init_guess_fn = std::make_shared<UseCurEstForInit>();
      init_guess_fn->vol_idx = 0;
      regi.init_mov_vol_poses = { init_guess_fn };
    }

    {
      auto ex_regi = std::make_shared<Intensity2D3DRegiExhaustive>();
      ex_regi->set_opt_vars(se3_vars);

      ConstSpacedMeshGrid::Range1DList ranges;

      if (true)
      {
        // Set of ranges to use in production
        ranges = { { -60 * kDEG2RAD, 60 * kDEG2RAD, 7.5 * kDEG2RAD },
                   { -40 * kDEG2RAD, 40 * kDEG2RAD,   5 * kDEG2RAD },
                   {   0 * kDEG2RAD,  0 * kDEG2RAD,   0 * kDEG2RAD },
                   { -200,           200,            20            },
                   { -200,           200,            20            },
                   { -250,           250,            25            } };
      }
      else
      {
        // Set of ranges for debugging on limited GPU
        ranges = { {   0 * kDEG2RAD,  0 * kDEG2RAD,   0 * kDEG2RAD },
                   { -40 * kDEG2RAD, 40 * kDEG2RAD,   5 * kDEG2RAD },
                   {   0 * kDEG2RAD,  0 * kDEG2RAD,   0 * kDEG2RAD },
                   {   0,             0,              0            },
                   {   0,             0,              0            },
                   {   0,             0,              0            } };
      }
      
      ex_regi->set_cam_wrt_vols(ConstSpacedMeshGrid(ranges), grid_batch_size);
      
      lvl.regis[0].regi = ex_regi;
    }

    {
      auto pso_regi = std::make_shared<Intensity2D3DRegiPSO>();
      pso_regi->set_opt_vars(se3_vars);
      
      pso_regi->set_num_particles(pso_num_parts);
      pso_regi->set_max_num_iters(50);

      pso_regi->set_box_contraint({ 7.5 * kDEG2RAD, 10 * kDEG2RAD, 10 * kDEG2RAD,
                                    20, 20, 25 });
      
      lvl.regis[1].regi = pso_regi;
    }
  }  // Level 1

  // Level 2
  {
    vout << "  setting up pelvis regi level 2..." << std::endl;

    auto& lvl = ml_mo_regi->levels[1];

    lvl.fixed_imgs_to_use = { 0 };
    
    lvl.ds_factor = 0.125;

    vout << "    setting up ray caster..." << std::endl;
    lvl.ray_caster = ray_cast_fact();
    
    vout << "    setting up sim metric..." << std::endl;
    lvl.sim_metrics = { MakeSimMetric(sim_metric_fact, lvl.ds_factor) };
    
    // BOBYQA
    lvl.regis.resize(1);

    auto& regi = lvl.regis[0];

    regi.mov_vols    = { 0 };  // pelvis vol pose is optimized over
    regi.ref_frames  = { 0 };  // use ref frame that is camera aligned
    regi.static_vols = { };    // No other objects exist

    auto init_guess_fn = std::make_shared<UseCurEstForInit>();
    init_guess_fn->vol_idx = 0;
    regi.init_mov_vol_poses = { init_guess_fn };

    auto bobyqa_regi = std::make_shared<Intensity2D3DRegiBOBYQA>();
    bobyqa_regi->set_opt_vars(se3_vars);
    bobyqa_regi->set_opt_x_tol(0.0001);
    bobyqa_regi->set_opt_obj_fn_tol(0.0001);
    bobyqa_regi->set_bounds({ 5 * kDEG2RAD, 5 * kDEG2RAD, 5 * kDEG2RAD,
                              10, 10, 20 });
    
    regi.regi = bobyqa_regi;
  }  // Level 2
  
  // Level 3
  {
    vout << "  setting up pelvis regi level 3..." << std::endl;

    auto& lvl = ml_mo_regi->levels[2];

    lvl.fixed_imgs_to_use = { 0 };
    
    lvl.ds_factor = 0.25;

    vout << "    setting up ray caster..." << std::endl;
    lvl.ray_caster = ray_cast_fact();
    
    vout << "    setting up sim metric..." << std::endl;
    lvl.sim_metrics = { MakeSimMetric(sim_metric_fact, lvl.ds_factor) };
    
    // BOBYQA
    lvl.regis.resize(1);

    auto& regi = lvl.regis[0];

    regi.mov_vols    = { 0 };  // pelvis vol pose is optimized over
    regi.ref_frames  = { 0 };  // use ref frame that is camera aligned
    regi.static_vols = { };    // No other objects exist

    auto init_guess_fn = std::make_shared<UseCurEstForInit>();
    init_guess_fn->vol_idx = 0;
    regi.init_mov_vol_poses = { init_guess_fn };

    auto bobyqa_regi = std::make_shared<Intensity2D3DRegiBOBYQA>();
    bobyqa_regi->set_opt_vars(se3_vars);
    bobyqa_regi->set_opt_x_tol(0.0001);
    bobyqa_regi->set_opt_obj_fn_tol(0.0001);
    bobyqa_regi->set_bounds({ 2.5 * kDEG2RAD, 2.5 * kDEG2RAD, 2.5 * kDEG2RAD,
                              5, 5, 10 });
    
    regi.regi = bobyqa_regi;
  }  // Level 3
}

void AddSetupForFemursAndPelvisRegi(MultiLevelMultiObjRegi* ml_mo_regi,
                                    const LandMap3& lands,
                                    RayCasterFactoryFn ray_cast_fact,
                                    SimMetricFactoryFn sim_metric_fact,
                                    std::ostream& vout)
{
  vout << "adding setup for femurs and pelvis regi (CMAES L. Femur, CMAES R. Femur, All)" << std::endl;

  // Assumes that Pelvis is volume 0, L. Femur is volume 1, and R. Femur is volume 2

  const Pt3 fh_l_wrt_vol = lands.find("FH-l")->second;
  const Pt3 fh_r_wrt_vol = lands.find("FH-r")->second;

  const FrameTransform app_to_vol_std = AnteriorPelvicPlaneFromLandmarksMap(lands);
 
  auto get_vol_to_app = [&app_to_vol_std] (const Pt3& origin_wrt_vol)
  {
    FrameTransform app_to_vol = app_to_vol_std;
    app_to_vol.matrix().block(0,3,3,1) = origin_wrt_vol;
    
    return FrameTransform(app_to_vol.inverse());
  };

  const FrameTransform vol_pelvis_to_app = get_vol_to_app((fh_l_wrt_vol + fh_r_wrt_vol) / 2);
  const FrameTransform vol_fhl_to_app    = get_vol_to_app(fh_l_wrt_vol);
  const FrameTransform vol_fhr_to_app    = get_vol_to_app(fh_r_wrt_vol);
  
  const size_type pelvis_app_ref_frame_idx = ml_mo_regi->ref_frames.size();
  const size_type fh_l_app_ref_frame_idx   = pelvis_app_ref_frame_idx + 1;
  const size_type fh_r_app_ref_frame_idx   = fh_l_app_ref_frame_idx + 1;

  ml_mo_regi->ref_frames.insert(ml_mo_regi->ref_frames.end(),
                                { MultiLevelMultiObjRegi::MakeStaticRefFrame(vol_pelvis_to_app, true),
                                  MultiLevelMultiObjRegi::MakeStaticRefFrame(vol_fhl_to_app, true),
                                  MultiLevelMultiObjRegi::MakeStaticRefFrame(vol_fhr_to_app, true) });
  
  
  // se(3) lie algebra vector space for optimization
  auto se3_vars = std::make_shared<SE3OptVarsLieAlg>();

  // so(3) lie algebra vector space for optimization
  auto so3_vars = std::make_shared<SO3OptVarsLieAlg>();

  if (ml_mo_regi->debug_info)
  {
    ml_mo_regi->debug_info->regi_names.insert(ml_mo_regi->debug_info->regi_names.end(),
                                              { { "CMAES-L", "CMAES-R" }, { "BOBYQA" } });
  }

  const size_type prev_num_levels = ml_mo_regi->levels.size();
  
  // new level 1: CMA-ES Left Femur, CMA-ES Right Femur
  // level 2: BOBYQA Simultaneous Pelvis, Both Femurs
  ml_mo_regi->levels.resize(prev_num_levels + 2);

  // New Level 1
  {
    vout << "  setting up femurs regi new level 1..." << std::endl;

    auto& lvl = ml_mo_regi->levels[prev_num_levels];

    lvl.fixed_imgs_to_use = { 0 };

    lvl.ds_factor = 0.125;
    
    vout << "    setting up ray caster..." << std::endl;
    lvl.ray_caster = ray_cast_fact();
    
    vout << "    setting up sim metric..." << std::endl;
    lvl.sim_metrics = { MakeSimMetric(sim_metric_fact, lvl.ds_factor) };
    
    // Two regis, L. Femur and R. Femur
    lvl.regis.resize(2);
    
    // left femur CMA-ES
    {
      auto& regi = lvl.regis[0];

      regi.mov_vols   = { 1 };  // left femur vol pose is optimized over
      regi.ref_frames = { fh_l_app_ref_frame_idx };  // use left femur APP ref frame
    }
    
    // right femur CMA-ES
    {
      auto& regi = lvl.regis[1];

      regi.mov_vols   = { 2 };  // right femur vol pose is optimized over
      regi.ref_frames = { fh_r_app_ref_frame_idx };  // use right femur APP ref frame
    }

    // common parameters for each femur
    for (auto& regi : lvl.regis)
    {
      regi.static_vols = { 0 };  // keep pelvis vol fixed in projections

      auto static_pelvis_fn = std::make_shared<UseCurEstForInit>();
      static_pelvis_fn->vol_idx = 0;
      regi.static_vol_poses = { static_pelvis_fn };

      // also use the pelvis pose as the initial estimate of the femur
      regi.init_mov_vol_poses = { static_pelvis_fn };
      
      auto cmaes_regi = std::make_shared<Intensity2D3DRegiCMAES>();
      cmaes_regi->set_opt_vars(so3_vars);
      cmaes_regi->set_opt_x_tol(0.01);
      cmaes_regi->set_opt_obj_fn_tol(0.01);
      cmaes_regi->set_pop_size(100);
      cmaes_regi->set_sigma({ 30 * kDEG2RAD, 25 * kDEG2RAD, 15 * kDEG2RAD });

      auto pen_fn = std::make_shared<Regi2D3DPenaltyFnSE3Mag>();

      pen_fn->rot_pdfs_per_obj = { std::make_shared<FoldNormDist>(45 * kDEG2RAD, 45 * kDEG2RAD) };
      pen_fn->trans_pdfs_per_obj = { std::make_shared<NullDist>(1) };

      cmaes_regi->set_penalty_fn(pen_fn);
      cmaes_regi->set_img_sim_penalty_coefs(0.9, 0.1);
      
      regi.regi = cmaes_regi;
    }
  }  // New Level 1
  
  // New Level 2
  {
    vout << "  setting up all objects regi new level 2..." << std::endl;

    auto& lvl = ml_mo_regi->levels[prev_num_levels + 1];

    lvl.fixed_imgs_to_use = { 0 };

    lvl.ds_factor = 0.25;
    
    vout << "    setting up ray caster..." << std::endl;
    lvl.ray_caster = ray_cast_fact();
    
    vout << "    setting up sim metric..." << std::endl;
    lvl.sim_metrics = { MakeSimMetric(sim_metric_fact, lvl.ds_factor) };
    
    lvl.regis.resize(1);
    
    auto& regi = lvl.regis[0];

    regi.mov_vols    = { 0, 1, 2 };  // All objects simultaneous opt 
    regi.ref_frames  = { pelvis_app_ref_frame_idx,  // use the APP ref frames
                         fh_l_app_ref_frame_idx,
                         fh_r_app_ref_frame_idx };
    regi.static_vols = { };    // No other objects remain

    // use current estimates for each object as their initialization
    auto pelvis_init_guess_fn = std::make_shared<UseCurEstForInit>();
    pelvis_init_guess_fn->vol_idx = 0;
    
    auto left_femur_init_guess_fn = std::make_shared<UseCurEstForInit>();
    left_femur_init_guess_fn->vol_idx = 1;
    
    auto right_femur_init_guess_fn = std::make_shared<UseCurEstForInit>();
    right_femur_init_guess_fn->vol_idx = 2;

    regi.init_mov_vol_poses = { pelvis_init_guess_fn,
                                left_femur_init_guess_fn,
                                right_femur_init_guess_fn };
    
    auto bobyqa_regi = std::make_shared<Intensity2D3DRegiBOBYQA>();
    bobyqa_regi->set_opt_vars(se3_vars);
    bobyqa_regi->set_opt_x_tol(0.0001);
    bobyqa_regi->set_opt_obj_fn_tol(0.0001);
    bobyqa_regi->set_bounds({ 2.5 * kDEG2RAD, 2.5 * kDEG2RAD, 2.5 * kDEG2RAD, 2.5, 2.5, 2.5,
                              2.5 * kDEG2RAD, 2.5 * kDEG2RAD, 2.5 * kDEG2RAD, 2.5, 2.5, 2.5,
                              2.5 * kDEG2RAD, 2.5 * kDEG2RAD, 2.5 * kDEG2RAD, 2.5, 2.5, 2.5 } );

    regi.regi = bobyqa_regi;
  }  // New Level 2
}

void RunProjPreProcAndSetRegi(MultiLevelMultiObjRegi* ml_mo_regi,
                              const ProjDataF32& pd,
                              std::ostream& vout,
                              const bool no_log_remap = false)
{
  ProjPreProc proj_preproc;
  proj_preproc.params.crop_width = kPROJ_CROP_PIX;
  proj_preproc.params.no_log_remap = no_log_remap;

  proj_preproc.input_projs = { pd };

  proj_preproc.set_debug_output_stream(*ml_mo_regi);
  
  vout << "preprocessing projection..." << std::endl;
  proj_preproc();

  ml_mo_regi->fixed_proj_data = proj_preproc.output_projs;
  
  if (ml_mo_regi->debug_info)
  {  
    ml_mo_regi->debug_info->fixed_projs = proj_preproc.input_projs;
    
    ml_mo_regi->debug_info->proj_pre_proc_info = proj_preproc.params;
  }
}

}  // un-named

std::string xreg::IPCAIVersionStr()
{
  return IPCAI_PROJ_VERSION;
}

xreg::RayCasterFactoryFn xreg::ProgOptsLineIntRayCasterFactory(ProgOpts& po)
{
  return [&po] () { return LineIntRayCasterFromProgOpts(po); };
}

xreg::RayCasterFactoryFn xreg::ProgOptsDepthRayCasterFactory(ProgOpts& po)
{
  return [&po] () { return DepthRayCasterFromProgOpts(po); };
}

xreg::SimMetricFactoryFn xreg::ProgOptsSimMetricFactory(ProgOpts& po)
{
  return [&po] () { return PatchGradNCCSimMetricFromProgOpts(po); };
}

void xreg::RunGlobalPelvisRegi(MultiLevelMultiObjRegi* ml_mo_regi,
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
                               std::ostream& vout)
{
  vout << "setting up multi-level regi object.." << std::endl;

  // calling this with save_debug==true will allocate the debug object
  ml_mo_regi->set_save_debug_info(save_debug);

  RunProjPreProcAndSetRegi(ml_mo_regi, pd, vout);

  ml_mo_regi->vol_names = { "Pelvis" };

  ml_mo_regi->vols = { vol_att };
  
  // setup the camera reference frame which we optimize in
  auto cam_align_ref = std::make_shared<MultiLevelMultiObjRegi::CamAlignRefFrameWithCurPose>();
  cam_align_ref->vol_idx = 0;
  cam_align_ref->center_of_rot_wrt_vol = ITKVol3DCenterAsPhysPt(vol_att.GetPointer());
  cam_align_ref->cam_extrins = ml_mo_regi->fixed_proj_data[0].cam.extrins;

  ml_mo_regi->ref_frames = { cam_align_ref };

  vout << "computing nominal AP view to use for initialization..." << std::endl;
  const FrameTransform app_to_vol_std = AnteriorPelvicPlaneFromLandmarksMap(lands);

  const bool pat_is_up = *pd.rot_to_pat_up == ProjDataRotToPatUp::kZERO;

  const FrameTransform cam_to_naive_ap = CreateAPViewOfAPP(pd.cam, 0.8, true, pat_is_up);

  const FrameTransform cam_to_vol_approx_ap = app_to_vol_std * cam_to_naive_ap;

  ml_mo_regi->init_cam_to_vols = { cam_to_vol_approx_ap };

  if (method_id == 1)
  {
    SetupGlobalPelvisRegiMethod1(ml_mo_regi,
                                 lands,
                                 ray_cast_fact,
                                 sim_metric_fact,
                                 de_pop_size,
                                 pat_is_up,
                                 grid_batch_size,
                                 vout);
  }  
  else if (method_id == 2)
  {
    SetupGlobalPelvisRegiMethod2(ml_mo_regi,
                                 ray_cast_fact,
                                 sim_metric_fact,
                                 grid_batch_size,
                                 pso_num_parts,
                                 vout);
  }
  else
  {
    xregThrow("Unsupported global pelvis regi method id: %d", method_id);
  }

  if (save_debug)
  {
    vout << "  setting regi debug info..." << std::endl;

    ml_mo_regi->debug_info->vols = { vol_hu };
  }

  vout << "running regi..." << std::endl;
  ml_mo_regi->run();
  
  vout << "regi complete!" << std::endl;
}

void xreg::RunFemursAndPelvisRegiUsingPrevRegi(MultiLevelMultiObjRegi* ml_mo_regi,
                                               const FrameTransform& prev_pelvis_regi_cam_to_vol,
                                               std::vector<itk::Image<float,3>::Pointer>& vols_hu,
                                               std::vector<itk::Image<float,3>::Pointer>& vols_att,
                                               const LandMap3& lands,
                                               const ProjDataF32& pd,
                                               const bool save_debug,
                                               RayCasterFactoryFn ray_cast_fact,
                                               SimMetricFactoryFn sim_metric_fact,
                                               std::ostream& vout)
{
  vout << "setting up multi-level regi object.." << std::endl;

  // calling this with save_debug==true will allocate the debug object
  ml_mo_regi->set_save_debug_info(save_debug);

  RunProjPreProcAndSetRegi(ml_mo_regi, pd, vout);

  ml_mo_regi->vol_names = { "Pelvis", "L. Femur", "R. Femur" };

  ml_mo_regi->vols = vols_att;
  
  ml_mo_regi->init_cam_to_vols = { prev_pelvis_regi_cam_to_vol,
                                   prev_pelvis_regi_cam_to_vol,
                                   prev_pelvis_regi_cam_to_vol };
  
  AddSetupForFemursAndPelvisRegi(ml_mo_regi, lands, ray_cast_fact, sim_metric_fact, vout);

  if (save_debug)
  {
    vout << "  setting regi debug info..." << std::endl;

    ml_mo_regi->debug_info->vols = { vols_hu[0], vols_hu[1], vols_hu[2] };
  }

  vout << "running regi..." << std::endl;
  ml_mo_regi->run();
  
  vout << "regi complete!" << std::endl;
}

namespace  // un-named
{

using namespace xreg;

void SingleChannelHipSegRule(const std::vector<cv::Mat>& depth_imgs,
                             Proj3DLabelsTo2D::LabelImg2D* seg_img)
{
  xregASSERT(depth_imgs.size() == 7);
 
  std::array<const float*,6> depth_row_bufs;

  cv::Mat seg_ocv = ShallowCopyItkToOpenCV(seg_img);
  
  const int nr = seg_ocv.rows;
  const int nc = seg_ocv.cols;

  for (int r = 0; r < nr; ++r)
  {
    for (int l = 0; l < 6; ++l)
    {
      depth_row_bufs[l] = &depth_imgs[l+1].at<float>(r,0);
    }

    unsigned char* seg_row = &seg_ocv.at<unsigned char>(r,0);

    for (int c = 0; c < nc; ++c)
    {
      auto& s = seg_row[c];

      if (depth_row_bufs[4][c] < kRAY_CAST_MAX_DEPTH)  // left femur
      {
        s = 5;
      }
      else if (depth_row_bufs[5][c] < kRAY_CAST_MAX_DEPTH)  // right femur
      {
        s = 6;
      }
      else
      {
        const auto& left_hemi_pelvis_depth  = depth_row_bufs[0][c];
        const auto& right_hemi_pelvis_depth = depth_row_bufs[1][c];
        
        const bool intersects_left_hemi_pelvis  = left_hemi_pelvis_depth  < kRAY_CAST_MAX_DEPTH;
        const bool intersects_right_hemi_pelvis = right_hemi_pelvis_depth < kRAY_CAST_MAX_DEPTH;
        
        if (intersects_left_hemi_pelvis && intersects_right_hemi_pelvis)
        {
          s = (left_hemi_pelvis_depth < right_hemi_pelvis_depth) ? 1 : 2;
        }
        else if (intersects_left_hemi_pelvis)
        {
          s = 1;
        }
        else if (intersects_right_hemi_pelvis)
        {
          s = 2;
        }
        else if (depth_row_bufs[2][c] < kRAY_CAST_MAX_DEPTH)  // vertebra
        {
          s = 3;
        }
        else if (depth_row_bufs[3][c] < kRAY_CAST_MAX_DEPTH)  // upper sacrum
        {
          s = 4;
        }
        else
        {
          s = 0;
        }
      }
    }
  }
}

}  // un-named

std::tuple<itk::Image<unsigned char,2>::Pointer, std::vector<cv::Mat>>
xreg::ProjectHipLabels(itk::Image<unsigned char,3>* ct_labels,
                       const CameraModel& cam,
                       const FrameTransform& pelvis_pose,
                       const FrameTransform& left_femur_pose,
                       const FrameTransform& right_femur_pose,
                       RayCasterFactoryFn ray_cast_fact,
                       const bool verbose,
                       std::ostream& vout)
{
  vout << "setting fields of label projector..." << std::endl;
  
  Proj3DLabelsTo2D proj_labels;
  proj_labels.set_debug_output_stream(vout, verbose);

  proj_labels.depth_ray_caster = ray_cast_fact();
  
  proj_labels.cam = cam;

  proj_labels.label_vol = ct_labels;
  proj_labels.labels_to_proj = { 1, 2, 3, 4, 5, 6 };  // not including lower sacrum
  proj_labels.obj_poses = { pelvis_pose,  // 1
                            pelvis_pose,  // 2
                            pelvis_pose,  // 3
                            pelvis_pose,  // 4
                            left_femur_pose,  // 5
                            right_femur_pose  // 6
                          };

  proj_labels.convert_to_single_chan_seg_fn = SingleChannelHipSegRule;

  proj_labels.init();

  proj_labels.run();
  
  return std::make_tuple(proj_labels.single_chan_seg, proj_labels.seg_channels);
}

void xreg::RunPelvisAndFemursIntraopRegi(MultiLevelMultiObjRegi* ml_mo_regi,
                                         const int method_id,
                                         std::vector<itk::Image<float,3>::Pointer>& vols_hu,
                                         std::vector<itk::Image<float,3>::Pointer>& vols_att,
                                         const LandMap3& lands_3d,
                                         const ProjDataF32& pd,
                                         const cv::Mat* pd_seg,
                                         const bool save_debug,
                                         RayCasterFactoryFn ray_cast_fact,
                                         SimMetricFactoryFn sim_metric_fact,
                                         std::ostream& vout,
                                         const bool no_log_remap)
{
  vout << "setting up multi-level regi object.." << std::endl;
  
  vout << "Method ID: " << method_id << std::endl;

  const size_type num_det_lands_2d = pd.landmarks.size();

  vout << "# 2D Landmarks: " << num_det_lands_2d << std::endl;

  const bool method_2_use_posit = num_det_lands_2d >= 4;

  const bool pat_is_up = *pd.rot_to_pat_up == ProjDataRotToPatUp::kZERO;

  const bool use_patch_wgt = pd_seg && ((method_id == 2) || (method_id == 3));

  const auto* pd_seg_to_use = pd_seg;

  cv::Mat pd_seg_tmp;
    
  if (use_patch_wgt && !pat_is_up)
  {
    vout << "rotating patient seg to remove rotation applied to get patient up..." << std::endl;

    xregASSERT(pd_seg);

    pd_seg_tmp = pd_seg->clone();

    FlipImageColumns(&pd_seg_tmp);
    FlipImageRows(&pd_seg_tmp);

    pd_seg_to_use = &pd_seg_tmp;
  }

  using pixel_patch_wgt_fn = RayCastPixelScalar (const unsigned char&);

  auto compute_wgt_img = [pd_seg_to_use] (pixel_patch_wgt_fn f)
  {
    xregASSERT(pd_seg_to_use);

    cv::Mat wgt_img(pd_seg_to_use->rows, pd_seg_to_use->cols, cv::DataType<RayCastPixelScalar>::type);

    for (int r = 0; r < wgt_img.rows; ++r)
    {
      auto* wgt_row = &wgt_img.at<RayCastPixelScalar>(r,0);
      
      const auto* seg_row = &pd_seg_to_use->at<unsigned char>(r,0);

      for (int c = 0; c < wgt_img.cols; ++c)
      {
        // Just weight hemipelvis and femurs equally as 1, everything else 0
        wgt_row[c] = f(seg_row[c]);
      }
    }

    return ITKImageDeepCopy(ShallowCopyOpenCVToItk<RayCastPixelScalar>(wgt_img).GetPointer());
  };

  // calling this with save_debug==true will allocate the debug object
  ml_mo_regi->set_save_debug_info(save_debug);

  RunProjPreProcAndSetRegi(ml_mo_regi, pd, vout, no_log_remap);

  ml_mo_regi->vol_names = { "Pelvis", "L. Femur", "R. Femur" };

  ml_mo_regi->vols = vols_att;

  vout << "Determining initialization..." << std::endl;

  FrameTransform app_to_vol_pelvis = AnteriorPelvicPlaneFromLandmarksMap(lands_3d);
  app_to_vol_pelvis.matrix().block(0,3,3,1) = (lands_3d.find("FH-r")->second +
                                               lands_3d.find("FH-l")->second) / 2;

  const FrameTransform vol_to_app_pelvis = app_to_vol_pelvis.inverse();

  Pt3 ap_view_3d_ref_pt;
  Pt2 ap_view_2d_ref_pt;

  if (method_id == 1)
  {
    // For method 1, align the centroid of the L/R ASIS, SPS with the center of the image

    vout << "  Computing method 1 alignment points..." << std::endl;

    ap_view_3d_ref_pt = vol_to_app_pelvis *
                          ((lands_3d.find("ASIS-l")->second +
                            lands_3d.find("ASIS-r")->second +
                            lands_3d.find("SPS-l")->second  +
                            lands_3d.find("SPS-r")->second) / 4);

    ap_view_2d_ref_pt[0] = pd.cam.num_det_cols / 2.0;
    ap_view_2d_ref_pt[1] = pd.cam.num_det_rows / 2.0;  
  }
  else
  {
    // For method 3 align a detected landmark, this is also used by method 2 when less than
    // 4 landmarks are detected in 2D

    vout << "  Methods 2,3: looking for the landmark to use for AP view alignment..." << std::endl;

    // ordered by what we generally think are more accurate
    const std::vector<std::string> lands_names_for_align = { "FH-l",  "FH-r",
                                                             "IOF-l", "IOF-r",
                                                             "IPS-l", "IPS-r",
                                                             "MOF-l", "MOF-r",
                                                             "SPS-l", "SPS-r",
                                                             "GSN-l", "GSN-r",
                                                             "ASIS-l", "ASIS-r" };
    
    std::string land_to_use_for_ap_align;
    for (const auto& land_name : lands_names_for_align)
    {
      auto land_it_2d = pd.landmarks.find(land_name);

      if (land_it_2d != pd.landmarks.end())
      {
        land_to_use_for_ap_align = land_name;
        break;
      }
    }
    xregASSERT(!land_to_use_for_ap_align.empty());

    vout << "    using: " << land_to_use_for_ap_align << std::endl;
    
    ap_view_3d_ref_pt = vol_to_app_pelvis * lands_3d.find(land_to_use_for_ap_align)->second;
    ap_view_2d_ref_pt = pd.landmarks.find(land_to_use_for_ap_align)->second;
  }
    
  vout << "     3D ref pt: "
       << fmt::format("[ {} , {} , {} ]",
                      ap_view_3d_ref_pt[0],
                      ap_view_3d_ref_pt[1],
                      ap_view_3d_ref_pt[2])
       << std::endl;
    
  vout << "    2D ref idx: "
       << fmt::format("[ {} , {} ]", ap_view_2d_ref_pt[0], ap_view_2d_ref_pt[1])
       << std::endl;

  const FrameTransform ref_cam_to_app = CreateAPViewOfAPP(pd.cam, 0.85, true, pat_is_up,
                                                          ap_view_3d_ref_pt, ap_view_2d_ref_pt);
  
  FrameTransform pelvis_init_cam_to_vol;
  
  if ((method_id == 1) || (method_id == 3) || !method_2_use_posit)
  {
    vout << "  nominal AP view init..." << std::endl;
    pelvis_init_cam_to_vol = app_to_vol_pelvis * ref_cam_to_app;
    
    vout << pelvis_init_cam_to_vol.matrix() << std::endl;
  }
  
  if (method_id == 2)
  {
    vout << "  method 2..." << std::endl;

    if (!method_2_use_posit)
    {
      vout << "    running iterative landmark re-proj dist. minimization "
              "starting from nominal AP view..." << std::endl;
      pelvis_init_cam_to_vol = PnPReprojCMAES(pd.cam, lands_3d, pd.landmarks, pelvis_init_cam_to_vol);
    }
    else
    {
      vout << "    using POSIT as initialization for iterative landmark re-proj "
              "dist. minimization..." << std::endl;
      pelvis_init_cam_to_vol = PnPPOSITAndReprojCMAES(pd.cam, lands_3d, pd.landmarks);
    }

    vout << pelvis_init_cam_to_vol.matrix() << std::endl;
  }

  ml_mo_regi->init_cam_to_vols = { pelvis_init_cam_to_vol,
                                   pelvis_init_cam_to_vol,
                                   pelvis_init_cam_to_vol };

  // setup the camera reference frame which we optimize in for the first two pelvis only registrations
  auto cam_align_ref = std::make_shared<MultiLevelMultiObjRegi::CamAlignRefFrameWithCurPose>();
  cam_align_ref->vol_idx = 0;
  cam_align_ref->center_of_rot_wrt_vol = ITKVol3DCenterAsPhysPt(ml_mo_regi->vols[0].GetPointer());
  cam_align_ref->cam_extrins = ml_mo_regi->fixed_proj_data[0].cam.extrins;

  ml_mo_regi->ref_frames = { cam_align_ref };
  
  // se(3) lie algebra vector space for optimization
  auto se3_vars = std::make_shared<SE3OptVarsLieAlg>();
 
  if (save_debug)
  {
    ml_mo_regi->debug_info->regi_names = { { "Pelvis1" }, { "Pelvis2" } };
  }

  ml_mo_regi->levels.resize(2);
  
  // Level 1
  {
    vout << "  setting up pelvis regi level 1..." << std::endl;

    auto& lvl = ml_mo_regi->levels[0];

    lvl.fixed_imgs_to_use = { 0 };
    
    lvl.ds_factor = 0.125;

    vout << "    setting up ray caster..." << std::endl;
    lvl.ray_caster = ray_cast_fact();
    
    vout << "    setting up sim metric..." << std::endl;
    lvl.sim_metrics = { MakeSimMetric(sim_metric_fact, lvl.ds_factor) };
    
    if (use_patch_wgt)
    {
      vout << "computing weight image and setting patch weights for pelvis regis..." << std::endl;

      auto* sm = dynamic_cast<ImgSimMetric2DPatchCommon*>(lvl.sim_metrics[0].get());
    
      auto pelvis_patch_wgt_fn = [] (const unsigned char& s)
      {
        return ((s == 1) || (s == 2)) ? RayCastPixelScalar(1) : RayCastPixelScalar(0);
      };
      
      sm->set_wgt_img(compute_wgt_img(pelvis_patch_wgt_fn));
    }

    // CMA-ES
    lvl.regis.resize(1);

    auto& regi = lvl.regis[0];

    regi.mov_vols    = { 0 };  // pelvis vol pose is optimized over
    regi.ref_frames  = { 0 };  // use ref frame that is camera aligned
    regi.static_vols = { };    // No other objects exist

    auto init_guess_fn = std::make_shared<UseCurEstForInit>();
    init_guess_fn->vol_idx = 0;
    regi.init_mov_vol_poses = { init_guess_fn };

    auto cmaes_regi = std::make_shared<Intensity2D3DRegiCMAES>();
    cmaes_regi->set_opt_vars(se3_vars);
    cmaes_regi->set_opt_x_tol(0.01);
    cmaes_regi->set_opt_obj_fn_tol(0.01);
    cmaes_regi->set_pop_size(100);
    cmaes_regi->set_sigma({ 15 * kDEG2RAD, 15 * kDEG2RAD, 30 * kDEG2RAD, 50, 50, 100 });

    if (method_id == 2)
    {
      auto pen_fn = std::make_shared<Regi2D3DPenaltyFnSE3EulerDecomp>();
      pen_fn->rot_x_pdf   = std::make_shared<NormalDist1D>(0, 15 * kDEG2RAD);
      pen_fn->rot_y_pdf   = std::make_shared<NormalDist1D>(0, 15 * kDEG2RAD);
      pen_fn->rot_z_pdf   = std::make_shared<NormalDist1D>(0, 10 * kDEG2RAD);
      pen_fn->trans_x_pdf = std::make_shared<NormalDist1D>(0, 30);
      pen_fn->trans_y_pdf = std::make_shared<NormalDist1D>(0, 30);
      pen_fn->trans_z_pdf = std::make_shared<NormalDist1D>(0, 150);

      cmaes_regi->set_penalty_fn(pen_fn);
      cmaes_regi->set_img_sim_penalty_coefs(0.9, 0.1);
    }
    else if (method_id == 3)
    {
      auto pen_fn = std::make_shared<Regi2D3DPenaltyFnLandReproj>();
      
      // Adjust the 2D landmarks for cropping and downsampling
      LandMap2 lands_2d_lvl = pd.landmarks;

      for (auto& lkv : lands_2d_lvl)
      {
        lkv.second[0] = (lkv.second[0] - kPROJ_CROP_PIX) * lvl.ds_factor;
        lkv.second[1] = (lkv.second[1] - kPROJ_CROP_PIX) * lvl.ds_factor;
      }

      pen_fn->set_lands(lands_3d, lands_2d_lvl);

      cmaes_regi->set_penalty_fn(pen_fn);
      cmaes_regi->set_img_sim_penalty_coefs(0.9, 0.1);
    }
    
    regi.regi = cmaes_regi;
  }  // Level 1
  
  // Level 2
  {
    vout << "  setting up pelvis regi level 2..." << std::endl;

    auto& lvl = ml_mo_regi->levels[1];

    lvl.fixed_imgs_to_use = { 0 };
    
    lvl.ds_factor = 0.25;

    vout << "    setting up ray caster..." << std::endl;
    lvl.ray_caster = ray_cast_fact();
    
    vout << "    setting up sim metric..." << std::endl;
    lvl.sim_metrics = { MakeSimMetric(sim_metric_fact, lvl.ds_factor) };
    
    // BOBYQA
    lvl.regis.resize(1);

    auto& regi = lvl.regis[0];

    regi.mov_vols    = { 0 };  // pelvis vol pose is optimized over
    regi.ref_frames  = { 0 };  // use ref frame that is camera aligned
    regi.static_vols = { };    // No other objects exist

    auto init_guess_fn = std::make_shared<UseCurEstForInit>();
    init_guess_fn->vol_idx = 0;
    regi.init_mov_vol_poses = { init_guess_fn };

    auto bobyqa_regi = std::make_shared<Intensity2D3DRegiBOBYQA>();
    bobyqa_regi->set_opt_vars(se3_vars);
    bobyqa_regi->set_opt_x_tol(0.0001);
    bobyqa_regi->set_opt_obj_fn_tol(0.0001);
    bobyqa_regi->set_bounds({ 2.5 * kDEG2RAD, 2.5 * kDEG2RAD, 2.5 * kDEG2RAD,
                              5, 5, 10 });
    
    regi.regi = bobyqa_regi;
  }  // Level 2

  AddSetupForFemursAndPelvisRegi(ml_mo_regi, lands_3d, ray_cast_fact, sim_metric_fact, vout);
  
  if (use_patch_wgt)
  {
    vout << "computing weight image and setting patch weights for femur regis..." << std::endl;

    auto* sm = dynamic_cast<ImgSimMetric2DPatchCommon*>(ml_mo_regi->levels[2].sim_metrics[0].get());
  
    auto pelvis_femurs_patch_wgt_fn = [] (const unsigned char& s)
    {
      return ((s == 1) || (s == 2) || (s == 5) || (s == 6)) ? RayCastPixelScalar(1) : RayCastPixelScalar(0);
    };
    
    sm->set_wgt_img(compute_wgt_img(pelvis_femurs_patch_wgt_fn));
  }

  if (save_debug)
  {
    vout << "  setting regi debug info..." << std::endl;

    ml_mo_regi->debug_info->vols = { vols_hu[0], vols_hu[1], vols_hu[2] };
  }

  vout << "running regi..." << std::endl;
  ml_mo_regi->run();
  
  vout << "regi complete!" << std::endl;
}

