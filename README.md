# Regi2D3D-IPCAI2020
This repository contains several implementations of the 2D/3D registration strategies described in the IPCAI/IJCARS 2020 paper "Automatic Annotation of Hip Anatomy in Fluoroscopy for Robust and Efficient 2D/3D Registration."
The paper may be found [here](https://arxiv.org/abs/1911.07042) or [here](https://doi.org/10.1007/s11548-020-02162-7).
This repository is a companion to [DeepFluoroLabeling-IPCAI2020](https://github.com/rg2/DeepFluoroLabeling-IPCAI2020), which provides implementations of the PyTorch models and references to the entire dataset used in the IPCAI/IJCARS paper.

The global, semi-automatic, approaches are inappropriate for intraoperative purposes as they often require several minutes to complete processing.
However, these *offline* techniques are useful for creating a dataset that may be used for training a convolutional neural network (CNN), capable of performing these annotations very quickly.
An *online* registration strategy may then use the CNN-derived annotations to calculate automatic and robust registrations with intraoperatively compatible runtimes.

The following tools are provided by this repository:
  * Data extraction from the HDF5 dataset:
    * [Extraction of 3D data (CT/CT segmentation/3D landmarks)](apps/extract_3d_data_from_large_h5)
    * [Extraction of 2D fluoroscopy](apps/extract_proj_data_from_large_h5)
  * Creating annotations for CNN-training:
    * [Global pelvis registration](apps/global_pelvis_regi)
    * [Pelvis and femurs registration](apps/global_femurs_pelvis_regi)
    * [Projection of 3D annotations into 2D](apps/proj_3d_labels_into_2d)
  * [Intraoperative Registration](apps/intraop_pelvis_femurs_regi)

Python scripts are provided ([`extract_fcsv_from_nn_csv.py`](extract_fcsv_from_nn_csv.py) and [`extract_seg_from_nn_h5.py`](extract_seg_from_nn_h5.py)) in order to extract segmentations and landmarks from CNN inferences.
Instructions for performing CNN training and testing may be found [here](https://github.com/rg2/DeepFluoroLabeling-IPCAI2020/tree/master/train_test_code).

The tools provided here rely on the [xReg library](https://github.com/rg2/xreg) and users should see the [xReg wiki](https://github.com/rg2/xreg/wiki) for details on building the software.

Demonstrations detailing the usage of the tools listed above are provided on the [wiki](https://github.com/rg2/Regi2D3D-IPCAI2020/wiki) of this repository. 

The tools provided here are intended to demonstrate simplified use-cases of these registration strategies.
As such, users are encouraged to extend, modify and adapt these programs in order to conduct large scale studies efficiently.

Please submit an [issue](https://github.com/rg2/Regi2D3D-IPCAI2020/issues) for any problems, feature requests, or suggestions.

## License and Attribution
The software is available for use under the [MIT License](LICENSE).

If you have found this software useful in your work, we kindly ask that you cite the IPCAI/IJCARS paper:
```
Grupp, Robert B., et al. "Automatic annotation of hip anatomy in fluoroscopy for robust and efficient 2D/3D registration." International Journal of Computer Assisted Radiology and Surgery (2020): 1-11.
----------------------------------------------------------------------
@article{grupp2020automatic,
  title={Automatic annotation of hip anatomy in fluoroscopy for robust and efficient {2D}/{3D} registration},
  author={Grupp, Robert B and Unberath, Mathias and Gao, Cong and Hegeman, Rachel A and Murphy, Ryan J and Alexander, Clayton P and Otake, Yoshito and McArthur, Benjamin A and Armand, Mehran and Taylor, Russell H},
  journal={International Journal of Computer Assisted Radiology and Surgery},
  pages={1--11},
  publisher={Springer}
}
```
