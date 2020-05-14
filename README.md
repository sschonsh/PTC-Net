# PTC-Net
THE TENSORFLOW CODES PRESENTED HERE ONLY WORK IN TF1.X. 

TranportDemo is a MatLab script which demonstrates the algorithmic process of parallel transporting compactly supported kernels in matlab.

ManifoldKernelTransformations shows transportations, rotations and dilations of a kernel on a closed surface.

DefineSurface.m is a  MatLab script which generates a surf structure which can be loaded into tensorflow-to use a different surfaces just change the surf.pts and surf.trg to match your data

TFMNIST1.py shows a simple MNIST classifier with archetecture x -> PTC(16) -> FC(10) -> y. This script requires the ExSurf.mat which DefineSurfaces creates to be in your working directory.

The MatLab fast marching toolbox is required for some matlab scripts (https://www.mathworks.com/matlabcentral/fileexchange/6110-toolbox-fast-marching)

