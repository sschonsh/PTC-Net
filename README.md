# PTC-Net
Some Demo Code for PTC Networks: from https://arxiv.org/abs/1805.07857

TranportDemo is a matlab script which demonstrates the algormitmic process of parallel transporting compactly supported kenels

DefineSurface is a matlab script which generates a surf structure which can be loaded into tensorflow-to use a differnt surfaces just change the surf.pts and surf.trg to match your data

TFMNIST1 shows a simple MNIST classifier with archetecture x -> PTC(16) -> FC(10) -> y