library(gkmSVM)
posfn='/home/hsher/encore_region_shrink_ri/RBFOX2_HepG2.rep1.fasta'
negfn='/home/hsher/encore_region_shrink_ri/RBFOX2_HepG2.rep1.ctrl.fasta'
testfn='/home/hsher/encore_region_shrink_ri/RBFOX2_HepG2.rep1.test.fasta'
kernelfn='/home/hsher/encore_region_shrink_ri/RBFOX2_HepG2.kernel.txt'
svmfnprfx= '/home/hsher/encore_region_shrink_ri/RBFOX2_HepG2.rep1.svm'
outfn='/home/hsher/encore_region_shrink_ri/RBFOX2_HepG2.rep1.test.out'
gkmsvm_kernel(posfn, negfn, kernelfn)
gkmsvm_train(kernelfn, posfn, negfn, svmfnprfx)
gkmsvm_classify(testfn, svmfnprfx, outfn); #scores test sequences