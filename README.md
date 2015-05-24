MARTIAL: MAnifold Reconciliation Through Iterative ALignment algorithm

----- Algorithm summary -----

This code provides an implementation of the MAnifold Reconciliation Through Iterative ALignment (MARTIAL) algorithm for multiclass domain adaptation.  

----- System Requirements -----

Tested using Matlab 7.14.0.739 (R2012a) on a Macbook Pro running OSX 10.6.

----- Disclaimer -----

Although this code has been reasonably well-tested, it is research code, and may contain bugs. Please refer to the FAQ below for info on commonly-occurring issues, and feel free to contact the author (bbue@alumni.rice.edu) if you have any difficulties using the code. 

----- Installation -----

1) First install the RelTrans framework, available at the following url:
   http://www.ece.rice.edu/~bdb1/#reltrans

2) Then run the following code to add the MARTIAL functions to your path:

  >> MARTIAL_ROOT='/path/to/MARTIAL/';
  >> addpath(genpath(MARTIAL_ROOT));
  >> savepath; % optional

----- External Libraries -----

The MARTIAL algorithm uses the following external libraries (provided in MARTIAL_ROOT/external):

1) Kabsch: Kabsch algorithm for computing the minimum Root Mean Square Deviation (RMSD) transformation between two point sets of equal cardinality and dimensionality.
  - Credit: Ehud Schreiber (schreib@compugen.co.il), 2012
  - Original implementation: http://www.mathworks.com/matlabcentral/fileexchange/25746-kabsch-algorithm/content/Kabsch.m
  - Modifications made by B. Bue, 2013: 
    + added zero_mean argument to toggle shifting point sets to zero mean
    + set default zero_mean argument to 0=disabled 

2) munkres: Hungarian algorithm for solving linear assignment problems. 
  - Credit: Yi Cao (y.cao@cranfield.ac.uk), 2008
  - Original implementation: http://www.mathworks.com/matlabcentral/fileexchange/20652-hungarian-algorithm-for-linear-assignment-problems-v2-3)


----- Example Usage and Output -----

The following code runs the MARTIAL algorithm to map labeled source domain data (S,SL) to the target domain represented by target test set (T,TL) using the unlabled source samples SU and target samples TU. Target labels TL for test set T only used in validation (and not for pivot selection).

  >> martial_options = struct('Nrand',Nrand,'Npool',Npool,'verbose',verbose)
  >> outparm = martial(S,SL,T,TL,SU,TU,Nkeep,martial_options);

Examine per-class RMSD and number of pivots used in each alignment step
  >> fprintf('Per class RMSD:\n')
  >> fprintf('init (%0.3f): ',mean(outparm.PSXrmsdinit)); fprintf('%0.3f ',outparm.PSXrmsdinit); fprintf('\n'); 
  >> fprintf('align (%0.3f): ',mean(outparm.PSXrmsdalign)); fprintf('%0.3f ',outparm.PSXrmsdalign); fprintf('\n'); 
  >> fprintf('improve (%0.3f): ',mean(outparm.PSXrmsdimprov)); fprintf('%0.3f ',outparm.PSXrmsdimprov); fprintf('\n'); 
  >> fprintf('Per class alignment:\n')
  >> fprintf('align (%0.3f): ',mean(outparm.PSXnalign)); fprintf('%0.3f ',outparm.PSXnalign); fprintf('\n'); 
  >> fprintf('improve (%0.3f): ',mean(outparm.PSXnimprove)); fprintf('%0.3f ',outparm.PSXnimprove); fprintf('\n'); 

Evaluate cross-validated prediction accuracy using the target-transformed source data produced by each of the alignment steps using a classifier of your hoosing (xvalidate function not included in this code): 
  >> STinit = outparm.SXinit; % source spectra after MARTIAL init step
  >> STalign = outparm.SXalign; % source spectra after MARTIAL align step
  >> STimprov = outparm.SXimprov;  % source spectra after MARTIAL improve ste
  >> STinitacc = xvalidate(STinit,SL,T,TL,xvalparm);
  >> STalignacc = xvalidate(STalign,SL,T,TL,xvalparm);
  >> STimprovacc = xvalidate(STimprov,SL,T,TL,xvalparm);  
  >> fprintf('Accuracy init: %0.3f align: %0.3f improve: %0.3f\n',...
                  STTinitacc,STTalignacc,STTimprovacc)

Output for 5-class Cuprite Hyp11->Av97 problem, Nkeep=25, Nrand=500, Npool=300 (baseline accuracy 0.9375): 
  Per class RMSD:
  init (0.370): 0.252 0.287 0.500 0.371 0.441 
  align (0.370): 0.252 0.287 0.500 0.371 0.441 
  improve (0.593): 0.494 0.529 0.690 0.589 0.663 
  Per class alignment:
  align (25.000): 25.000 25.000 25.000 25.000 25.000 
  improve (25.000): 25.000 25.000 25.000 25.000 25.000 
  Accuracy init: 0.976 align: 0.976 improve: 0.977

----- FAQ -----

Please refer to the FAQ provided in the RelTrans README (included in the RelTrans package, and also available at: http://www.ece.rice.edu/~bdb1/code/README_RelTrans.txt) for the most up-to-date list of frequently asked questions. 

----- Citation -----

Please cite the following publications if you publish any works using this
code:

  B. D. Bue and C. Jermaine. "Multiclass Domain Adaptation with Iterative Manifold Alignment."Workshop on Hyperspectral Image and Signal Processing: Evolution in Remote Sensing (WHISPERS), June 2013.

----- Additional References -----

Much of the early development of the MARTIAL algorithm is based on the RelTrans (RELational class knowledge TRANSfer) algorithm, detailed in the following publications:

  B. D. Bue and D. R. Thompson, “Multiclass Continuous Correspondence Learning,” NIPS Domain Adaptation Workshop, Dec. 2011.
  B. D. Bue, E. Merényi, and B. Csathó, “An Evaluation of Class Knowledge Transfer from Real to Synthetic Imagery,” Workshop on Hyperspectral Image and Signal Processing: Evolution in Remote Sensing (WHISPERS), Jun. 2011.
  B. D. Bue and E. Merényi, “Using spatial correspondences for hyperspectral knowledge transfer: evaluation on synthetic data,” Workshop on Hyperspectral Image and Signal Processing: Evolution in Remote Sensing (WHISPERS), Jun. 2010.

----- Changelog -----

06/17/13 - initial release.

----- Contact -----

Please contact the author (bbue@alumni.rice.edu) if you have any questions
regarding this code.

----- Copyright -----

All code described above with the exception of the external libraries copyright 2012-2013 Rice University.
