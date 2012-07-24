Identifying fMRI Model Violations 
=================================
with Lagrange Multiplier Tests
------------------------------

(fmriVIOLA software version 0.2)

Author: Ben Cassidy

------------------------

This software comprises tools for identifying model violation in Functional Magnetic Resonance Imaging analysis. 

The diagnostic tests are available as a MATLAB toolbox, and can be run from the MATLAB command line or incorporated in batch scripts as part of a full fMRI analysis pipeline.  

The current version includes tests for violation of the following model assumptions:

- Double Gamma (canonical) hemodynamic response function
- Double Gamma and temporal derivative HRF (derivative as defined in SPM software)
- Double Gamma and temporal derivative HRF (derivative as defined in FSL software)
- Time-varying HRF
- Non-linear HRF

------------------------
Basic instructions:
1- Add folders to matlab path via 
    -> File -> Set Path -> Add With Subfolders

2- Load LM_testData.mat from the 'examples' folder

3- Run the example script "example_fmriTest.m"

4- Test your own data!

------------------------
Publications
============

Full details of the tests are available in:

B. Cassidy, CJ. Long, C. Rae, V. Solo, 
"Identifying fMRI Model Violations with Lagrange Multiplier Tests", 
IEEE Transactions on Medical Imaging, in press.

http://dx.doi.org/10.1109/TMI.2012.2195327

B. Cassidy, V. Solo, 
"fMRI Model Diagnostics for the Double Gamma and Temporal Derivative", 
IEEE International Symposium on Biomedical Imaging, 2012.

http://neura.edu.au/research/facilities/neura-imaging-centre/software

------------------------

Technical notes
===============

The OOP software structure allows simple addition of new tests for other model violations (fMRI or otherwise).

The object hierarchy is
	
	LM_test
	   |
    LM_test_fmri
	   |
    anomalyTest

where anomalyTest is either LM_DG, LM_NL or LM_TV depending on the
desired test.  The tests are implemented as objects, so the user 
first defines an anomaly test object with default properties, 
modifies those properties for individual preference, then runs the 
test. 

Each different test comprises a new test object, though the data 
and input stimulus signal within each test can be modified easily 
to facilitate batch processing using otherwise identical model 
assumptions.

Default test definitions take an input of form:

    TEST_OBJ = SPECIFIC_TEST_OBJECT
                (data, ...
                 input_stim_init, ...
                 TR, ...
                 <LM_DG test: optional> DG_DERIV, ...
                 <ALL: optional> extraDriftVars)

where

    data            = type=double  : dim (#TimePoints x #Voxels)
                         -- fMRI data matrix (2D)
    input_stim_init = type=logical : dim (#TimePoints x 1)
                         -- input stimulus signal
    TR              = type=double  : scalar
                         -- Sampling repetition time (seconds)

The Double Gamma HRF test has an additional optional input to test the
Double Gamma & Temporal Derivative models from either FSL or SPM:
    
    DG_DERIV        = type=string  : 'SPM', 'FSL', false
                         -- optional, test validity of DG+derivative model

however for all tests only data, input_stim_init and TR are mandatory.

Known issues/ improvements needed:
--------------------------------------------------------

Requires MATLAB 2010b or newer.

Currently requires fMRI data input as a 2D matrix (time x voxels). Future versions will include support for SPM data structures and additional tools for reporting and interpreting results.  

Output anomaly test statistic is in the same shape as the input data, so will need reshaping to view the anomaly map (although the map is instructive, the presence of model anomaly anywhere within an fMRI data volume is cause for concern; see http://dx.doi.org/10.1109/TMI.2012.2195327).

This version selects model orders (AR and regression) separately, as a first approximation to a cyclic descent procedure. A git development branch is available (origin/simulParaEst/) which estimates these parameters concurrently (although with added performance penalty).



License
=======

This software is released under the BSD 3-clause license.

Copyright (c) 2012, authors of fmriVIOLA
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
* Neither the name of the fmriVIOLA authors nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

