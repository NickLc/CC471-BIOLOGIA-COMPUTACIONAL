
################################################################################

HMMVE 1.2 

Developed By:   Bioinformatics, Data Mining, and Machine Learning Lab (BDM)
		Department of Computer Science
		University of Missouri, Columbia

If you have questions or suggestions, please contact:        
		Jianlin Cheng, PhD
		Bioinformatics, Data Mining, and Machine Learning Lab (BDM)
		Department of Computer Science
		University of Missouri, Columbia
		Email: chengji@missouri.edu

################################################################################

HMMVE 1.2    Nov, 2009
-----------------

Introduction
-------------
HMMVE is a visual editor for profile Hidden Markov Model. With HMMVE, you can
  * View profile Hidden Markov Model in different views, save it in JPG of PNG format
  * Edit profile Hidden Markov Model
  * View viterbi path of a sequence against a profile Hidden Markov Model


Install and Run
-------------
1. Download and extract the package HMMVE_1.2.zip
2. Build HMMVE using the command
    ant
   You must have ant installed before building HMMVE, for more information about ant, see:
        http://ant.apache.org
   You must use JDK 1.6+ to compile and run HMMVE 1.2
3. Run HMMVE using
    HMMVE.sh (Unix)
    HMMVE.bat (Windows)


Change
-------------
* added HMMER3 support (model manipulation only, no alignments)

Documentation and Support
-------------
Please check http://bioserver.cs.ucf.edu/hmmeditor


Known Limitations in this Release
-------------
* Multi-domain match is not supported in this version
* HMMVE 1.x is only for protein sequence, DNA sequence is not fully tested
