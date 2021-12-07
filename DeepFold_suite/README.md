# DeepFold

DeepFold is a program for ab initio protein structure prediction using deep learning restraints
developed by Robin Pearce (robpearc@umich.edu) at the Zhang Lab 

The DeepFold online server can be found at https://zhanggroup.org/DeepFold/

           INSTALLATION AND IMPLEMENTATION OF THE DEEPFOLD SUITE
   (Copyright 2021 by Zhang Lab, University of Michigan, All rights reserved)
                    (Version 1.0, 12/1/2021)

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT 
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


   #######################################################
   #                                                     #
   #  1. What is the DeepFold Suite?                     #
   #                                                     #
   #######################################################
   
   The DeepFold Suite is a composite package of programs for ab initio 
   protein structure prediction. The Suite includes the following programs:

   a) DeepFold: A rapid protein folding program 
   b) DeepMSA: A program for multiple sequence alignment generation. Note, in the
      DeepFold manuscript, DeepMSA2 was used, but given the size of the sequence
      databases required to run DeepMSA2, we have included DeepMSA version 1, which
      is faster and can be run using smaller sequence databases.
   c) DeepPotential: A program for deep learning spatial restraint generation
   d) ModRefiner: A program for protein structure refinement
   e) FASPR: A program for rapid side-chain packing

   #######################################################
   #                                                     #
   #  2. How to install the DeepFold Suite?              #
   #                                                     #
   #######################################################
   
   a) download the DeepFold Suite 'DeepFold_suite.tar.bz2' 
      and unpack 'DeepFold_suite.tar.bz2 by
      > tar -xvf DeepFold_suite.tar.bz2
      The root path of this package is called the $pkgdir, e.g. 
      /home/yourname/DeepFold_suite. You should have all the programs in this 
      directory. You can install the package at any location on your computer.
   
   b) Third-party software installation:

      While the majority of the programs in the package were developed in the Zhang Lab, 
      wherein the permission for use is released, there are some programs and databases 
      (including the nr, uniclust30, uniref90 and metaclust databases) that were developed 
      by third-party groups. It is the user's obligation to obtain license permission from 
      the developers for all the third-party software before using them. In addition, your 
      system needs to have Java, Perl, and python3 (which supports pytorch 0.3.0) with numpy, 
      scipy and PyTorch installed.

      To use DeepMSA, you need to download uniclust30, uniref90 and metaclust from 
      http://gwdu111.gwdg.de/~compbiol/uniclust/2017_04/uniclust30_2017_04_hhsuite.tar.gz ,      
      ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz ,
      and https://metaclust.mmseqs.org/2017_05/metaclust_2017_05.fasta.gz. You should unpack the folders.
      Then use $pkgdir/bin/hhsuite2/bin/esl-sfetch to create .ssi index files for uniref90 and metaclust, here 
      $pkgdir is the path where you put the DeepFold suite. For example, if the uniref90 database in the 
      uniref90 folder is named uniref90.fasta, then go to the uniref90 folder and run 
      $pkgdir/bin/hhsuite2/bin/esl-sfetch --index uniref90.fasta. You will find a new file named uniref90.fasta.ssi 
      after the command is finished. Then do the same thing for the metaclust database. When you run the program 
      runDeepFoldPipeline.pl, it takes as argument the path to the sequence databases as outlined in section 4.2.

      A script 'download_lib.pl' is provided in the package to automatically download these databases.

   #######################################################
   #                                                     #
   # 3. Bug report:                                      #
   #                                                     #
   #######################################################
                                                  
   Please report bugs and suggestions to robpearc@umich.edu


   #######################################################
   #                                                     #
   #  4. Installation and implementation of DeepFold     #
   #                                                     #
   #######################################################
   
4.1. Introduction to DeepFold
   
   DeepFold is a new method for ab initio protein structure prediction. The pipeline uses DeepMSA to search
   multiple whole-genome and metagenomic sequence databases to generate an MSA. From the generated MSA,
   DeepPotential is applied to predict spatial restraints including distance and contact maps as well as
   inter-residue orientations. Then, DeepFold generates predicted models using rapid LBFGS folding simulations
   based on the predicted restraints from deep learning. 

4.2. How to run DeepFold?
   
   a) The main script for running DeepFold is $pkgdir/runDeepFoldPipeline.pl. 
      Running it directly without arguments will output the help information.

   b) The following arguments must be set to run DeepFold: 

      "$pkgdir/runDeepFoldPipeline.pl -datadir=<path to data directory> -generate_msa=<msa flag> -run_refine=<refine flag> -nr_path=<path to the nr database> -uniclust_path=<path to the uniclust database> -uniref_path=<path to the uniref database> -metaclust_path=<path to the metaclust database>"

      -datadir          the directory that contains your input sequence (seq.txt)
                        and where the output files will be stored
      
      -generate_msa     (True/False) whether or not to run DeepMSA. This is
                        set to True by default. If you want to use a pre-generated
                        MSA, it should in a file at <datadir>/MSA/protein.aln.
                
      -run_refine       (True/False) whether or not to run ModRefiner to refine
                        the output structure. This is set to False by default
    
      -nr_path          path to the nr database required for secondary structure
                        and backbone torsion angle prediction.
      
      -uniclust_path    path to the directory where the uniclust database (uniclust30_2017_04_a3m.ffdata) 
                        is installed. This argument is required if generate_msa is set to True.

      -uniref_path      path to the directory where the uniref database (uniref90.fasta) is installed.
                        This argument is required if generate_msa is set to True.
   
      -metaclust_path   path to the directory where the metaclust database (metaclust.fasta) is installed.
                        This argument is required if generate_msa is set to True.

   NOTE:
   a) Outline of the steps involved in running DeepFold by 'runDeepFoldPipeline.pl':
      a1) $pkgdir/scripts/1_predPhiPsi.pl - predict phi/psi angles by Anglor and the secondary structure
          by PSSpred.
      a2) $pkgdir/scripts/2_runDeepMSA.pl - run DeepMSA to generate the input MSA.
      a3) $pkgdir/scripts/3_runDeepPotential.pl - predict spatial restraints by DeepPotential.
      a4) $pkgdir/scripts/4_runDeepFold.pl - run DeepFold to generate a protein structure.

4.3 System requirement:

   a) x86_64 machine, Linux kernel OS.
   b) Perl and java interpreters should be installed. 
   c) Basic compress and decompress package should be installed to support: 
      tar and bunzip2.
   d) Python3 with numpy, scipy and PyTorch

NOTE:
   There are several executable files under the '$pkgdir/distance/DeepPotential_04132020/bin' path.
   You may also need to set executable permission for the bin file. 
   'chmod +x -R $pkgdir/distance/DeepPotential_04132020/bin'

4.4. How to cite the DeepFold Suite?

   1. Robin Pearce, Yang Li, Gilbert S. Omenn, Yang Zhang. Fast and Accurate Ab Initio 
      Protein Structure Prediction Using Deep Learning Spatial Restraints (submitted).

