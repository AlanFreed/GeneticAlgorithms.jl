(* ------|---------|---------|---------|---------|---------|---------|------- *)
(*       BBBB      EEEEE     L         The                                    *)
(*       B   B     E         L           BIOLOGICAL                           *)
(*       BBBB      EEE       L           ENGINEERING                          *)
(*       B    B    E         L           LABORATORY                           *)
(*       BBBBB     EEEEEE    LLLLLL        @ Saginaw Valley State University  *)
(* ------|---------|---------|---------|---------|---------|---------|------- *)

This is version 4.1 of April, 2013, for BEL written by Alan D. Freed.

BEL was compiled on a computer running:
Mac OS X version 10.7.2
Mono JIT compiler, version 2.10.8    (see http://www.mono-project.com).
Zonnon compiler,   version 1.3.0.0   (see http://www.zonnon.ethz.ch).

BEL is a suite of libraries for the .NET/Mono platforms written in Zonnon.
BEL is released under Version 3 of the GNU Lesser General Public License
(http://www.gnu.org/licences/gpl.html) for libraries.  A copy of this license
agreement can be found in the <docs> directory of this distribution, and is
also attached as an appendix to the Users Guide.  The author has tested this
suite on the both the .NET and Mono platforms.

BEL is a platform for developing computational programs and applications. 
It is similar to the Oberon packages CAT and CAPO that the author wrote. 
Contributions that are in the spirit of both Zonnon and BEL are welcomed and, 
if you choose, may be considered for inclusion in future releases.

--------------------------------------------------------------------------------

As of version 3.0, BEL uses NPlot to provide a graphical plotting capability.
NPlot is an open source .NET/mono software released under the following terms
and conditions.

NPlot home page:
   http://netcontrols.org/nplot/wiki/.
It is released under the terms of a 3-clause-BSD license.  Specifically:
   NPlot - A charting library for .NET
   Copyright (C) 2003-2006 Matt Howlett and others.  All rights reserved.
Redistribution and use in source and binary forms, with or without modifi-
cation, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of NPlot nor the names of its contributors may be used
   to endorse or promote products derived from this software without
   specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRI-
BUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

--------------------------------------------------------------------------------

If you find the BEL software suite to be useful, the author would like to hear
form you.

AUTHOR

Prof. Alan D. Freed
Clifford H. Spicer Chair in Engineering
Saginaw Valley State University
202 Pioneer Hall, 7400 Bay Road
University Center, MI 48710

Tel:  1-989-964-2288
Fax:  1-989-964-2717
Work: adfreed@svsu.edu
Home: alan.d.freed@gmail.com


DIRECTORY STRUCTURE

I started using the command-line compiler (to which I have since returned). 
As an aid I have written highlighting schemes for the Linux KDE Kate editor 
and the multi-platform jEdit editor, both of which are open source editors.
At present, I prefer the jEdit editor in conduction with the command-line 
compiler as my development platform.  jEdit requires a Java VM to be installed.
For both adaptations I use the Mono version of the Zonnon compiler, even when
running on a Windows machine.

In the past, I have used the Zonnon Builder and Eclipse IDEs to teach from.
Zonnon Builder comes with the Zonnon compiler for use on a Windows machine.
Eclipse is an open-source Java IDE developed by Oracle for which ETH-Z has 
written a Zonnon plugin for, which comes with their mono distribution of the
compiler.  I started teaching with Zonnon Builder, which is a very nice IDE. 
The switch was made from Zonnon Builder to Eclipse because Eclipse runs on
all three major platforms that I use: Mac OSX, Linux, and Windows.  This way 
my students were learning an IDE product that could be used on virtually any 
machine for programmatic development in almost any programming language. On each
of these platforms it requires (in addition to Eclipse, which can be downloaded
from:  http://www.eclipse.org) a Java RTE and Mono (even on a Windows machine).
Eclipse prefers you place your programming projects in a directory structure
~/workspace, where ~/ is your root directory.  It is strongly advised that this
root directory path be free from spaces of any kind, as this has been known to
confuse the Eclipse environment.  (Do not install Eclipse under the Windows 
default directory of C:\Program Files, where there is a space between Program 
and Files.)

Documentation is found in the <bel>/docs directory of the zip file, including
this file.  Also located there are several BEL user guides (one for each BEL
library), the GNU LGPL license, two highlighting files (called zonnon.xml) 
that are located in the zonnonEditorHighlighters directory: one for Kate 
Kate and the other for jEdit.  Zonnon.xml is now distributed with Kate, so 
most likely it is already installed on your Linux box, but the zonnon.xml 
file for jEdit has not been turned over the jEdit developers yet as it is 
still being developed.

THE BEL LIBRARY CONTAINS:

::Core Library::

   Bel.Entity
   Bel.Object
   Bel.Typeset
   Bel.Version
   Bel.Log
   Bel.DataFiles
   Bel.TextFiles
   Bel.Types
   Bel.Series
   Bel.Math
   Bel.Curves
   Bel.Plots
   Bel.Queue
   Bel.Stack
   Bel.Keys
   Bel.List
   Bel.Tree

::Genetic Algorithm::

   BelGAlg.Statistics
   BelGAlg.Genes
   BelGAlg.Chromosomes
   BelGAlg.Genome
   BelGAlg.Creatures
   BelGAlg.Colonies
   BelGAlg.GeneticAlgorithm

::Advanced Math Functions/Procedures::

   BelMath.Interpolations
   BelMath.Distributions
   BelMath.Regression
   BelMath.Roots
   BelMath.Derivatives
   BelMath.Integrals
   BelMath.RungeKutta
   BelMath.Irks
   BelMath.ConvolutionIntegrals


COMPILING USING THE COMMAND LINE (my preferred way to compile)

My directory structure typically looks something like:

MyWorkSpace
   bin
      System.Drawing.dll, NPlot.dll, NPlot.xml, Zonnon.RTL.dll, and
      (my projects libraries and executables)
   Zonnon
      compiler
         (the six files that come from ETH-Z under directory: compiler)
   Version###
      Bel
         (source code)
      BelGAlg
         (source code)
      BelMath
         (source code)
      É  (my personal directories)
      docs
   testCode
      (code written to test/validate my software)
   É

NOTICE:  You will need to place a copy the Zonnon runtime library Zonnon.RTL.dll
in your bin directory in order for your programs to run.

To compile, I locate myself in the directory that contains the source code that I'm 
about to compile.  From there I execute the command-line compiler from a terminal
window.  To compile the code on my machine (a Mac) I run commands that look like:

From the Bel directory:
myPrompt$ mono ../../Zonnon/compiler/zc.exe version.znn log.znn dataFiles.znn textFiles.znn types.znn objects.znn curves.znn plots.znn math.znn keys.znn lists.znn queues.znn stacks.znn series.znn trees.znn /ref:../../bin/System.Drawing.dll /ref:../../bin/NPlot.dll /out:../../bin/Bel

Note: I soft link System.Drawing.dll to its default location, which on a Mac is at:
/Library/Frameworks/Mono.framework/Versions/Current/lib/mono/2.0/System.Drawing.dll
NPlot.dll is third-party software that I distribute with Bel (see license info).

From the BelGAlg directory:
myPrompt$ mono ../../Zonnon/compiler/zc.exe statistics.znn genes.znn chromosomes.znn genome.znn creatures.znn colonies.znn geneticAlg.znn /ref:../../bin/Bel.dll /out:../../bin/BelGAlg

From the BelMath directory:
myPrompt$ mono ../../Zonnon/compiler/zc.exe interpolations.znn distributions.znn randomVariates.znn roots.znn derivatives.znn integrals.znn convolutionIntegrals.znn rungeKutta.znn regression.znn irks.znn /ref:../../bin/Bel.dll /out:../../bin/BelMath


INTERPRETING ERROR MESSAGES

If you get an error message, e.g., 
myPrompt$  mono ../../../Zonnon/compiler/zc.exe test.znn main.znn /ref:../../../bin/Bel.dll /entry:main /out:../../../bin/test
Zonnon Compiler, Version 1.3.0.0 of Tuesday, September 11, 2012, 10:22:24 PM
(c) 2003-2009 ETH Zurich
2019: /Users/adfreed/WorkspaceMac/testCode/Bel/testPlot/test.znn(7,7): Cannot import 'NPlot'. If it is an external name check if you have referenced the library
2021: /Users/adfreed/WorkspaceMac/testCode/Bel/testPlot/test.znn(46,22): Name 'LinePlot' does not denote a type as expected. 
2070: /Users/adfreed/WorkspaceMac/testCode/Bel/testPlot/test.znn(46,7): Entity 'NP.LinePlot;' is of an unknown type
É
Go to the first line in the error stream and resolve it first, then recompile.  
Often, but not always, the errors that follow are an artifact of previous errors, 
i.e., the compiler can get confused from prior bugs in your code.  To determine 
what is wrong, read the error message!  The error code (e.g., 2019) belongs to 
the compiler.  I'm not sure where to get their descriptions in English but, 
no matter, I have never needed to know what they say.  The next bit of information
tells you the file where the error occurred (in this case test.znn) followed by
(line number, column number) so you can go to your editor and pinpoint the 
location where the compiler says your error is.  Note, I have often found that 
the error may reside in the prior line, as will be the case if you forget to 
finalize that line with a semicolon.  Following this is a short statement of 
what the compiler is complaining about, in this case, it cannot import NPlot.  
Obviously not, I forgot to reference it, so correcting my error

myPrompt$  mono ../../../Zonnon/compiler/zc.exe test.znn main.znn /ref:../../../bin/Bel.dll /ref:../../../bin/NPlot.dll /entry:main /out:../../../bin/test

gives the desired compilation output of

myPrompt$  mono ../../../Zonnon/compiler/zc.exe test.znn main.znn /ref:../../../bin/Bel.dll /ref:../../../bin/NPlot.dll /entry:main /out:../../../bin/test
Zonnon Compiler, Version 1.3.0.0 of Tuesday, September 11, 2012, 10:22:24 PM
(c) 2003-2009 ETH Zurich
Compilation completed successfully

One can now move to the bin directory and execute your program, in this case as

myPrompt$ mono test.exe


COMPILE USING ECLIPSE

To compile these files as a .NET library, i.e., as a .DLL file using Eclipse,
copy the <source> directory from the Bel<version>.zip file into your Eclipse
<workspace> directory.  Rename the directory to whatever you choose; as for me,
I name it Bel_<versionNumber> where <versionNumber> is the version number of 
the software release, e.g., 3.3.  Now, open Eclipse, go to the menu item 
"File/New/Project.../Zonnon/Zonnon Project" to create the project.  When you
click on the "next" button a dialog will open where you can enter the name of
the directory you just created.  If you replicate it faithfully, when you go
into the project panel in Eclipse (which is on the left) and double click on
your project name, all the files will be automatically imported.  This is the
easiest way to create a project like this.  

Now that you've created your new project, you need to tell the compiler what 
it will be compiling, and what libraries it will require.  This is done by 
highlighting your project in the project panel, and then following the menu 
path "Project/Properties/Zonnon Compiler" which will open a dialog window.  
In this dialog, you will need to name the output file.  I use the path:  
./../../bin/Bel_<versionNumber>.dll  which is a relative reference (that is why 
it starts with ./).  From here I choose to go back up the directory path twice 
(which is done by  ../../  like from a DOS or Unix command line), and then place 
the to-be-created .DLL into my  ~/workspace/bin  directory.  Notice that 
Eclipse will also create a bin directory under each project.  This is normal.
It is used by Eclipse, so do not delete them.  I choose to place all if my 
.DLL and .EXE files in a common bin directory, but this is not necessary.  
In the case of compiling a library like Bel, the "Startup module" entry needs 
to be empty, and the output type selected should be "Class Library".  In the 
case of Bel, one also needs to reference two external libraries, so press the 
"Add Reference" button and enter "System.Drawing.dll" into the text window. 
Select this and press the "Add Reference" button one more time, this time 
entering "./../../bin/NPlot.dll", which assumes that you have placed it into
your bin directory (NPlot.zip is supplied with Bel).  NPlot comes with a 
variety of prebuilt libraries.  Because Eclipse runs Mono, choose the one
for Mono running framework 2.  Press the "Apply" button.  At this juncture,
Eclipse will likely compile Bel for you, unless you toggled off the "Build
Automatically" menu button in the "Project" Menu.  If Eclipse did not build
the library automatically for you, about four or five icons from the left in 
the menu bar is an icon with 0101 in it.  Press it to compile your library.

If errors occur, they will be reported via the Console window at the bottom
of the Eclipse IDE.  Red squares will locate the error in the editor window
of the affected software module.  A debugger has not been written for Zonnon,
so this capability of Eclipse cannot be used; nevertheless, I find the error
messages that are reported into the Console window, with their associated 
links to where the errors occur, to be sufficient for my own debugging
purposes.  


DEBUGGER

If someone is interested in writing a debugger for Zonnon, please contact 
this author and he will put you in touch with the Zonnon compiler developers.


::Test Code::

The test code used to validate BEL is included in the subdirectory <bel>/test.
These examples are provided so that you can learn how to program, both in
Zonnon and with the BEL package.  I hope you find them useful.  They include:

   i)     testLog.znn
   ii)    testFiles.znn
   iii)   testPlots.znn
   iv)    testNumbers.znn
   v)     testMath.znn
   vi)    testArrays.znn
   vii)   testMatrices.znn
   viii)  testLinearAlgebra.znn
   ix)    testDataStructures.znn
   x)     testRoots.znn
   xi)    testInterpolation.znn
   xii)   testDistributions.znn
   xiii)  testRegression.znn
   xiv)   testDerivatives.znn
   xv)    testIntegrals.znn
   xvi)   testRungeKutta.znn

To compile an application, you will need to first compile BEL as a library.
Then you can create a project for any of the above test cases in which you
load BEL.dll (or something similarly named) as an external library.


HISTORY OF CHANGES

Version 4.1
Sep 2012 - The plotting capability has been enhanced so that the resulting graphs 
           can be tailored more to the individuals liking.  Some minor bugs were 
           fixed.

Version 4.0
Jan 2012 - This was another major restructuring of the code to prepare if for the
           next release of the compiler.  In that release, all constants will have
           their associated default numeric type.  This is huge, although on the
           surface it seems trivial.  Because of this change, having to deal with
           all possible .NET numeric types has reduced down to two: integer and
           real.  Type explosion is no longer a dominant driving force behind the
           design of this code.  Types Number, Array and Matrix were able to be
           eliminated, because it is now possible to rely on the built-in types.
           This change effects virtually all modules in BEL, and in many of their
           designs and interfaces; hence, the need to start a new version sequence.

Version 3.3
Oct 2011 - There was a major restructuring of the three BEL types: Numbers,
           Arrays, and Matrices.  These were largely internal changes done for
           the purpose of being able to take advantage of runtime efficiencies
           that are to be part of the next compiler release.  Minimal interface
           changes took place.  The modules were also collected into a different
           packaging of DLLs.  This was done to simplify things for my students.

Version 3.2
Jul 2011 - Rely on .NET and Mono implementations of IEEE 754 standard for
           floating point arithmetic to handle infinities and NaNs properly in
           an effort to speed up the basic arithmetic operations, viz., removed
           all precondition testing to address all possible special cases in
           software.  Now relying on this to be done in hardware.

Version 3.0
Jun 2011 - Reworked numbers.znn, in particular, made NumberType a wrapper for
           the underlying type that is a Bel.Numbers.Number.  This will make it
           much easier to implement complex as the base type at a future time.
           Cleaned up the code here and there in many of the modules.
Apr 2011 - Added modules curves.znn and plots.znn to BEL so that 2D graphical
           plots can be made.  These are created as bitmaps and written to
           file where the user can access them.  They are placed in the iofiles
           subdirectory below the directory where your program was executed.

Version 2.3 released, January 2011
Jan 2011 - A simplified stack tracing method used when logging an error message.
           Procedure LnDeterminant was added to Bel.LinearAlgebra.
Dec 2010 - Corrected R^2 statistic being reported back in Bel.MATH.Regression
           procedures.  It was reporting the coefficient of correction, i.e.,
           R instead of R^2, which is the coefficient of determination.
           Relaxed error tolerance on Gauss-Kronrod integration a bit to make
           it more robust.  Moved modules Irks.znn and convolution.znn to the
           Bel.SS library, as this is what they were written for.  Got rid of
           secondary name spaces of IO, DATA, MF and MATH so all modules in
           BEL now have the name space Bel.<moduleName>.  Applications, like
           Bel.GA, will still have secondary name spaces assigned to them.

Version 2.1 released, November 2010
Nov 2010 - Fixed roundoff error (occurred in .NET, not Mono) causing overflow
           from type conversion in procedure GetKey in Bel.MATH.Integrators.
Oct 2010 - Ported over my solver for convolution integrals from CAPO.  Rewrote
           module Bel.MATH.Series replacing definition Bel.EvaluateSeries with
           the procedure type Bel.MATH.GetCoef.  This forced some internal
           changes in Bel.MATH.Functions, but not in its interface.  The inter-
           face of Bel.MATH.Series was significantly altered; it is simpler now.
           This was necessary so that series could be evaluated outside the
           core Bel.dll.  Something in the old Bel.EvaluateSeries interface
           prohibited this - not sure what, though.  Added methods IsInteger
           and GetInteger to a Bel.MF.Numbers.Number.

Version 2.0 released, August 2010
Apr 2010 - Put the applications into seperate libraries to break apart BEL into
           more logical components.  Rewrote the base definitions for objects.
           Consequently, BEL-2.0 may not be backward compatible to code you've
           written for earlier versions of BEL.  Definition Bel.Object has
           changed, and definitions Bel.Datum and Bel.Typeset are new.  These
           changes were done in an effort to improve their teaching ability,
           as I use this library in my classes.

Version 1.3 internal.
Apr 2010 - Fixed an error in RungeKutta.Integrator.Solve so that it correctly
           handles integration in the negative direction, accomplished by
           assigning a negative valued 'h' to RungeKutta.Integrator.Initialize.
           Added an Integrate procedure for user control over stepsize h.
Mar 2010 - Fixed a bug in Colony.Update to correctly handle the elite
           individual in the ancestory list, and in GeneticAlgorithm.Optimize
           in its writing of the parameters to the various files.  A Parse
           method was added to type Bel.MF.Arrays.Array.
Feb 2010 - Fixed a bug in an index-overflow error check in Array.SetSubArray
           and in Matrix.SetSubMatrix.  Introduced Bel.MATH.Irks to solve
           stiff ODEs using intrinsic Runge-Kutta stable integrators, and
           altered Runge-Kutta interface slightly to make it compatible with
           Irks interface.  Introduced a module for least-squares regression.
Jan 2010 - The Runge-Kutta integrator was reworked to introduce Richardson
           extrapolation per John Butcher's suggestion.  The PI.4.2 controller
           of Soderlind was implemented therein for adaptive step control.
           Butcher's explicit and implicit IRKS methods 555b and 556a were
           coded into IRKS.  Removed newtonRaphson.znn from the distribution.
           The genetic algorithm, Bel.GA, is a very effective replacement for
           it.  Bel.MATH.Roots was added for finding the roots of a function.

Version 1.2 released.
Jan 2010 - Reworked the hypoelastic constitutive model.
Dec 2009 - Genetic algorithm fitness function changed to R^2 statistic, etc.
           Confidence regions added to the genetic algorithm's capability.
           Added method ToString to type Bel.PF.Scalars.Scalar.
Nov 2009 - Created BelExtra in an effort to pair down the number of files in
           Bel.  Moved out of Bel were the more sophisticated material models.
Oct 2009 - Incorporated Nina's changes into Bel.MF.Arrays, Bel.MF.Matrices,
           Bel.PF.Vectors2, Bel.PF.Tensors2 and Bel.PF.QuadTensors2 taking
           advantage of the new math extensions added to the Zonnon language.
Sep 2009 - Revised definition for computing the residuals in the statistics
           module for the GA.  Rewrote the hypoelastic modules to reflect a
           change in its structure suggested by a reviewer of my paper.
Aug 2009 - Corrected an error in the traction-to-stress mapping in Kinetics.
           Revised the genetic algorithm according to Goldberg's 2nd book.
Jul 2009 - Corrected a couple of errors in the hypoelastic formulation and
           added an discussion about hypoelasticity to the users manual.
           Created a test case for hypoelasticity with the genetic algorithm.
Jun 2009 - Prof. Blake Johnson at SVSU finished the art work for the BEL icon.
           Wrote modules for isotropic, anisotropic and composite, isochoric,
           hypoelastic, constitutive equations for tissue mechanics modeling,
           and added module Physical2 as a bottleneck for constitutive models.
May 2009 - The major rewrite of the core, that became vs 1.2, was completed.
           Added the genetic algorithm suite of modules to the distribution.
           Wrote the first documentation for the library.
Apr 2009 - A major rewrite of the library was begun.  Many changes were made.
           The purpose was to unify the various DLLs into a single DLL because
           of some meta programming issues that arise when using multiple DLLs
           using the current state of the compiler.  Roman Mitin is aware of
           these issues, and they will be fixed in a future compiler release.

Version 1.1 - an internal release only.
Mar 2009 - Unified the interpolation algorithms into a single procedure call
           in Math.Interpolations.
Feb 2009 - Changed licensing from GNU General Public License to the GNU Lesser
           General Public License, the latter being designed for libraries.
           Cleaned up some modules due to an improvement in the compiler that
           checks for declared but unused variables.  Added 'Bel' prefix to all
           module and dll names to avoid potential name conflicts.  Removed
           material models from CCM; it was felt that CCM should pertain to
           just the general linear algebra parts of continuum tensor analysis.
           Added Eigenvalues, Eigenvectors and SpectralDecomposition procedures
           to CCM Tensors2 module.
Dec 2008 - Rewrote Math.Derivatives and Math.Integrals to allow for both fixed
           step-size algorithms, and ones that internally adjust the step size
           and use Richardson extrapolation to achieve a predefined level of
           precision, i.e., course (roughly 8 significant figures), medium
           (roughly 11 significant figures) and fine (roughly 14 significant
           figures).
Nov 2008 - Fixed bug and introduced enumerated type in Math.Distributions.znn
           Restructured the error codes in Core.Log so that it would be easier
           for a user/programmer to locate an error message or add a new one.
           This required all other files to have their log calls remapped.
           Extracted NewtonRaphson from Optimization.Static and moved it to
           Math.NewtonRaphson, allowing the Opt library to be removed
Oct 2008 - fixed a bug in the Core.Math.ArcTan2 function
           Change type conversion procedures from private to public in module
           Core.Numbers.  This was done so that the definitions of overloaded
           operators in other modules could be rewritten in terms of proper
           procedure calls only, i.e., free from other overloaded operators.
           This avoids nesting of overloaded operators, and should improve on
           overall efficiency.  This forced changes in modules: Core.Arrays,
           Core.Matrices, Ccm.Units, Ccm.Scalars, Ccm.Vectors, Ccm.Tensors and
           Ccm.QuadTensors.
           Renamed modules biVectors.znn -> vectors.znn, biTensors.znn ->
           tensors.znn and biQuadTensors.znn -> quadTensors.znn in CCM.
           CCM is now considered to be just for membrane analysis.
           Added modules kinematics.znn, hypoelastic.znn, displacementBVP.znn
           and tractionBVP.znn to the CCM library.
Sep 2008 - Added some new error codes to Core.Log.
           Overloaded arithmetic '+' and '-' and boolean '=', '#', '>', '>=',
           '<', and '<=' operators for Ccm.Scalars.Scalar introduced to handle
           operations between dimensionless scalars and the core number types.
           Allow negative tolerance in Static and Dynamic optimizers to handle
           case where statistics of fit can be got without updating parameters.
           Bug fixes in the Static and Dynamic optimizers.

Version 1.0 released.
Aug 2008 - Two additional libraries have been added: Ccm.dll and Opt.dll. These
           are more on the application side.
May 2008 - The first three BEL libraries were added: Core.dll, Data.dll and
           Math.dll.  These were ports from the Oberon package CAPO.