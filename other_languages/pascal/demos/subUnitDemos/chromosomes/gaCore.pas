(* *****************************************************************************
   Author         - Alan D. Freed, Ph.D
   License        - GNU Lesser General Public License, vs. 3 or later
   Copyright      - (c) Alan D. Freed 2014-2015
   Pascal Version - 1.3.3
--------------------------------------------------------------------------------
   GenAlg is a free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the Free
   Software Foundation, either version 3 of the License, or (at your option)
   any later version.

   GenAlg is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
   more details.

   You should have received a copy of the GUN Lesser General Public License
   along with GenAlg.  If not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------------------
   Port began on  - March    18, 2014
   Last modified  - November 21, 2015
--------------------------------------------------------------------------------
   Version 1.2
   Adopted 'Arrays' for creating and managing vector arrays.  
   This library calls librarys 'Arrays' and 'rng' written by the author.
--------------------------------------------------------------------------------
   This unit establishes the global data used by the genetic algorithm.
***************************************************************************** *)

// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

unit gaCore;

interface

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    Classes, SysUtils;

// -----------------------------------------------------------------------------

  // constants required by the 'rng' library, which is heavily used by genalg
  const
    SB = 1;      // denotes a Johnson SB (bounded)   distribution
    SL = 2;      // denotes a Johnson SL (lognormal) distribution
    SN = 3;      // denotes a Johnson SN (normal)    distribution
    SU = 4;      // denotes a Johnson SU (unbounded) distribution
    ST = 5;      // denotes a Johnson SB (bounded)   distribution
                 //   near the boundary of inadmissibility

  // constants used by the genetic algorithm

  // machine-dependent constants for 64-bit floating point real numbers
  const
    machEp = 2.220446049250313E-16;     // Floating point precision: 2^(-52)
    maxNum = 1.797693134862315E+308;    // Max. floating point number: 2^128
    minNum = 2.225073858507202E-308;    // Min. floating point number: 2^(-126)
    maxLn  =  709.7827128933840;        // Max. argument for Exp = Ln(maxNum)
    minLn  = -708.3964185322641;        // Min. argument for Exp = Ln(minNum)

  // other constants required by the genetic algorithm
  const
    pi                           = 3.141592653589793;
    evenOdds                     = 0.5;
    meanProbabilityOfMutation    = 0.05;
    meanProbabilityOfCrossover   = 0.8;
    standardDeviationOfMutation  = 0.015;
    standardDeviationOfCrossover = 0.05;
    tolerance                    = 0.0001;
    dimensionOfSchemata          = 7;
    maxSignificantFigures        = 7;
    maxSteps                     = 250;

  // basic array types imported from and managed by the 'Arrays' library

  type
    // vector arrays
    TBVector  = array of Boolean;
    TIVector  = array of Integer;
    TRVector  = array of Real;
    TSVector  = array of String;
    // rectangular data arrays
    TIMatrix  = array of TIVector;
    TRMatrix  = array of TRVector;
    // staggered data arrays
    TRData    = array of TRMatrix;
    // vector function
    TVectorFn = function (var v : TRVector) : TRVector;
    
  // functions and procedures imported from the 'Arrays' library
  
  function NewBVector (len : Integer) : TBVector;
  // Creates a boolean vector indexing from 1 to len with elements of FALSE.

  function NewIVector (len : Integer) : TIVector;
  // Creates an integer vector indexing from 1 to len whose elements are all 0.

  function NewRVector (len : Integer) : TRVector;
  // Creates a real vector indexing from 1 to len whose elements are all 0.0.

  function NewSVector (len : Integer) : TSVector; 
  // Creates a string vector indexing from 1 to len whose elements are all ''.

  function NewRMatrix (rows, columns : Integer) : TRMatrix; 
  // Creates a real matrix indexing as [1..rows][1..columns]
  // whose elements are assigned values of 0.0.

  function NewRData (experiments         : Integer;
                     var variablesPerExp : TIVector;
                     var samplesPerExp   : TIVector) : TRData;
  // Creates a new real-valued data structure indexing as
  // [1..experiments][1..variablesPerExp[experiment]][1..samplesPerExp[experiment]
  // whose elements are all assigned values of 0.0.
  
  function CopyBVector (var v : TBVector) : TBVector;
  // Creates a deep copy of boolean vector v.
  
  function CopyRVector (var v : TRVector) : TRVector;
  // Creates a deep copy of real vector v.
  
  function CopySVector (var v : TSVector) : TSVector;
  // Creates a deep copy of string vector v.
  
  function CopyRMatrix (var m : TRMatrix) : TRMatrix; 
  // Creates a deep copy of real matrix m.
  
  function CopyRData (var d : TRData) : TRData;
  // Creates a deep copy of the data structure d.
  
  procedure DelBVector (var v : TBVector); 
  // Deletes the information held by dynamic boolean vector v.
  
  procedure DelIVector (var v : TIVector); 
  // Deletes the information held by dynamic integer vector v.
  
  procedure DelRVector (var v : TRVector); 
  // Deletes the information held by dynamic real vector v.
  
  procedure DelIMatrix (var m : TIMatrix);
  // Deletes the information held by integer matrix m.
  
  procedure DelRMatrix (var m : TRMatrix); 
  // Deletes the information held by real matrix m.
  
  procedure DelRData (var d : TRData); 
  // Deletes all information in data structure d.

  function LenBVector (var v : TBVector) : Integer; 
  // Returns the length of boolean vector v.

  function LenRVector (var v : TRVector) : Integer; 
  // Returns the length of real vector v.

  
implementation

  function NewBVector (len : Integer) : TBVector; external 'Arrays';

  function NewIVector (len : Integer) : TIVector; external 'Arrays';

  function NewRVector (len : Integer) : TRVector; external 'Arrays';

  function NewSVector (len : Integer) : TSVector; external 'Arrays';

  function NewRMatrix (rows, columns : Integer) : TRMatrix; external 'Arrays';

  function NewRData (experiments         : Integer;
                     var variablesPerExp : TIVector;
                     var samplesPerExp   : TIVector) : TRData; external 'Arrays';
  
  function CopyBVector (var v : TBVector) : TBVector; external 'Arrays';
  
  function CopyRVector (var v : TRVector) : TRVector; external 'Arrays';
  
  function CopySVector (var v : TSVector) : TSVector; external 'Arrays';
  
  function CopyRMatrix (var m : TRMatrix) : TRMatrix; external 'Arrays';
  
  function CopyRData (var d : TRData) : TRData; external 'Arrays';
  
  procedure DelBVector (var v : TBVector); external 'Arrays';
  
  procedure DelIVector (var v : TIVector); external 'Arrays';
  
  procedure DelRVector (var v : TRVector); external 'Arrays';
  
  procedure DelIMatrix (var m : TIMatrix); external 'Arrays';
  
  procedure DelRMatrix (var m : TRMatrix); external 'Arrays';
  
  procedure DelRData (var d : TRData); external 'Arrays';
  
  function LenBVector (var v : TBVector) : Integer; external 'Arrays';

  function LenRVector (var v : TRVector) : Integer; external 'Arrays';
  
begin

end.
