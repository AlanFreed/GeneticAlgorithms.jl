(* *****************************************************************************
  genAlg Author         - Alan D. Freed
  genAlg License        - GNU Lesser General Public License, vs 3 or later
  genAlg Copyright      - (c) Alan D. Freed 2014
  genAlg Pascal Version - 1.3.3
--------------------------------------------------------------------------------
  Use a genetic algorithm for the purpose of parameter estimation.
--------------------------------------------------------------------------------
  This application uses three external libraries, all written by the author:
     'Arrays' provides a collection of array type and procedures to manage them
     'rng'    provides a random number generator and statistical functions.
     'GenAlg' provides a genetic algorithm for use in parameter estimation.
--------------------------------------------------------------------------------
   This is one of five parts to the template.  Here the data types are entered.
--------------------------------------------------------------------------------
   References
      Goldberg, D.E., Genetic Algorithms in Search, Optimization, and Machine
         Learning, Addison-Wesley, Boston, 1989.
      Goldberg, D.E., The Design of Innovation: Lessons learned from and for
         competent genetic algorithms.  In: Genetic algorithms and evolutionary
         computation, Vol. 7, Klewer, Boston, 2002.
      Johnson, N.L., "Systems of Frequency Curves Generated by Methods of
         Translation," Biometrika, Vol. 36 (1949), 149-176.       
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
    
  // Array types defined in 'Arrays' and used by the 'GenAlg' library

  type
    TBVector = array of Boolean;
    TIVector = array of Integer;
    TRVector = array of Real;
    TSVector = array of String;
    TRMatrix = array of TRVector;
    TRData   = array of TRMatrix;

  // Functions used to create and manage the array types of 'Arrays'.
  
  function NewBVector (len : Integer) : TBVector;
  // Creates a boolean vector indexing from 1 to len with elements of FALSE.
  
  function NewIVector (len : Integer) : TIVector; 
  // Creates an integer vector indexing from 1 to len whose elements are all 0.

  function NewRVector (len : Integer) : TRVector;
  // Creates a real vector indexing from 1 to len whose elements are all 0.0.

  function NewSVector (len : Integer) : TSVector; 
  // Creates a string vector indexing from 1 to len whose elements are all ''.
  
  function LenBVector (var v : TBVector) : Integer; 
  // Returns the length of boolean vector v.

  function LenIVector (var v : TIVector) : Integer; 
  // Returns the length of integer vector v.

  function LenRVector (var v : TRVector) : Integer;
  // Returns the length of real vector v.
  
  function LenSVector (var v : TSVector) : Integer; 
  // Returns the length of string vector v.

  function NewRMatrix (rows, columns : Integer) : TRMatrix; 
  // Creates a real matrix indexing as [1..rows][1..columns]
  // whose elements are assigned values of 0.0.

  function NewRData (experiments         : Integer;
                     var variablesPerExp : TIVector;
                     var samplesPerExp   : TIVector) : TRData;
  // Creates a new real-valued data structure indexing as
  // [1..experiments][1..variablesPerExp[experiment]][1..samplesPerExp[experiment]
  // whose elements are all assigned values of 0.0.
  
  // You may need other functions exported by 'Arrays'.  These are most common.
  
implementation
  
  function NewBVector (len : Integer) : TBVector; external 'Arrays';
  
  function NewIVector (len : Integer) : TIVector; external 'Arrays';

  function NewRVector (len : Integer) : TRVector; external 'Arrays';

  function NewSVector (len : Integer) : TSVector; external 'Arrays';
  
  function LenBVector (var v : TBVector) : Integer; external 'Arrays';

  function LenIVector (var v : TIVector) : Integer; external 'Arrays';

  function LenRVector (var v : TRVector) : Integer; external 'Arrays';
  
  function LenSVector (var v : TSVector) : Integer; external 'Arrays';

  function NewRMatrix (rows, columns : Integer) : TRMatrix; external 'Arrays';

  function NewRData (experiments : Integer; var variablesPerExp : TIVector;
                     var samplesPerExp : TIVector) : TRData; external 'Arrays';

begin
end.
