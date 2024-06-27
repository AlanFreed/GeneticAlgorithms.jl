(* *****************************************************************************
  genAlg Author         - Alan D. Freed
  genAlg License        - GNU Lesser General Public License, vs 3 or later
  genAlg Copyright      - (c) Alan D. Freed 2014
  genAlg Pascal Version - 1.2
--------------------------------------------------------------------------------
  Use a genetic algorithm for the purpose of parameter estimation.
--------------------------------------------------------------------------------
  This application uses three external libraries, all written by the author:
     'Arrays' provides a collection of array type and procedures to manage them
     'rng'    provides a random number generator and statistical functions.
     'GenAlg' provides a genetic algorithm for use in parameter estimation.
--------------------------------------------------------------------------------
   This is one of three parts to the template.  Here your model is defined.
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

unit gaModel;

interface

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    Classes, SysUtils, gaData;
    
  // Exported constants for the genetic algorithm pertinent to this application

  const
    dimP  = 2;          // number of parameters, both fixed and varied
    dimPV = 2;          // number of parameters that are allowed to vary
    
  // Array types defined in 'Arrays' and used by the 'GenAlg' library
  
  type
    TIVector = array of Integer;
    TSVector = array of String;
    
  // Additional array type needed by 'GenAlg' library
  
  type
    TBVector = array of Boolean;

  // A model is implemented for parameter estimation via an instance of type

  type
    TModel = procedure (experiment           : Integer;      // input
                        var modelParameters  : TRVector;     // input
                        var controlData      : TRMatrix;     // input
                        var responseData     : TRMatrix);    // output
{ where
    experiment       : used to determine which experiment to get response data
                       an interger between 1 and dimE
    modelParameters  : all model parameters (both fixed and adjustable)
                       indices range over [1..dimP]
    controlData      : a sequence of experimentally controlled values
                       [controlVariable][valueAtASamplingInstant]
                       indices range over [1..dimNC[e]][1..dimNS[e]]
    responseData     : predicted responses for these controls & parameters
                       [responseVariable][valueAtASamplingInstant]
                       indices range over [1..dimNR[e]][1..dimNS[e]]           }
                       
  // Exported procedures and functions used by the driver
                       
  procedure MyModel (experiment           : Integer;      // input
                     var modelParameters  : TRVector;     // input
                     var controlData      : TRMatrix;     // input
                     var responseData     : TRMatrix);    // output
  { an instance of type TModel }

  // Associated functions that are sent to the genetic algorithm
  
  function VaryParameters : TBVector;
  { specify which parameters are to be varied (TRUE) and which are not (FALSE) }
  
  // The remaining procedure provide information for the varied parameters only
  
  function NameParameters : TSVector;
  { assign names to the varied parameters that will appear in the report }
  
  function FixedParameters : TRVector;  
  { assign values to the fixed parameters - if all vary return nil }
    
  function AlienParameters : TRVector;  
  { best guess at what the varied parameters might be }
    
  function MaximumParameters : TRVector;  
  { upper boundary of the search domain for all of the varied parameters }
    
  function MinimumParameters : TRVector;  
  { lower boundary of the search domain for all of the varied parameters }
  
implementation

  var
    emptyBVector : TBVector;
    
  // Functions imported from the 'Arrays' library.  Others are available.
  
  // Vectors index from [1..len].

  function NewRVector (len : Integer) : TRVector; external 'Arrays';
  // Creates a real vector indexing from 1 to len whose elements are all 0.0.

  function NewSVector (len : Integer) : TSVector; external 'Arrays';
  // Creates a string vector indexing from 1 to len whose elements are all ''.

  // Matrices index from [1..rows][1..cols]

  function NewRMatrix (rows, columns : Integer) : TRMatrix; external 'Arrays';
  // Creates a real matrix indexing as [1..rows][1..columns]
  // whose elements are assigned values of 0.0.

  // Boolean arrays are not supplied by library 'Arrays' so we define them here.

  function NewBVector (len : Integer) : TBVector;
    var
      i : Integer;
      v : TBVector;
  begin
    // Allocate vector
    try
      SetLength(v, Succ(len))
    except on E : EOutOfMemory do
      begin
        WriteLn('Memory error. Details: ' + E.ClassName + '/' + E.Message);
        Exit(emptyBVector);
      end
    end;
    // Initialize vector
    for i := 0 to len do
      v[i] := False;      // location [0] is not to be used in applications
    NewBVector := v
  end;

  // Create the model to be used to fit the data against

  // Model parameters to be solved for are:
  //   b  is the y-intercept of the line
  //   m  is the slope of the line
  // which describe  y = mx + b, the model sent to the genetic algorithm.

  procedure MyModel (experiment           : Integer;      // input
                     var modelParameters  : TRVector;     // input
                     var controlData      : TRMatrix;     // input
                     var responseData     : TRMatrix);    // output
    var
      b, m : Real;
      s    : Integer;
  begin
    // assign model pararmeters to global variables
    m := modelParameters[1];    // slope
    b := modelParameters[2];    // y intercept
    // create the response data predicted by the model
    responseData := NewRMatrix(dimNR[experiment], dimS);
    for s := 1 to dimS do
      responseData[dimNR[experiment]][s]
        := m * controlData[dimNC[experiment]][s] + b
  end;
  
  // user-defined procedures required in every GenAlg application
  
  // the length of VaryParameters is for all parameter, both varied and fixed

  function VaryParameters : TBVector;
    var
      vary : TBVector;
  begin
    vary := NewBVector(dimP);
    vary[1] := True;    // vary m
    vary[2] := True;    // vary b
    VaryParameters := vary
  end;
  
  // all remaining functions return vectors of the varied parameter length

  function NameParameters : TSVector;
    var
      names : TSVector;
  begin
    names := NewSVector(dimPV);
    // assign names to both the fixed and varied parameters
    names[1] := 'slope m';
    names[2] := 'intercept b';
    NameParameters := names
  end;

  function FixedParameters : TRVector;
    var
      fixed : TRVector;
  begin
    fixed := nil;     // all parameters will vary, none are fixed
    FixedParameters := fixed
  end;

  function AlienParameters : TRVector;
    var
      alien : TRVector;
  begin
    // this is a guess at what one expects the varied parameters should be
    alien := NewRVector(dimPV);
    // because this is a demo, choose not so good initial guesses
    alien[1] := 1.2 * mm;
    alien[2] := 0.8 * bb;
    AlienParameters := alien
  end;

  function MaximumParameters : TRVector;
    var
      maxP : TRVector;
  begin
    // the greatest (or most positive) parameter values to be considered
    maxP := NewRVector(dimPV);
    if mm > 0.0 then            // greatest value for m
      maxP[1] := 5.0 * mm
    else if mm = 0.0 then
      maxP[1] := 1.0
    else
      maxP[1] := mm / 5.0;
    if bb > 0.0 then            // greatest value for b
      maxP[2] := 5.0 * bb
    else if bb = 0.0 then
      maxP[2] := 1.0
    else
      maxP[2] := bb / 5.0;
    MaximumParameters := maxP
  end;

  function MinimumParameters : TRVector;
    var
      minP : TRVector;
  begin
    // the smallest (or most negative) parameter values to be considered
    minP := NewRVector(dimPV);
    if mm > 0.0 then            // smallest value for m
      minP[1] := mm / 5.0
    else if mm = 0.0 then
      minP[1] := -1.0
    else
      minP[1] := 5.0 * mm;
    if bb > 0.0 then            // smallest value for b
      minP[2] := bb / 5.0
    else if bb = 0.0 then
      minP[2] := -1.0
    else
      minP[2] := 5.0 * bb;
    MinimumParameters := minP
  end;
    
begin

  SetLength(emptyBVector, 0);
  
  // initialize additional variables used in this application   

end.
