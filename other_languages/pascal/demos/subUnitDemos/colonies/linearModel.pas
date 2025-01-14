// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

program linearModel;

// this example considers fitting random data to a straight line where the 
// random data are noise about a line at a specified standard deviation

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    Classes, SysUtils, gaCore, gaProbabilities, gaStatistics, gaGenes, 
    gaChromosomes, gaGenome, gaCreatures, gaSpecies, gaColonies;

  const
    nbrC = 1;           // number of controlled variables
    nbrE = 1;           // number of experiments
    nbrN = 100;         // number of experimental sample points in the data set
    nbrM = 2;           // number of model points created for drawing a curve
    nbrP = 2;           // number of parameters
    nbrR = 1;           // number of response variables
    xMin = 0.0;         // lower bound for the independent variable
    xMax = 10.0;        // upper bound for the independent variable
  
  type
    // vector types
    TIVector = array of Integer;
    TRVector = array of Real;
    TSVector = array of String;
    // matrix types
    TIMatrix = array of TIVector;
    TRMatrix = array of TRVector;
    TSMatrix = array of TSVector;
    // data type
    TRData   = array of TRMatrix;

  var
    a0    : Real = 0.0; // y-intercept for data from linear regression
    a1    : Real = 0.0; // slope of data gotten from linear regression
    bb    : Real;       // mean y-intercept of noisy experimental data
    mm    : Real;       // mean slope of noisy experimental data
    xData : TRVector;   // independent variables in the experimental data set
    yData : TRVector;   //  dependent  variables in the experimental data set

  // ---------------------------------------------------------------------------
  
  // Functions supplied by the 'Arrays' library

  function NewIVector (len : Integer) : TIVector; external 'Arrays';
  // Creates an integer vector indexing from 1 to len whose elements are all 0.
  
  function NewRVector (len : Integer) : TRVector; external 'Arrays';
  // Creates a real vector indexing from 1 to len whose elements are all 0.0.

  function NewSVector (len : Integer) : TSVector; external 'Arrays';
  // Creates a string vector indexing from 1 to len whose elements are all ''.

  function NewRData (experiments         : Integer;
                    var variablesPerExp : TIVector;
                    var samplesPerExp   : TIVector) : TRData; external 'Arrays';
  // Creates a new real-valued data structure indexing as
  // [1..experiments][1..variablesPerExp[experiment]][1..samplesPerExp[experiment]
  // whose elements are all assigned values of 0.0.
  
  procedure DelIVector (var v : TIVector); external 'Arrays';
  // Deletes the information held by dynamic integer vector v.
  
  procedure DelRVector (var v : TRVector); external 'Arrays';
  // Deletes the information held by dynamic real vector v.

  function NewRMatrix (rows, columns : Integer) : TRMatrix; external 'Arrays';
  // Creates a real matrix indexing as [1..rows][1..columns]
  // whose elements are assigned values of 0.0.
  
  procedure DelRMatrix (var m : TRMatrix); external 'Arrays';
  // Deletes the information held by real matrix m.

  // Function supplied by the 'rng' library

  function RandomNormal (mean, standardDeviation : Real) : Real; external 'rng';

  // ---------------------------------------------------------------------------

  // Create the data to be used by the genetic algorithm.

  procedure ReadInData;
  var
    detM, sd  : Real;
    n, p, q   : Integer;
    mInv      : TRMatrix;
    mMtx      : TRMatrix;
    rhs       : TRVector;
    summation : Real;
    y         : Real;
    zMtx      : TRMatrix;
  begin
    WriteLn;
    Write('What is the y-intercept of the line to be? ');
    ReadLn(bb);
    Write('What is the slope of the line to be?       ');
    ReadLn(mm);
    Write('What is the standard deviation for noise?  ');
    sd := 0.0;
    ReadLn(sd);
    for n := 1 to nbrN do begin
      xData[n] := xMin + (n - 1) * (xMax - xMin) / (nbrN - 1);
      y        := mm * xData[n] + bb;
      yData[n] := RandomNormal(y, sd)
    end;
    // perform a linear regression of these data
    zMtx := NewRMatrix(nbrN, nbrP);
    for n := 1 to nbrN do begin
      zMtx[n][1] := 1.0;
      zMtx[n][2] := xData[n]
    end;
    mMtx := NewRMatrix(nbrP, nbrP);
    for p := 1 to nbrP do
      for q := p to nbrP do begin
        summation := 0.0;
        for n := 1 to nbrN do
          summation := summation + zMtx[n][p] * zMtx[n][q];
        if p = q then
          mMtx[p][p] := summation
        else begin
          mMtx[p][q] := summation;
          mMtx[q][p] := summation
        end
      end;
    rhs := NewRVector(nbrP);
    for p := 1 to nbrP do begin
      summation := 0.0;
      for n := 1 to nbrN do
        summation := summation + zMtx[n][p] * yData[n];
      rhs[p] := summation
    end;
    mInv := NewRMatrix(nbrP, nbrP);
    detM := mMtx[1][1] * mMtx[2][2] - mMtx[1][2] * mMtx[2][1];
    mInv[1][1] :=  mMtx[2][2] / detM;
    mInv[1][2] := -mMtx[1][2] / detM;
    mInv[2][1] :=  mInv[1][2];
    mInv[2][2] :=  mMtx[1][1] / detM;
    a0 := mInv[1][1] * rhs[1] + mInv[1][2] * rhs[2];
    a1 := mInv[2][1] * rhs[1] + mInv[2][2] * rhs[2];
    // clean up
    DelRMatrix(mInv);
    DelRMatrix(mMtx);
    DelRVector(rhs);
    DelRMatrix(zMtx)
  end;

  // Reconfigure these experimentat data for use with the genetic algorithm.

  function GetXData : TRData;
  var
    n          : Integer;
    cVec, nVec : TIVector;
    xMtxData   : TRData;
  begin
    // there is just 1 experiment with 1 control
    cVec     := NewIVector(nbrE);
    cVec[1]  := nbrC;
    nVec     := NewIVector(nbrE);
    nVec[1]  := nbrN;
    xMtxData := NewRData(nbrE, cVec, nVec);
    for n := 0 to nbrN do
      xMtxData[nbrE][nbrC][n] := xData[n];
    // clean up
    DelIVector(cVec);
    DelIVector(nVec);
    // return the data
    GetXData := xMtxData
  end;

  function GetYData : TRData;
  var
    n          : Integer;
    nVec, rVec : TIVector;
    yMtxData   : TRData;
  begin
    // there is just 1 experiment with 1 response
    nVec     := NewIVector(nbrE);
    nVec[1]  := nbrN;
    rVec     := NewIVector(nbrE);
    rVec[1]  := nbrR;
    yMtxData := NewRData(nbrE, rVec, nVec);
    for n := 0 to nbrN do
      yMtxData[nbrE][nbrR][n] := yData[n];
    // clean up
    DelIVector(nVec);
    DelIVector(rVec);
    // return the data
    GetYData := yMtxData
  end;

  // -----------------------------------------------------------------------------

  // Model parameters to be solved for are:
  //   b  is the y-intercept of the line
  //   m  is the slope of the line
  // which describe  y = mx + b, the model sent to the genetic algorithm.

  procedure ResponseFn (experiment          : Integer;
                        var modelParameters : TRVector;
                        var controlData     : TRMatrix;
                        var responseData    : TRMatrix);
  var
    b, m : Real;
    n    : Integer;
  begin
    n := experiment; // gets rid of warning message, experiment is not used
    // there is just 1 experiment with 1 control and 1 response
    // assign model pararmeters to global variables
    m := modelParameters[1];    // slope
    b := modelParameters[2];    // y intercept
    // create the response data predicted by the model
    responseData := NewRMatrix(nbrR, nbrN);
    for n := 1 to nbrN do
      responseData[nbrR][n] := m * controlData[nbrC][n] + b
  end;

  // -----------------------------------------------------------------------------

  // Functions and procedures used to initialize the genetic algorithm.

  function ParameterNames : TSVector;
  var
    names : TSVector;
  begin
    names := NewSVector(nbrP);
    names[1] := 'slope m';
    names[2] := 'intercept b';
    // return the parameter names
    ParameterNames := names
  end;

  function VaryParameters : TBVector;
  var
    data : TBVector;
  begin
    SetLength(data, Succ(nbrP));
    data[1] := True;            // vary m
    data[2] := True;            // vary b
    // return the parameter properties
    VaryParameters := data
  end;

  function MyFixedParameters (var variedParams : TRVector) : TRVector;
  var
    data : TRVector;
  begin
    SetLength(data, 0);            // all parameters will vary, none are fixed
    // return the parameters
    MyFixedParameters := data
  end;

  function MaximumParameters : TRVector;
  var
    data : TRVector;
  begin
    // the greatest (or most positive) parameter values to be considered
    data := NewRVector(nbrP);
    if mm > 0.0 then            // greatest value for m
      data[1] := 5.0 * mm
    else if mm = 0.0 then
      data[1] := 1.0
    else
      data[1] := mm / 5.0;
    if bb > 0.0 then            // greatest value for b
      data[2] := 5.0 * bb
    else if bb = 0.0 then
      data[2] := 1.0
    else
      data[2] := bb / 5.0;
    // return the parameters
    MaximumParameters := data
  end;

  function MinimumParameters : TRVector;
  var
    data : TRVector;
  begin
    // the smallest (or most negative) parameter values to be considered
    data := NewRVector(nbrP);
    if mm > 0.0 then            // smallest value for m
      data[1] := mm / 5.0
    else if mm = 0.0 then
      data[1] := -1.0
    else
      data[1] := 5.0 * mm;
    if bb > 0.0 then            // smallest value for b
      data[2] := bb / 5.0
    else if bb = 0.0 then
      data[2] := -1.0
    else
      data[2] := 5.0 * bb;
    // return the parameters
    MinimumParameters := data
  end;

  function AlienParameters : TRVector;
  var
    data : TRVector;
  begin
    // this is a guess at what one expects the parameters to be
    data := NewRVector(nbrP);
    data[1] := 0.75 * mm;              // guess at m
    data[2] := 1.25 * bb;              // guess at b
    // return the parameters
    AlienParameters := data
  end;

  // ---------------------------------------------------------------------------

  procedure Run;
  var
    alienPa     : TRVector;
    c           : Char;
    colony      : TColony;
    converged   : Boolean;
    elite       : TCreature;
    expCtrl     : TRData;
    expResp     : TRData;
    finished    : Boolean;
    fitness     : Real;
    fixedPa     : TVectorFn;
    generation  : Integer;
    kurt        : Real;
    maxPara     : TRVector;
    mean        : Real;
    minPara     : TRVector;

    model       : TModel;
    myResp      : TRMatrix;

    param       : TRVector;
    paramNames  : TSVector;
    sigFigs     : Integer;
    skew        : Real;
    stdDev      : Real;
    varyPar     : TBVector;
  begin
    converged  := False;
    mean       := 0.0;
    stdDev     := 0.0;
    skew       := 0.0;
    kurt       := 0.0;
    fitness    := 0.0;
    param      := NewRVector(nbrP);
    sigFigs    := 5;
    paramNames := ParameterNames;
    ReadInData;
    expCtrl := GetXData;
    expResp := GetYData;
    alienPa := AlienParameters;
    fixedPa := @MyFixedParameters;
    maxPara := MaximumParameters;
    minPara := MinimumParameters;
    varyPar := VaryParameters;

    // I get an access violation if I do not call model first - don't know why
    model  := @ResponseFn;
    myResp := expResp[1];
    model(1, alienPa, expCtrl[1], myResp);

    // start optimization
    colony := TColony.Create(model, expCtrl, expResp, varyPar, fixedPa,
                             alienPa, minPara, maxPara, sigFigs, paramNames);
    // the first generation
    elite   := colony.BestCreature;
    fitness := elite.GetFitness;
    elite.GetAllParameters(param);
    WriteLn;
    WriteLn('The following parameters are being fit:');
    WriteLn('   m    is the slope of a straight line.');
    WriteLn('   b    is the y intercept for the line.');
    WriteLn;
    WriteLn('For the first generation:');
    WriteLn('   The elite creature had a fitness of ', fitness:8:5, '.');
    WriteLn('   Its parameters were:');
    WriteLn('     m = ', param[1]:6:4, ' in [',
            minPara[1]:6:4, ', ', maxPara[1]:6:4, ']');
    WriteLn('     b = ', param[2]:6:4, ' in [',
            minPara[2]:6:4, ', ', maxPara[2]:6:4, ']');
    colony.FitnessMoments(mean, stdDev, skew, kurt);
    WriteLn('   The population as a whole had fitness statistics of: ');
    WriteLn('      mean = ', mean:8:5, '  std devivation = ', stdDev:8:5);
    WriteLn('      skew = ', skew:8:5, '  and kurtosis   = ', kurt:8:5);
    finished   := False;
    generation := 1;
    // the following generations
    while not finished do begin
      Inc(generation);
      WriteLn;
      colony.NextGeneration(converged);
      elite   := colony.BestCreature;
      fitness := elite.GetFitness;
      elite.GetAllParameters(param);
      WriteLn('For generation ', generation, ':');
      WriteLn('   The elite creature had a fitness of ', fitness:8:5, '.');
      WriteLn('   Its parameters were:');
      WriteLn('     m = ', param[1]:6:4, ' in [',
              minPara[1]:6:4, ', ', maxPara[1]:6:4, ']');
      WriteLn('     b = ', param[2]:6:4, ' in [',
              minPara[2]:6:4, ', ', maxPara[2]:6:4, ']');
      colony.FitnessMoments(mean, stdDev, skew, kurt);
      WriteLn('   The population as a whole had fitness statistics of: ');
      WriteLn('      mean = ', mean:8:5, '  std devivation = ', stdDev:8:5);
      WriteLn('      skew = ', skew:8:5, '  and kurtosis   = ', kurt:8:5);
      WriteLn;
      Write('Generation ', generation, ' of the genetic algorithm has ');
      if converged then begin
        writeLn('converged.');
        finished := True
      end
      else begin
        writeLn('not converged.');
        Write('Do you want to construct the next generation? [Y,N]  ');
        Readln(c);
        if (c = 'y') or (c = 'Y') then
          finished := False
        else 
          finished := True
      end
    end
  end;

begin
  a0    := 0.0;
  a1    := 0.0;
  bb    := 0.0;
  mm    := 0.0;
  xData := NewRVector(nbrN);
  yData := NewRVector(nbrN);
  Run
end.

