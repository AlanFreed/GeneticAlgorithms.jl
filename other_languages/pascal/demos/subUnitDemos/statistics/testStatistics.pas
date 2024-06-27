// global directives
{$MODE OBJFPC}    // allows objects, classes, interfaces and exception handling
{$H+}             // use ANSI strings
{$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
{$PIC ON}         // create Position Independent Code - for Unix libraries

program testStatistics;

  uses
    gaCore, gaProbabilities, gaStatistics;

  const
    SB = 1;      // denotes a Johnson SB (bounded)   distribution
    SL = 2;      // denotes a Johnson SL (lognormal) distribution
    SN = 3;      // denotes a Johnson SN (normal)    distribution
    SU = 4;      // denotes a Johnson SU (unbounded) distribution
    ST = 5;      // denotes a Johnson SB (bounded)   distribution
                 //   near the boundary of skewness/kurtosis inadmissibility

  type
    TRVector = array of Real;
    TRMatrix = array of TRVector;

  // imported procedures, etc.

  function NewRVector (len : Integer) : TRVector; external 'Arrays';
  // Creates a real vector indexing from 1 to len whose elements are all 0.0.

  procedure DelRVector (var v : TRVector); external 'Arrays';
  // Deletes the information held by dynamic real vector v.

  function NewRMatrix (rows, columns : Integer) : TRMatrix; external 'Arrays';
  // Creates a real matrix indexing as [1..rows][1..columns]
  // whose elements are assigned values of 0.0.

  procedure DelRMatrix (var m : TRMatrix); external 'Arrays';
  // Deletes the information held by real matrix m.

  procedure JohnsonStatistics (var v            : TRVector;
                               out distribution : Integer;
                               out gamma        : Real;
                               out delta        : Real;
                               out xi           : Real;
                               out lambda       : Real); external 'rng';
  { given a vector or reals, returns the Johnson distribution
      and its statistics, where the distributions are:
        SB  =>  Z = gamma + delta*ln[(X - xi)/(xi + lambda - X)],
        SL  =>  Z = gamma + delta*ln(X - xi),
        SN  =>  Z = gamma + delta*(X - xi)/lambda
        SU  =>  Z = gamma + delta*arcsinh[(X - xi)/lambda]
      including
        ST  =>  Z = delta*ln[(X - xi)/(xi + lambda - X)],
          is SB distribution in transition zone boardering inadmissible domain.
          Instances of this category may not do good job of representing data.
    with statistics
        gamma   is a shape parameter, akin to:   1.0 / standard deviation
        delta   is a shape parameter, akin to: -mean / standard deviation
        xi      is a location parameter
        lambda  is a scale parameter  }

  // Function exported by the random library.

  function RandomNormal (mean, standardDeviation : Real) : Real; external 'rng';

   procedure Test;
   var
      i      : Integer;
      delta  : Real;
      dist   : Integer;
      disNam : String;
      exp    : TRVector;
      gamma  : Real;
      kurt   : Real;
      lambda : Real;
      model  : TRVector;
      mean   : Real;
      mu     : Real;
      random : TRVector;
      sigma  : Real;
      skew   : Real;
      stdDev : Real;
      r      : Real;
      xi     : Real;
   begin
      dist   := SB;
      exp    := nil;
      model  := nil;
      random := nil;
      kurt   := 0.0;
      mu     := 0.0;
      sigma  := 0.0;
      skew   := 0.0;
      gamma  := 0.0;
      delta  := 0.0;
      xi     := 0.0;
      lambda := 0.0;
      WriteLn;
      WriteLn('The program tests the statistics unit.');
      WriteLn;
      random := NewRVector(5);
      for i := 1 to 5 do
         random[i] := i;
      SampleStatistics(random, mu, sigma);
      WriteLn('The data set {1, 2, 3, 4, 5} has a');
      WriteLn('   mean value of  3 which was computed to be ', mu);
      WriteLn('   variance of  2.5 which was computed to be ', sigma*sigma);
      HigherStatistics(random, skew, kurt);
      WriteLn('   skewness of  0.0 which was computed to be ', skew);
      WriteLn('   kurtosis of -1.3 which was computed to be ', kurt);
      WriteLn;
      Write('Enter a mean value for a population within [0,1]:     ');
      mu   := 0.0;
      mean := 0.0;
      ReadLn(mean);
      Write('Enter the associated standard deviation less than 1:  ');
      sigma  := 0.0;
      stdDev := 0.0;
      ReadLn(stdDev);
      exp := NewRVector(1000);
      for i := 1 to 1000 do
        exp[i] := RandomProbability(mean, stdDev);
      WriteLn;
      WriteLn('Two randomly generated populations with these statistics have:');
      SampleStatistics(exp, mu, sigma);
      WriteLn('   population 1: a mean value of:         ', mu);
      WriteLn('                 a standard deviation of: ', sigma);
      HigherStatistics(exp, skew, kurt);
      WriteLn('                 a skewness of:           ', skew);
      WriteLn('                 a kurtosis of:           ', kurt);
      WriteLn('be patient');
      JohnsonStatistics(exp, dist, gamma, delta, xi, lambda);
      if dist = SB then
        disNam := 'SB'
      else if dist = SL then
        disNam := 'SL'
      else if dist = SN then
        disNam := 'SN'
      else if dist = SU then
        disNam := 'SU'
      else
        disNam := 'ST';
      WriteLn('                 a Johnson ', disNam, ' gamma of    ', gamma);
      WriteLn('                 a Johnson ', disNam, ' delta of    ', delta);
      WriteLn('                 a Johnson ', disNam, ' xi of       ', xi);
      WriteLn('                 a Johnson ', disNam, ' lambda of   ', lambda);
      model := NewRVector(1000);
      for i := 1 to 1000 do
        model[i] := RandomNormal(mean, stdDev);
      WriteLn;
      SampleStatistics(model, mu, sigma);
      WriteLn('   population 2: a mean value of:         ', mu);
      WriteLn('                 a standard deviation of: ', sigma);
      HigherStatistics(model, skew, kurt);
      WriteLn('                 a skewness of:           ', skew);
      WriteLn('                 a kurtosis of:           ', kurt);
      WriteLn('be patient');
      JohnsonStatistics(model, dist, gamma, delta, xi, lambda);
      if dist = SB then
        disNam := 'SB'
      else if dist = SL then
        disNam := 'SL'
      else if dist = SN then
        disNam := 'SN'
      else if dist = SU then
        disNam := 'SU'
      else
        disNam := 'ST';
      WriteLn('                 a Johnson ', disNam, ' gamma of    ', gamma);
      WriteLn('                 a Johnson ', disNam, ' delta of    ', delta);
      WriteLn('                 a Johnson ', disNam, ' xi of       ', xi);
      WriteLn('                 a Johnson ', disNam, ' lambda of   ', lambda);
      WriteLn;
      r := RMSE(exp, model);
      WriteLn('The RMSE between these two distributions is: ', r);
      WriteLn('Population 1 is derived from random Johnson SB variates.');
      WriteLn('Pouluation 2 is derived from random normal variates.');
      WriteLn
   end;

  procedure NoiseAboutLine;
  var
    a0, a1, bb, detM,
    mm, rM, rT, sd,
    summation, xMax,
    xMin             : Real;
    n, nbrN,
    nbrP, p, q       : Integer;
    rhs, xExp,
    yExp, yMod, yThy : TRVector;
    mInv, mMtx, zMtx : TRMatrix;
  begin
    bb := 0.0;
    mm := 0.0;
    sd := 0.0;
    nbrN := 1000;
    nbrP := 2;
    xMin := 0.0;
    xMax := 10.0;
    WriteLn;
    Write('What is the y-intercept of the line to be? ');
    ReadLn(bb);
    Write('What is the slope of the line to be?       ');
    ReadLn(mm);
    Write('What is the standard deviation for noise?  ');
    ReadLn(sd);
    xExp := NewRVector(nbrN);
    yExp := NewRVector(nbrN);
    yThy := NewRVector(nbrN);
    for n := 1 to nbrN do begin
      xExp[n] := xMin + (n - 1) * (xMax - xMin) / (nbrN - 1);
      yThy[n] := mm * xExp[n] + bb;
      yExp[n] := RandomNormal(yThy[n], sd)
    end;
    // perform a linear regression of these data
    zMtx := NewRMatrix(nbrN, nbrP);
    for n := 1 to nbrN do begin
      zMtx[n][1] := 1.0;
      zMtx[n][2] := xExp[n]
    end;
    // M = Z^T.Z
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
    // rhs = Z^T.y
    rhs := NewRVector(nbrP);
    for p := 1 to nbrP do begin
      summation := 0.0;
      for n := 1 to nbrN do
        summation := summation + zMtx[n][p] * yExp[n];
      rhs[p] := summation
    end;
    // to solve  Z^T.Z.a = Z^T.y  for  a compute  (Z^T.Z)^-1
    mInv := NewRMatrix(nbrP, nbrP);
    detM := mMtx[1][1] * mMtx[2][2] - mMtx[1][2] * mMtx[2][1];
    mInv[1][1] :=  mMtx[2][2] / detM;
    mInv[1][2] := -mMtx[1][2] / detM;
    mInv[2][1] :=  mInv[1][2];
    mInv[2][2] :=  mMtx[1][1] / detM;
    // solve for the model coefficients, i.e.,  y = a0 + a1 x
    a0 := mInv[1][1] * rhs[1] + mInv[1][2] * rhs[2];
    a1 := mInv[2][1] * rhs[1] + mInv[2][2] * rhs[2];
    // determine the y values for this linear regression
    yMod := NewRVector(nbrN);
    for n := 1 to nbrN do
       yMod[n] := a0 + a1 * xExp[n];
    WriteLn('The equation  asigned  was:  y = ', mm:8:6, '*x + ', bb:8:6);
    WriteLn('The regressed equation was:  y = ', a1:8:6, '*x + ', a0:8:6);
    WriteLn('Root mean square error (RMSE) between data a straight lines:');
    rT := RMSE(yExp, yThy);
    WriteLn('   random noise about the  asigned  straight line, RMSE = ', rT:8:6);
    rM := RMSE(yExp, yMod);
    WriteLn('   random noise about the regressed straight line, RMSE = ', rM:8:6);
    // clean up
    DelRMatrix(mInv);
    DelRMatrix(mMtx);
    DelRVector(rhs);
    DelRVector(xExp);
    DelRVector(yExp);
    DelRVector(yMod);
    DelRVector(yThy);
    DelRMatrix(zMtx)
  end;

begin
  Test;
  NoiseAboutLine
end.

