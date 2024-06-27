(* ****************************************************************************
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
   Port began on - March    19, 2014
   Last modified - November 21, 2015
--------------------------------------------------------------------------------
   Version 1.2 - 1.3
   Adopted 'Arrays' for creating and managing vector arrays.  
   This library calls librarys 'Arrays' and 'rng' written by the author.
--------------------------------------------------------------------------------
   This unit establishes the basic statistics used by the genetic algorithm.
***************************************************************************** *)

// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

unit gaStatistics;

interface

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    gaCore;  // genetic algorithm's basic type definitions

  // argument arrays are passed by reference
  // values held in element v[0] are NOT entered into any of these calculations

  procedure SampleStatistics (var v : TRVector; out mu, sigma : Real);
  { Supplies the sample mean 'mu' and standard deviation 'sigma' given by
      mu = E(v)  and  sigma = sqrt[(N/(N-1))[E(v^2) - E^2(v)]] }

  procedure HigherStatistics (var v : TRVector; out skewness, kurtosis : Real);
  { Supplies the: skewness = E[(X-mu)^3]/(E[(X-mu)^2])^(3/2)
                  kurtosis = E[(X-mu)^4]/(E[(X-mu)^2])^2 - 3
     where        skewness < 0 implies heavily tailed to the left
                  skewness = 0 implies evenly tailed like normal distribution
                  skewness > 0 implies heavily tailed to the right
     and          kurtosis < 0 implies peaked less than a normal distribution
                  kurtosis = 0 implies peaked just like a normal distribution
                  kurtosis > 0 implies peaked more than a normal distribution
    The kurtosis returned is the so-called 'excess kurtosis'. }

  function RMSE (var experiment, model : TRVector) : Real;
  { The Root Mean Square Error (RMSE) is defined as
                  /  1   N  /  mod     exp \ 2 \
       RMSE = SQRT| --- SUM | Y    -  Y    |   |
                  \ N-2 n=1 \  n       n   /   /
    it is the distance between model and experiment in the units supplied.
    It is scaled by N-2 because the mean is also determined from these data. }

implementation

  uses
    Classes,    // allows creation of objects
    SysUtils;   // exception handler
    
  var
    emptyRVector : TRVector;
    
  procedure NormalStatistics (var v     : TRVector;
                              out mu    : Real;
                              out sigma : Real); external 'rng';
  { given a vector or reals, returns the Gauss statistics
      mu    = E(v)                           sample mean
      sigma = sqrt[(N/(N-1)) E[(X - mu)^2]   sample standard deviation }

  procedure MomentStatistics (var v      : TRVector;
                              out gamma1 : Real;
                              out gamma2 : Real); external 'rng';
  { given a vector of reals, returns the first two moment statistics
      gamma1 = E[(X - mu)^3] / (E[(X - mu)^2])^1.5    skewness
      gamma2 = E[(X - mu)^4] / (E[(X - mu)^2])^2 - 3  excess kurtosis }

  // exported procedures and function of this unit

  procedure SampleStatistics (var v : TRVector; out mu, sigma : Real);
  begin
    if (v = nil) or (v = emptyRVector) then raise Exception.Create
      ('supplied vector sent to SampleStatistics has no length');
    NormalStatistics(v, mu, sigma)
  end;

  procedure HigherStatistics (var v : TRVector; out skewness, kurtosis : Real);
  begin
    if (v = nil) or (v = emptyRVector) then raise Exception.Create
      ('supplied vector sent to HigherStatistics has no length');
    MomentStatistics(v, skewness, kurtosis)
  end;

  function RMSE (var experiment, model : TRVector) : Real;
    var
      dim, n   : Integer;
      dif, mse : Real;
  begin
    if ((experiment = emptyRVector) or (experiment = nil)) or 
       ((model = emptyRVector) or (model = nil)) then raise
      Exception.Create('a supplied array to RMSE has no length');
    if Length(experiment) <> Length(model) then raise
      Exception.Create('arrays must have equal length to compute RMSE');
    if Length(model) <= 3 then raise
      Exception.Create('the supplied arrays to RMSE are not long enough');
    dim := LenRVector(model);
    mse := 0.0;
    for n := 1 to dim do begin  
      dif := model[n] - experiment[n];
      mse := mse + dif * dif
    end;
    if mse > 0.0 then
      RMSE := Sqrt(mse / (dim - 2))
    else
      RMSE := machEp
  end;
  
begin

  SetLength(emptyRVector, 0)

end.
