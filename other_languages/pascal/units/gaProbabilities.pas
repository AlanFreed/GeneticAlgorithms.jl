(* *****************************************************************************
   Author         - Alan D. Freed, Ph.D
   License        - GNU Lesser General Public License, vs. 3 or later
   Copyright      - (c) Alan D. Freed 2014-2016
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
   Port began on - March    30, 2014
   Last modified - November 21, 2015
--------------------------------------------------------------------------------
   Version 1.2 - 1.3
   Adopted 'Arrays' for creating and managing vector arrays.  
   This library calls librarys 'Arrays' and 'rng' written by the author.
--------------------------------------------------------------------------------
   This unit provides the probability routines used by the genetic algorithm.
***************************************************************************** *)

// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

unit gaProbabilities;

interface

  function IsHeads (probabilityOfHeads  : Real) : Boolean;
  { returns outcome of coin flip with specified odds
    heads -> true and tails -> false }

  function RandomInteger (loInt, hiInt : Integer) : Integer;
  { returns a random variate from the interval [loInt, hiInt]
    which must lie within the interval [-2^30 + 1, 2^30 - 1] }

  function RandomProbability (mu, sigma : Real) : Real;
  { returns a Johnson SB distributed random variate for a bounded distribution
      with lower bound xi set to 0
      upper bound lambda  set to 1
      with effective mean  gamma set to (0.5 - mu) / sigma
      and effective stdDev delta set to   1 / sigma
    parameters mu and sigma are to be interpreted in a 'normal' sense }

implementation

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    Classes,   // allows creation of objects
    SysUtils,  // exception handler
    gaCore;    // genetic algorithm's core routines

  // functions exported from the external library 'rng' that are used herein

  function RandomInterval (loInt, hiInt : Integer) : Integer; external 'rng';
  { Generates a random number on the [loInt, hiInt] interval }

  function RandomPosReal : Real; external 'rng';
  { Generates a random number on the (0.0, 1.0) line segment }

  function RandomJohnson (distribution : Integer;
                          gamma        : Real;
                          delta        : Real;
                          xi           : Real;
                          lambda       : Real) : Real; external 'rng';
  { Generates a random number for specified distribution type and statistics.
    Admissible values to pass as variable "distribution" include:
      SB   Johnson bounded distribution
      SL   Johnson lognormal distribution
      SN   Johnson normal distribution
      SU   Johnson unbounded distribution
      ST   Johnson SB distribution along boundary of inadmissible region
    Johnson distributions can describe the full range of admissible values
    for the first four moments of a data set:
      mean, standard deviation, skewness, and kurtosis. }

  // functions that are exported via this interface

  function IsHeads (probabilityOfHeads : Real) : Boolean;
  begin
    if ((probabilityOfHeads < 0.0) or (probabilityOfHeads > 1.0)) then raise
      Exception.Create('probability lies outside its range of [0,1]');
    if probabilityOfHeads >= RandomPosReal then
      IsHeads := True
    else
      IsHeads := False
  end;

  function RandomInteger (loInt, hiInt : Integer) : Integer;
  begin
    RandomInteger := RandomInterval(loInt, hiInt)
  end;

  function RandomProbability (mu, sigma : Real) : Real;
    var
      arg    : Real;
      delta  : Real;
      gamma  : Real;
      lambda : Real;
      xi     : Real;
      random : Real;
  begin
    if (mu <= 0.0) or (mu >= 1.0) then raise
      Exception.Create ('mean must lie within the unit interval (0,1)');
    if sigma <= 0.0 then raise
      Exception.Create ('sigma must be positive valued');
    arg := mu / (1.0 - mu);
    if arg > maxLn then
      arg := maxLn;
    if arg < minLn then
      arg := minLn;
    gamma  := -Ln(arg) / sigma;
    delta  := 1.0 / sigma;
    xi     := 0.0;
    lambda := 1.0;
    random := RandomJohnson(SB, gamma, delta, xi, lambda);
    RandomProbability := random
  end;

end.
