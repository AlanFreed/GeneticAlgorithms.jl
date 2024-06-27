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
   Port began on - April     1, 2014
   Last modified - November 21, 2015
--------------------------------------------------------------------------------
   Version 1.2 - 1.3
   Adopted 'Arrays' for creating and managing vector arrays.  
   This library calls librarys 'Arrays' and 'rng' written by the author.
--------------------------------------------------------------------------------
   This file provides the model interface for the genetic algorithm.
***************************************************************************** *)

// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

unit gaSpecies;

interface

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    Classes,       // allows creation of objects
    SysUtils,      // the system's exception handler
    gaCore;        // genetic algorithm's basic type definitions

  //  argument arrays are passed by reference for C compliance

  type
    TModel = procedure(
      experiment          : Integer;      // IN     exp is in range 1..dimE
      var modelParameters : TRVector;     // IN     [1..dimP]
      var controlData     : TRMatrix;     // IN     [1..dimNC[e]][1..dimNS[e]]
      var responseData    : TRMatrix);    // OUT    [1..dimNR[e]][1..dimNS[e]]
    { experiment      : used to determine which experiment is to be modeled
      modelParameters : all model parameters (both fixed and adjustable)
      controlData     : a sequence of experimentally controlled values
                        [controlVariable][valueAtASamplingInstant]
      responseData    : predicted responses for these controls and parameters
                        [responseVariable][valueAtASamplingInstant] }

implementation

  // This is a virtual interface.  There is no implementation.

end.
