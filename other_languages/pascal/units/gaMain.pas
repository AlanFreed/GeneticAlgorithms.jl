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
   Port began on - June     24, 2014
   Last modified - November 21, 2015
--------------------------------------------------------------------------------
   Version 1.2 - 1.3
   Adopted 'Arrays' for creating and managing vector arrays.  
   This library calls librarys 'Arrays' and 'rng' written by the author.
   Version 1.3.3
   'fixedParameters' made a fn call - allow them to depend on varied parameters.
--------------------------------------------------------------------------------
   The main user interface to GenAlg.
***************************************************************************** *)

// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

unit gaMain;

interface

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    gaCore,        // genetic algorithm's basic type definitions
    gaSpecies;     // genetic algorithm's interface for the model

  // argument arrays are passed by reference in order to be C compliant

  // create the optimizer, establish first generation, all arguments are inputs
  procedure StartGeneticAlgorithm (numericalModel        : TModel;
                                   var expControlData    : TRData;
                                   var expResponseData   : TRData;
                                   var varyParameters    : TBVector;
                                   fixedParameters       : TVectorFn;
                                   var alienParameters   : TRVector;
                                   var minParameters     : TRVector;
                                   var maxParameters     : TRVector;
                                   nbrSignificantFigures : Integer;
                                   var parameterNames    : TSVector;
                                   reportFileName        : String);

  // advance the colony to next generation and test for algorithmic convergence
  procedure AdvanceToNextGeneration (out converged : Boolean);

  // discontinue the analysis
  procedure StopGeneticAlgorithm;
  
  // remove the genetic algorithm's data structure from memory
  procedure DeleteGeneticAlgorithm;

  // retreive cumulative sums for algorithmic data over all of the generations
  procedure GeneticAlgorithmData (out nbrCreatures   : Integer;
                                  out nbrGenerations : Integer;
                                  out nbrImmigrants  : Integer;
                                  out nbrCrossovers  : Integer;
                                  out nbrMutations   : Integer);

  // retreive statistics for the elite creature in the current generation
  procedure EliteCreatureData (out fitness    : Real;
                               var parameters : TRVector;
                               var genome     : TSVector);
                               

  // retreive curves for the fitted model and the data used for parameterization
  procedure EliteCreatureFits (selectExperiment  : Integer;           // 1..dimE
                               var expControl    : TRMatrix;          // output
                               var expResponse   : TRMatrix;          // output
                               var modelControl  : TRMatrix;          // input
                               var modelResponse : TRMatrix);         // output

  // moment statistics for population fitness over current generation
  procedure FitnessMoments (out mean     : Real;
                            out stdDev   : Real;
                            out skewness : Real;
                            out kurtosis : Real);

  // Johnson statistics for population fitness over current generation
  procedure FitnessStatistics (out distribution : Integer;
                               out gamma        : Real;
                               out delta        : Real;
                               out xi           : Real;
                               out lambda       : Real);
                               
implementation

  uses
    Classes,       // allows creation of objects
    SysUtils,      // the system's exception handler
    gaColonies,    // genetic algorithm's manager of colonies
    gaCreatures,   // genetic algorithm's manager of creatures
    gaStatistics,  // genetic algorithm's statistics package
    gaWriteReport; // genetic algorithm's report writer

  var
    colony       : TColony;
    expCtrlData  : TRData;
    expRespData  : TRData;
    fileName     : String;
    model        : TModel;
    
  // procedures local to this unit

  procedure StartGeneticAlgorithm (numericalModel        : TModel;
                                   var expControlData    : TRData;
                                   var expResponseData   : TRData;
                                   var varyParameters    : TBVector;
                                   fixedParameters       : TVectorFn;
                                   var alienParameters   : TRVector;
                                   var minParameters     : TRVector;
                                   var maxParameters     : TRVector;
                                   nbrSignificantFigures : Integer;
                                   var parameterNames    : TSVector;
                                   reportFileName        : String);
  begin
    colony := TColony.Create(numericalModel,
                             expControlData,
                             expResponsedata,
                             varyParameters,
                             fixedParameters,
                             alienParameters,
                             minParameters,
                             maxParameters,
                             nbrSignificantFigures,
                             parameterNames);
    model       := numericalModel;
    expCtrlData := expControlData;
    expRespData := expResponseData;
    fileName    := reportFileName + '.txt';
    ReportHeader(colony, fileName);
    ReportBody  (colony, fileName)
  end;

  procedure AdvanceToNextGeneration (out converged : Boolean);
  begin
    colony.NextGeneration(converged);
    ReportBody(colony, fileName);
    Write('.')
  end;

  procedure StopGeneticAlgorithm;
  begin
    ReportFooter(colony, fileName)
  end;
  
  procedure DeleteGeneticAlgorithm;
  begin
    DelRData(expCtrlData);
    DelRData(expRespData);
    fileName := '';
    colony.DeleteCreatures;
    colony := nil;
    model  := nil
  end;

  procedure GeneticAlgorithmData (out nbrCreatures   : Integer;
                                  out nbrGenerations : Integer;
                                  out nbrImmigrants  : Integer;
                                  out nbrCrossovers  : Integer;
                                  out nbrMutations   : Integer);
  begin
    colony.PopulationData(nbrCreatures, nbrGenerations,
                          nbrImmigrants, nbrCrossovers, nbrMutations)
  end;

  procedure EliteCreatureData (out fitness    : Real;
                               var parameters : TRVector;
                               var genome     : TSVector);
  var
    elite : TCreature;
  begin
    elite   := colony.BestCreature;
    fitness := elite.GetFitness;
    SetLength(parameters, 0);
    elite.GetAllParameters(parameters);
    elite.GetGenome(genome)
  end;

  procedure EliteCreatureFits (selectExperiment  : Integer;
                               var expControl    : TRMatrix;    // output
                               var expResponse   : TRMatrix;    // output
                               var modelControl  : TRMatrix;    // input
                               var modelResponse : TRMatrix);   // output
  var
    c, r, s : Integer;
    elite   : TCreature;
    param   : TRVector;
  begin
    if modelControl = nil then raise
      Exception.Create('argument "modelControl" had not been created');
    if (selectExperiment < 1) or (selectExperiment > colony.dimE) then raise
      Exception.Create('the selected experiment does not exist');
    // extract the pertenent experimental data for returning
    SetLength(expControl, 0);
    expControl := NewRMatrix(colony.dimNC[selectExperiment],
                             colony.dimNS[selectExperiment]);
    for c := 1 to colony.dimNC[selectExperiment] do
      for s := 1 to colony.dimNS[selectExperiment] do
        expControl[c][s] := expCtrlData[selectExperiment][c][s];
    SetLength(expResponse, 0);
    expResponse := NewRMatrix(colony.dimNR[selectExperiment],
                              colony.dimNS[selectExperiment]);
    for r := 1 to colony.dimNR[selectExperiment] do
      for s := 1 to colony.dimNS[selectExperiment] do
        expResponse[r][s] := expRespData[selectExperiment][r][s];
    // prepare to run the model
    elite := colony.BestCreature;
    SetLength(param, 0);
    elite.GetAllParameters(param);
    SetLength(modelResponse, 0);
    model(selectExperiment, param, modelControl, modelResponse);
    // clean up
    DelRVector(param)
  end;

  procedure FitnessMoments (out mean     : Real;
                            out stdDev   : Real;
                            out skewness : Real;
                            out kurtosis : Real);
  begin
    colony.FitnessMoments(mean, stdDev, skewness, kurtosis)
  end;

  procedure FitnessStatistics (out distribution : Integer;
                               out gamma        : Real;
                               out delta        : Real;
                               out xi           : Real;
                               out lambda       : Real);
  begin
    colony.FitnessStatistics(distribution, gamma, delta, xi, lambda)
  end;
  
begin
  colony   := nil;
  fileName := '';
  model    := nil;
  SetLength(expCtrlData, 0);
  SetLength(expRespData, 0)
end.

