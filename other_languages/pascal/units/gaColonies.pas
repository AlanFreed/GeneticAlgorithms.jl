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
   Port began on - April     3, 2014
   Last modified - November 21, 2015
--------------------------------------------------------------------------------
   Version 1.2 - 1.3
   Adopted 'Arrays' for creating and managing vector arrays.  
   This library calls librarys 'Arrays' and 'rng' written by the author.
   Version 1.3.3
   'fixedParameters' made a fn call - allow them to depend on varied parameters.
--------------------------------------------------------------------------------
   This belongs to the colony part of this genetic algorithm, which includes:
      gaCreatures   - this is the interface between genetics and optimization
      gaColonies    - establishes a collection of creatures
      gaMain        - the basic driver of the gentic algorithm
      gaParameters  - calculates statistics for the estimated parameters
      gaWriteReport - writes out a report for the user post analysis
***************************************************************************** *)

// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

unit gaColonies;

interface

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    gaCore,       // genetic algorithm's basic type definitions
    gaCreatures,  // genetic algorithm's creature class
    gaSpecies;    // genetic algorithm's virtual model interface

  type
    TColony = class
    private
      expCtrlData : TRData;       // exp data used for controlling a model
      expRespData : TRData;       // exp data used to compare against a model
      varyParam   : TBVector;     // specifies which parameters can be varied
      fixedParam  : TVectorFn;    // values for those parameters held fixed
      minFitParam : TRVector;     // lower bounds on the adjusted parameters
      maxFitParam : TRVector;     // upper bounds on the adjusted parameters
      sigFigures  : Integer;      // number of significant figures being sought

      probImgrant : Real;         // probability an immigrant enters population

      convergeAt  : Integer;      // number of generations to convergence
      crossovers  : Integer;      // tracks the number of crossover events
      generations : Integer;      // tracks the number of generations
      immigrants  : Integer;      // tracks the number of immigrants
      mutations   : Integer;      // tracks the number of mutation events
      population  : Integer;      // number of creatures in any given generation

      fitnessVec  : TRVector;     // fitness of all living creatures

      elite       : TCreature;    // most fit creature in the population
      adults      : array of TCreature;   // current population of creatures
      children    : array of TCreature;   // next population of creatures
      contestants : array of TCreature;   // tournament participants

      model       : TModel;       // model whose parameters are being sought

      procedure Evaluate (var c : TCreature);
      procedure TournamentPlay (var mostFit : TCreature);
      procedure Mate;

    public
      // dimensioning of the data structures
      dimE  : Integer;  // number of experiments
      dimP  : Integer;  // number of parameters, varied and fixed
      dimPV : Integer;  // number of parameters that are allowed to vary
      dimR  : Integer;  // number of response variables over all experiments
      dimS  : Integer;  // number of fitted points over all experiments

      // all dimXX vectors index on the experiment number which starts at 1
      dimNC : TIVector;     { # of control variables per experiment
                                    index: experiment number
                                    value: # of control variables: 1,2,... }
      dimNR : TIVector;     { # of response variables per experiment
                                    index: experiment number
                                    value: # of response variables: 1,2,... }
      dimNS : TIVector;     { # of sampling points taken per experiment
                                    index: experiment number
                                    value: # of sampling points plus 1 for IC
                                    initial condition is stored in element 0}

      names : TSVector; // string descriptions of the varied parameters

      function BestCreature : TCreature;
      { retreive the most fit member of the current generation }

      procedure PopulationData (out nbrCreatures   : Integer;
                                out nbrGenerations : Integer;
                                out nbrImmigrants  : Integer;
                                out nbrCrossovers  : Integer;
                                out nbrMutations   : Integer);
      { retreive run data from the optimizer }

      procedure FitnessMoments (out mean     : Real;
                                out stdDev   : Real;
                                out skewness : Real;
                                out kurtosis : Real);
      { moment statistics for fitness of the overall population }

      procedure FitnessStatistics (out distribution : Integer;
                                   out gamma        : Real;
                                   out delta        : Real;
                                   out xi           : Real;
                                   out lambda       : Real);
      { Johnson statistics for fitness of the overall population  }
      
      
      
      procedure GetBoundaries (var lower, upper : TRVector);
      { retrieve the lower and upper boundaries for parameters search domain }
      
      function VaryParameters : TBVector;
      { specifies if each parameter is fixed or allowed to vary }
        
      procedure DeleteCreatures;
      { removes the dynamically allocated data structures within colony }

      //  argument arrays are passed by reference
      constructor Create (numericalModel        : TModel;
                          var expControlData    : TRData;
                          var expResponseData   : TRData;
                          var varyThisParameter : TBVector;
                          fixedParameters       : TVectorFn;
                          var alienParameters   : TRVector;
                          var minFitParameters  : TRVector;
                          var maxFitParameters  : TRVector;
                          significantFigs       : Integer;
                          var parameterNames    : TSVector);
      { all parameter vectors are to index from 1 }

      procedure NextGeneration (out converged : Boolean);
      { advance the colony to its next generation and test for convergence }
    end;

implementation

  uses
    Classes,         // for constructing objects
    SysUtils,        // exception handler
    Math,            // basic math routines
    gaGenes,         // genetic algorithm's gene class
    gaProbabilities, // genetic algorithm's probability routines
    gaStatistics;    // genetic algorithm's statistical routines
  
  procedure Quartiles (var v       : TRVector;
                       out minimum : Real;
                       out first   : Real;
                       out second  : Real;
                       out third   : Real;
                       out maximum : Real); external 'rng';
  { given a vector of reals, returns the three quartiles of the data set
      first  = median of the lower half of the data set
      second = median of the data set
      third  = median of the upper half of the data set }

  procedure JohnsonStatistics (var v            : TRVector;
                               out distribution : Integer;
                               out gamma        : Real;
                               out delta        : Real;
                               out xi           : Real;
                               out lambda       : Real); external 'rng';
  { given a vector or reals, returns the Johnson statistics for
      distribution  is a CONST with permissible values:
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

  // ---------------------------------------------------------------------------

  { private routines}

  { the purpose of this procedure is to determine the fitness of a creature }
  procedure TColony.Evaluate (var c : TCreature);
    var
      e, r, s    : Integer;
      expCtrl    : TRMatrix;
      expResp    : TRMatrix;
      fitness    : Real;
      maxResp    : Real;
      minResp    : Real;
      modResp    : TRMatrix;
      parameters : TRVector;
      rmseMtx    : TRData;
      sVec, vVec : TIVector;
  begin
    if c = nil then raise 
      Exception.Create('the supplied creature does not exist');
    sVec := NewIVector(dimE);
    vVec := NewIVector(dimE);
    for e := 1 to dimE do begin
      sVec[e] := 1;
      vVec[e] := dimNR[e]
    end;
    rmseMtx := NewRData(dimE, vVec, sVec);
    parameters := nil;
    c.GetAllParameters(parameters);
    fitness := 0.0;
    for e := 1 to dimE do begin
      expCtrl := CopyRMatrix(expCtrlData[e]);
      expResp := CopyRMatrix(expRespData[e]);
      modResp := NewRMatrix(dimNR[e],dimNS[e]);
      model(e, parameters, expCtrl, modResp);
      for r := 1 to dimNR[e] do begin
        // get the root mean squared errors
        rmseMtx[e][r][1] := RMSE(expResp[r], modResp[r]);
        // calculate the range for the experimental response
        minResp := expResp[r][1];
        for s := 2 to dimNS[e] do
          if expResp[r][s] < minResp then
            minResp := expResp[r][s];
        maxResp := expResp[r][1];
        for s := 2 to dimNS[e] do
          if expResp[r][s] > maxResp then
            maxResp := expResp[r][s];
        // assign fitness as a harmonic mean of RMSEs normalized by their ranges
        fitness := fitness + (maxResp - minResp) / rmseMtx[e][r][1];
      end;
      // cleanup
      DelRMatrix(expCtrl);
      DelRMatrix(expResp);
      DelRMatrix(modResp)
    end;
    c.SetFitness(fitness / dimR);
    c.SetRMSE(rmseMtx);
    // clean up
    DelRVector(parameters);
    DelIVector(sVec);
    DelIVector(vVec)
  end;

  procedure TColony.TournamentPlay (var mostFit : TCreature);
    var
      creature, i : Integer;
  begin
    if IsHeads(probImgrant) then    // an immigrant migrates into the population
      begin
        mostFit := TCreature.Create(varyParam,   fixedParam,
                                    minFitParam, maxFitParam, sigFigures);
        mostFit.Procreate;
        Inc(immigrants)
      end
    else  // use tournament play to select a creature for mating
      begin
        for i := 1 to High(contestants) do begin
          creature       := RandomInteger(1, population);
          contestants[i] := adults[creature]
        end;
        mostFit := contestants[1];
        for i := 2 to High(contestants) do
          if contestants[i].GetFitness > mostFit.GetFitness then
            mostFit := contestants[i]
      end
  end;

  procedure TColony.Mate;
    var
      creature : TCreature;
      i        : Integer;
      parentA  : TCreature;
      parentB  : TCreature;
  begin
    parentA := nil;
    parentB := nil;
    // elite creature from the prior generation lives into the next generation
    for i := 1 to Pred(population) do begin
      TournamentPlay(parentA);
      TournamentPlay(parentB);
      while parentA.IsEqualTo(parentB) do begin
        parentB := nil;
        TournamentPlay(parentB)
      end;
      creature := children[i];
      creature.Conceive(parentA, parentB, mutations, crossovers);
      Evaluate(creature);
      children[i] := creature;
      parentA := nil;
      parentB := nil
    end
  end;

  { public routines }

  function TColony.BestCreature : TCreature;
  begin
    BestCreature := elite
  end;

  procedure TColony.PopulationData (out nbrCreatures   : Integer;
                                    out nbrGenerations : Integer;
                                    out nbrImmigrants  : Integer;
                                    out nbrCrossovers  : Integer;
                                    out nbrMutations   : Integer);
  begin
    nbrCreatures   := population;
    nbrGenerations := generations;
    nbrImmigrants  := immigrants;
    nbrCrossovers  := crossovers;
    nbrMutations   := mutations
  end;

  procedure TColony.FitnessMoments (out mean     : Real;
                                    out stdDev   : Real;
                                    out skewness : Real;
                                    out kurtosis : Real);
  begin
    SampleStatistics(fitnessVec, mean, stdDev);
    HigherStatistics(fitnessVec, skewness, kurtosis)
  end;

  procedure TColony.FitnessStatistics (out distribution : Integer;
                                       out gamma        : Real;
                                       out delta        : Real;
                                       out xi           : Real;
                                       out lambda       : Real);
  begin
    JohnsonStatistics(fitnessVec, distribution, gamma, delta, xi, lambda)
  end;
      
  procedure TColony.GetBoundaries (var lower, upper : TRVector);
  begin
    lower := CopyRVector(minFitParam);
    upper := CopyRVector(maxFitParam)
  end;
      
  function TColony.VaryParameters : TBVector;
  begin
    VaryParameters := CopyBVector(varyParam)
  end;
  
  procedure TColony.DeleteCreatures;
    var
      i : Integer;
  begin
    DelBVector(varyParam);
    DelRVector(minFitParam);
    DelRVector(maxFitParam);
    DelRVector(fitnessVec);
    DelRData(expCtrlData);
    DelRData(expRespData);
    for i := population downto 1 do
      adults[i].DeleteGenome;
    for i := population downto 0 do
      adults[i] := nil;
    SetLength(adults, 0);
    // do not call children[i].DeleteGenome in a loop - will crash
    for i := population downto 1 do
      children[i] := nil;
    SetLength(children, 0);
    for i := High(contestants) downto 1 do
      contestants[i].DeleteGenome;
    SetLength(contestants, 0);
    elite.DeleteGenome
  end;

  constructor TColony.Create (numericalModel        : TModel;
                              var expControlData    : TRData;
                              var expResponseData   : TRData;
                              var varyThisParameter : TBVector;
                              fixedParameters       : TVectorFn;
                              var alienParameters   : TRVector;
                              var minFitParameters  : TRVector;
                              var maxFitParameters  : TRVector;
                              significantFigs       : Integer;
                              var parameterNames    : TSVector);
    var
      chromosomes : Integer;
      combatants  : Integer;
      creature    : TCreature;
      e, genes, n : Integer;
      experiments : Integer;
      locate      : Integer;
      maxFit      : Real;
  begin
    if (expControlData = nil) or (expResponseData = nil) then raise
      Exception.Create('the supplied experiemental data had not been created');
    if (varyThisParameter = nil) or (alienParameters  = nil) or
       (minFitParameters  = nil) or (maxFitParameters = nil) then raise
      Exception.Create
        ('of the parameter data supplied, only fixedParameters can be nil');
    experiments := High(expControlData);
    if (High(expResponseData) <> experiments) then raise
      Exception.Create
        ('the number of experiments were not compatible between data sets');

    dimE  := experiments;
    dimS  := 0;
    dimNC := NewIVector(dimE);
    dimNR := NewIVector(dimE);
    dimNS := NewIVector(dimE);
    dimR  := 0;
    for e := 1 to dimE do begin
      dimNC[e] := High(expControlData[e]);
      dimNR[e] := High(expResponseData[e]);
      dimNS[e] := High(expResponseData[e][1]);
      if High(expControlData[e][1]) <> High(expResponseData[e][1]) then raise
        Exception.Create
          ('experimental control and response data are not same dimension');
      dimS := dimS + dimNS[e];
      dimR := dimR + dimNR[e]
    end;
    dimPV := 0;
    dimP  := High(varyThisParameter);
    for n := 1 to dimP do
      if varyThisParameter[n] then
        Inc(dimPV);
    if (High(alienParameters) <> dimPV) or (High(minFitParameters) <> dimPV)
      or (High(maxFitParameters) <> dimPV) then raise
      Exception.Create
        ('the parameters, those that can vary, have the wrong dimension');
    if High(parameterNames) <> dimP then raise
      Exception.Create
        ('the string array to name the parameters has the wrong dimension');

    expCtrlData := CopyRData(expControlData);
    expRespData := CopyRData(expResponseData);
    fixedParam  := fixedParameters;
    maxFitParam := CopyRVector(maxFitParameters);
    minFitParam := CopyRVector(minFitParameters);
    model       := numericalModel;
    sigFigures  := significantFigs;
    varyParam   := CopyBVector(varyThisParameter);
    names       := CopySVector(parameterNames);

    // introduce an 'alien' as the elite creature of the population
    elite := TCreature.Create(varyThisParameter, fixedParameters,
                              minFitParameters, maxFitParameters,
                              significantFigs);
    elite.Alien(alienParameters);
    Evaluate(elite);
    chromosomes := 0;
    genes       := 0;
    elite.GetGenetics(chromosomes, genes);

    // algorithm of D. Goldberg (2002) for estimating population size
    population  := Ceil(IntPower(alphabet, dimensionOfSchemata)
                 * dimensionOfSchemata * Ln(alphabet) + Ln(genes));
                 
    // create basic data vectors used by the genetic algorithm
    probImgrant := 1.0 / population;    // odds at 1 per generation
    try
      begin
        SetLength(adults, Succ(population));
        adults[0] := nil;
        adults[1] := elite;
        SetLength(children, population);
        children[0] := nil
      end
    except on E : EOutOfMemory do
      begin
        WriteLn('Memory error. Details: ' + E.ClassName + '/' + E.Message);
        Exit
      end
    end;
    fitnessVec    := NewRVector(population);
    fitnessVec[1] := elite.GetFitness;

    combatants := population div 50;
    if combatants < 3 then
      combatants := 3;
    try
      SetLength(contestants, Succ(combatants));
      contestants[0] := nil
    except on E : EOutOfMemory do
      begin
        WriteLn('Memory error. Details: ' + E.ClassName + '/' + E.Message);
        Exit
      end
    end;
    
    // algorithm  of D. Goldberg (2002) for estimating convergence
    convergeAt := Ceil(Sqrt(genes) * Ln(population) / Ln(combatants));

    // procreate the original colony
    crossovers  := 0;
    generations := 1;
    immigrants  := 0;
    mutations   := 0;
    for n := 2 to population do begin
      creature := TCreature.Create(varyThisParameter, fixedParameters,
                                   minFitParameters, maxFitParameters,
                                   significantFigs);
      creature.Procreate;
      Evaluate(creature);
      adults[n]     := creature;
      fitnessVec[n] := creature.GetFitness
    end;
    
    // find the elite procreated creature and place it in location [1];
    locate   := 1;
    creature := adults[1];
    maxFit   := creature.GetFitness;
    for n := 2 to population do begin
      creature := adults[n];
      if creature.GetFitness > maxFit then begin
        maxFit := creature.GetFitness;
        locate := n
      end
    end;
    elite          := adults[locate];
    adults[locate] := adults[1];
    adults[1]      := elite;
    for n := 1 to population do begin
      creature      := adults[n];
      fitnessVec[n] := creature.GetFitness;
    end;
      
    // create the children from which the next generation will arise
    for n := 1 to Pred(population) do begin
      creature := TCreature.Create(varyThisParameter, fixedParameters,
                                   minFitParameters, maxFitParameters,
                                   significantFigs);
      children[n] := creature
    end;
    
    // create the contestants that will be used in tournament play
    for n := 1 to combatants do begin
      creature := TCreature.Create(varyThisParameter, fixedParameters,
                                   minFitParameters, maxFitParameters,
                                   significantFigs);
      contestants[n] := creature
    end
  end;

  procedure TColony.NextGeneration (out converged : Boolean);
    var
      creature : TCreature;
      i        : Integer;
      locate   : Integer;
      maxFit   : Real;
  begin
    Mate;
    // prior generation's elite creature
    elite  := adults[1];
    maxFit := elite.GetFitness;
    locate := 1;
    // locate this generation's elite creature
    for i := 1 to Pred(population) do begin
      creature := children[i];
      if creature.GetFitness > maxFit then begin
        maxFit := creature.GetFitness;
        locate := i
      end
    end;
    if locate = 1 then 
      for i := 2 to population do 
        adults[i] := children[Pred(i)]
    else begin
      // assign the new elite creature
      adults[1] := children[locate];
      // begin filling in the population
      for i := 2 to Pred(locate) do 
        adults[i] := children[Pred(i)];
      // insert the prior elite creature into the new generation
      adults[locate] := elite;
      // handle the discontinuity caused by the insertion
      adults[Succ(locate)] := children[Pred(locate)];
      // assign the remaining creatures to the new generation
      for i := locate + 2 to population do 
        adults[i] := children[Pred(i)]
    end;
    // assign the fitness vector
    for i := 1 to population do begin
      creature      := adults[i];
      fitnessVec[i] := creature.GetFitness
    end;
    elite := adults[1];
    if generations < convergeAt then
      converged := False
    else
      converged := True;
    Inc(generations)
  end;

end.
