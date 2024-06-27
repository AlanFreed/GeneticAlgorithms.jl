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
   Port began on - March    19, 2014
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

unit gaCreatures;

interface

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    gaCore,    // genetic algorithm's basic type definitions
    gaGenome;  // genetic algorithm's genome or dna class

  //  argument arrays are passed by reference to be C compliant

  type
    TCreature = class
    private
      dimP               : Integer;
      dimPV              : Integer;
      fixedParameters    : TRVector;
      fixedParameterFn   : TVectorFn;
      genome             : TGenome;
      maxFitParameters   : TRVector;
      minFitParameters   : TRVector;
      quality            : Real;
      rmseValues         : TRData;
      significantFigures : Integer;
      varyThisParameter  : TBVector;
    public
      function IsEqualTo (c : TCreature) : Boolean;
      { compares two creatures to see if their genes are equivalent }
      
      function Copy : TCreature;
      { retrieves a deep copy of the creature (except for its probabilities) }

      procedure GetAllParameters (var parameters : TRVector);
      { retrieves all parameters that describe the creature, they index from 1 }
      
      procedure GetGenome (var genomes : TSVector);
      { retrieve the Gray encodings belonging to chromosomes of the creature }

      procedure GetGenetics (out chromosomes, genes : Integer);
      { retrieves number of genes and chromosomes in the creature }

      procedure GetProbabilities (out mutation, crossover : Real);
      { retrieves probabilities applied to chromosomes & genes in creature }

      function GetFitness : Real;
      { retrieves a fitness or quality statistic for the creature }

      procedure SetFitness (fitness : Real);
      { assigns a fitness or quality statistic for the creature }

      function GetRMSE (experiment       : Integer;
                        responseVariable : Integer) : Real;
      { retrieve RMSE for a specified response variable of an experiment }

      procedure SetRMSE (var rmseMtx : TRData);
      { assigns the RMSE for all response variables and all experiments }
      
      procedure DeleteGenome;
      { remove the dynamically allocated information held by the creature }
      
      procedure AssignFixedParameters;
      { after the genome becomes known, the fixed parameters can be assigned }

      procedure Procreate;
      { creates a creature whose genetic information is randomly assigned }

      procedure Alien (var fitParameters : TRVector);
      { creates an 'Adam' creature for the specified genotype/parameters }

      procedure Conceive (parentA, parentB       : TCreature;
                          var numberOfMutations  : Integer;
                          var numberOfCrossovers : Integer);
      { creates a child creature from the genetic material of two parents }

      constructor Create (var varyPar   : TBVector;
                          fixedPar      : TVectorFn;
                          var minFitPar : TRVector;
                          var maxFitPar : TRVector;
                          sigFigures    : Integer);
    end;

implementation

  uses
    Classes,       // for constructing objects
    SysUtils,      // exception handler
    gaChromosomes; // genetic algorithm's chromosome type

  { public routines }

  function TCreature.IsEqualTo (c : TCreature) : Boolean;
  begin
    IsEqualTo := genome.IsEqualTo(c.genome)
  end;
  
  function TCreature.Copy : TCreature;
    var
      c : TCreature;
  begin
    c := TCreature.Create(varyThisParameter, fixedParameterFn, 
                          minFitParameters, maxFitParameters, 
                          significantFigures);
    c.genome     := genome.Copy;
    c.quality    := quality;
    c.rmseValues := CopyRData(rmseValues);
    Copy := c
  end;

  procedure TCreature.GetAllParameters (var parameters : TRVector);
    var
      chromosome : TChromosome;
      i, m, n    : Integer;
  begin
    if (parameters = nil) or (LenRVector(parameters) <> dimP) then
      parameters := NewRVector(dimP); 
    m := 1;
    n := 1;
    for i := 1 to dimP do
      if varyThisParameter[i] then begin
        chromosome    := genome.Pop(m);
        parameters[i] := chromosome.Decode;
        chromosome.DeleteGenes;
        Inc(m)
      end
      else begin
        parameters[i] := fixedParameters[n];
        Inc(n)
      end
  end;
      
  procedure TCreature.GetGenome (var genomes : TSVector);
    var
      i, m    : Integer;
      sVector : TSVector;
  begin
    sVector := genome.Print;
    genomes := NewSVector(dimP);
    m := 1;
    for i := 1 to dimP do 
      if not varyThisParameter[i] then
        genomes[i] := ''
      else begin
        genomes[i] := sVector[m];
        Inc(m)
      end
  end;

  procedure TCreature.GetGenetics (out chromosomes, genes : Integer);
  begin
    chromosomes := genome.Chromosomes;
    genes       := genome.Genes
  end;

  procedure TCreature.GetProbabilities (out mutation, crossover : Real);
  begin
    genome.GetProbabilities(mutation, crossover)
  end;

  procedure TCreature.SetFitness (fitness : Real);
  begin
    quality := fitness
  end;

  function TCreature.GetFitness : Real;
  begin
    GetFitness := quality
  end;

  procedure TCreature.SetRMSE (var rmseMtx : TRData);
  begin
    rmseValues := CopyRData(rmseMtx)
  end;

  function TCreature.GetRMSE (experiment       : Integer;
                              responseVariable : Integer) : Real;
  begin
    GetRMSE := rmseValues[experiment][responseVariable][1]
  end;
  
  procedure TCreature.DeleteGenome;
  begin
    DelRVector(fixedParameters);
    DelRVector(maxFitParameters);
    DelRVector(minFitParameters);
    DelBVector(varyThisParameter);
    DelRData(rmseValues);
    if genome <> nil then
      genome.DeleteChromosomes;
    genome := nil
  end;

  procedure TCreature.AssignFixedParameters;
    var
      fitParameters : TRVector;
  begin
    if fixedParameterFn = nil then
      fixedParameters := NewRVector(0)
    else begin
      fitParameters   := genome.Decode; 
      fixedParameters := fixedParameterFn(fitParameters)
    end
  end;
      
  procedure TCreature.Procreate;
  begin
    genome := TGenome.Create(minFitParameters, maxFitParameters, 
                             significantFigures);
    AssignFixedParameters
  end;

  procedure TCreature.Alien (var fitParameters : TRVector);
    var
      c          : Integer;
      chromosome : TChromosome;
  begin
    if (fitParameters = nil) or (LenRVector(fitParameters) = 0) then raise 
      Exception.Create('supplied parameter vector had not been created');
    if LenRVector(fitParameters) <> dimPV then raise Exception.Create
      ('alien parameters must have length of # varied parameters');
    genome := TGenome.Create(minFitParameters, maxFitParameters,
                             significantFigures);
    for c := 1 to dimPV do begin
      chromosome := genome.Pop(c);
      chromosome.Encode(fitParameters[c]);
      genome.Put(c, chromosome)
    end;
    AssignFixedParameters;
    // clean up
    chromosome.DeleteGenes
  end;

  procedure TCreature.Conceive (parentA, parentB       : TCreature;
                                var numberOfMutations  : Integer;
                                var numberOfCrossovers : Integer);
  begin
    if (parentA = nil) or (parentB = nil) then raise Exception.Create
      ('a supplied parent had not been created and cannot conceive');
    gaGenome.GeneticRecombination(parentA.genome, parentB.genome,
                                  numberOfCrossovers, genome);
    genome.Mutate(numberOfMutations);
    AssignFixedParameters
  end;

  { the following procedure is the class constructor }

  constructor TCreature.Create (var varyPar   : TBVector;
                                fixedPar      : TVectorFn;
                                var minFitPar : TRVector;
                                var maxFitPar : TRVector;
                                sigFigures    : Integer);
  begin
    if ((varyPar = nil)   or (LenBVector(varyPar) = 0))   or
       ((minFitPar = nil) or (LenRVector(minFitPar) = 0)) or 
       ((maxFitPar = nil) or (LenRVector(maxFitPar) = 0)) then raise
      Exception.Create
        ('a supplied vector for initialization had not been created');
    if LenRVector(maxFitPar) <> LenRVector(minFitPar) then raise 
      Exception.Create
        ('lower and upper parameter bounds must have equal lengths');
    dimP  := LenBVector(varyPar);
    dimPV := LenRVector(minFitPar);
    varyThisParameter := CopyBVector(varyPar);
    maxFitParameters  := CopyRVector(maxFitPar);
    minFitParameters  := CopyRVector(minFitPar);
    fixedParameterFn  := fixedPar;
    if sigFigures < 1 then
      significantFigures := 1
    else if sigFigures > maxSignificantFigures then
      significantFigures := maxSignificantFigures
    else
      significantFigures := sigFigures
  end;

end.
