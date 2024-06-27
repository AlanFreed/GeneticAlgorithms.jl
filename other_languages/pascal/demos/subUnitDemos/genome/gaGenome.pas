(* *****************************************************************************
   Author         - Alan D. Freed, Ph.D
   License        - GNU Lesser General Public License, vs. 3 or later
   Copyright      - (c) Alan D. Freed 2014
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
   Port began on - March    16, 2014
   Last modified - November 21, 2015
--------------------------------------------------------------------------------
   Version 1.2 - 1.3
   Adopted 'Arrays' for creating and managing vector arrays.  
   This library calls librarys 'Arrays' and 'rng' written by the author.
--------------------------------------------------------------------------------
   This unit belongs to the genetic part of this genetic algorithm, including:
      gaProbabilities - the probability functions utilized
      gaGenes         - defines a gene class based upon a haploid gene
      gaChromosomes   - defines a chromosome class as a container of genes
      gaGenome        - defines a genome class as a container of chromosomes
--------------------------------------------------------------------------------
   This unit introduces the genome - the DNA of a creature.
***************************************************************************** *)

// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

unit gaGenome;

interface

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    gaCore,          // basic types used in the genetic algorithm
    gaChromosomes;   // genetic algorithm's chromosome class

  type
    TGenome = class
    private
      genome     : array of TChromosome;
      probMutate : Real;
      probXover  : Real;
    public
      function IsEqualTo (g : TGenome) : Boolean;
      { compares two genome to see if their genes are equivalent }

      function Copy : TGenome;
      { returns a deep copy of the genome }

      function Print : TSVector;
      { writes the genome into an array of strings of chromosomes }

      function Genes : Integer;
      { returns the number of individual genes in the genome }

      function Chromosomes : Integer;
      { returns the number of individual chromosomes in the genome }
      
      procedure DeleteChromosomes;
      { removes all chromosomes from the genome }

      procedure GetProbabilities (out mutation, crossover : Real);
      { returns probabilities shared by all chromosomes & genes therein }

      { genome POP and PUT chromosomes index from 1 }

      function Pop (chromosome : Integer) : TChromosome;
      { extracts a chromosome from the genome at the specified location }

      procedure Put (chromosome : Integer; c : TChromosome);
      { inserts a chromosome into the genome at the specified location }

      procedure Mutate (var numberOfMutations : Integer);
      { assesses whether or not to mutate each gene, and if there are,
        the number of mutatations increment the counter accordingly }

      { all procedures that deal with parameter vectors index from 1 }

      function Decode : TRVector;
      { extract the phenotypes held in the genetic code of the genome }

      //  argument arrays are passed by reference

      procedure Encode (var phenotypes : TRVector);
      { imsert pheontypes into the genetic code of the genome }

      constructor Create (var minParameters, maxParameters : TRVector;
                          numberOfSignificantFigures       : Integer);
      { the two parameter vectors index from 1 }
    end;

  procedure GeneticRecombination(parentA, parentB       : TGenome;
                                 var numberOfCrossovers : Integer;
                                 out child              : TGenome);
  { This takes two parent genomes and gives birth to a child genome,
    adjusting the counter for crossovers accordingly.  The location of
    chromosome splittings occur at random locations along their lengths.}

implementation

  uses
    Classes,         // allows creation of objects
    SysUtils,        // exception handler
    gaProbabilities; // genetic algorithm's probability functions

  { probabilities are not tested for equality, only the allele of the DNA }
  function TGenome.IsEqualTo (g : TGenome) : Boolean;
    var
      chromosomeL : TChromosome;
      chromosomeR : TChromosome;
      i           : Integer;
      truth       : Boolean;
  begin
    if Chromosomes = 0 then
      if (g = nil) or (g.Chromosomes = 0) then
        truth := True
      else
        truth := False
    else begin
      i  := 1;
      repeat
        chromosomeL := genome[i].Copy;
        chromosomeR := g.genome[i].Copy;
        if not chromosomeL.IsEqualTo(chromosomeR) then
          Exit(False);
        Inc(i)
      until i > Chromosomes;
      truth := True
    end;
    // clean up
    chromosomeL.DeleteGenes;
    chromosomeR.DeleteGenes;
    // return the result
    IsEqualTo := truth
  end;

  { write the chromosomes into a string vector starting at location [1] }
  function TGenome.Print : TSVector;
    var
      i    : Integer;
      sVec : TSVector;
  begin
    sVec := NewSVector(Chromosomes);
    for i := 1 to Chromosomes do
      sVec[i] := genome[i].Print;
    Print := sVec
  end;

  { probabilities are not copied, but assigned anew and randomly }
  function TGenome.Copy : TGenome;
    var
      g    : TGenome;
      i    : Integer;
      maxP : TRVector;
      minP : TRVector;
      sigF : Integer;
  begin
    maxP := NewRVector(Chromosomes);  
    minP := NewRVector(Chromosomes);  
    sigF := 0;
    // make an exact copy of the chromosome information
    for i := 1 to Chromosomes do
      genome[i].GetPhenotypeInfo(minP[i], maxP[i], sigF);
    g := TGenome.Create(minP, maxP, sigF);
    for i := 1 to Chromosomes do
      g.genome[i] := genome[i].Copy;
    // the probailities for mutation and crossover are not duplicated
    g.probMutate := RandomProbability(meanProbabilityOfMutation,
                                      standardDeviationOfMutation);
    g.probXover  := RandomProbability(meanProbabilityOfCrossover,
                                      standardDeviationOfCrossover);
    // clean up
    DelRVector(maxP);
    DelRVector(minP);
    // return the result
    Copy := g
  end;

  function TGenome.Genes : Integer;
    var
      c, g : Integer;
  begin
    g := 0;
    for c := 1 to Chromosomes do
      g := g + genome[c].Genes;
    // return the result
    Genes := g
  end;

  function TGenome.Chromosomes : Integer;
  begin
    if (genome = nil) or (Length(genome) < 2) then
      Chromosomes := 0
    else
      Chromosomes := High(genome)
  end;
  
  procedure TGenome.DeleteChromosomes;
    var
      c : Integer;
  begin
    for c := Chromosomes downto 1 do
      genome[c].DeleteGenes;
    for c := Chromosomes downto 0 do
      genome[c] := nil;
    SetLength(genome, 0)
  end;

  procedure TGenome.GetProbabilities (out mutation, crossover : Real);
  begin
    mutation  := probMutate;
    crossover := probXover
  end;

  function TGenome.Pop (chromosome : Integer) : TChromosome;
  begin
    if Chromosomes = 0 then raise 
      Exception.Create('cannot pop a chromosome from an empty genome');
    if (chromosome < 1) or (chromosome > Length(genome)) then raise
      Exception.Create('pop-from location is out of range');
    Pop := genome[chromosome].Copy
  end;

  procedure TGenome.Put (chromosome : Integer; c : TChromosome);
  begin
    if Chromosomes = 0 then raise 
      Exception.Create('cannot push a chromosome into an empty genome');
    if (c = nil) or (c.Genes = 0) then raise 
      Exception.Create('cannot push an empty chromosome into a genome');
    if (genome[chromosome].Genes <> c.Genes) then raise
      Exception.Create('cannot push a chromosome with the wrong # genes');
    if (chromosome < 1) or (chromosome > Chromosomes) then raise
      Exception.Create('push-to location is out of range');
    genome[chromosome] := c.Copy
  end;

  procedure TGenome.Mutate (var numberOfMutations : Integer);
    var
      c : TChromosome;
      i : Integer;
      n : Integer;
  begin
    if Chromosomes > 0 then
      for i := 1 to Chromosomes do begin
        c := Pop(i);
        n := numberOfMutations;
        c.Mutate(probMutate, numberOfMutations);
        if numberOfMutations > n then
          Put(i, c);
      end;
    // clean up
    c.DeleteGenes
  end;

  function TGenome.Decode : TRVector;
    var
      i          : Integer;
      phenotypes : TRVector;
  begin
    phenotypes := NewRVector(Chromosomes);
    for i := 1 to Chromosomes do 
      phenotypes[i] := genome[i].Decode;
    // return the result
    Decode := phenotypes
  end;

  procedure TGenome.Encode (var phenotypes : TRVector);
    var
      i : Integer;
  begin
    if (phenotypes = nil) or (Chromosomes <> LenRVector(phenotypes)) then raise
      Exception.Create('# of phenotypes to encode must equal # of chromosomes');
    for i := 1 to Chromosomes do 
      genome[i].Encode(phenotypes[i])
  end;

  { probabilities are assigned randomly for each creation }
  constructor TGenome.Create (var minParameters, maxParameters : TRVector;
                              numberOfSignificantFigures       : Integer);
    var
      i : Integer;
  begin
    if (minParameters = nil) or (maxParameters = nil) then raise
      Exception.Create('a parameter array sent to the constructor was nil');
    if Length(minParameters) <> Length(maxParameters) then raise
      Exception.Create('lengths of minParamters and maxParameters must equal');
    try
      SetLength(genome, LenRVector(minParameters) + 1);
      genome[0] := nil   // element [0] is not used
    except on E : EOutOfMemory do
      begin
        WriteLn('Memory error. Details: ' + E.ClassName + '/' + E.Message);
        Exit
      end
    end;
    for i := 1 to LenRVector(minParameters) do
      genome[i] := TChromosome.Create(minParameters[i], maxParameters[i],
                                      numberOfSignificantFigures);
    probMutate := RandomProbability(meanProbabilityOfMutation,
                                    standardDeviationOfMutation);
    probXover  := RandomProbability(meanProbabilityOfCrossover,
                                    standardDeviationOfCrossover)
  end;

  { end of method declarations for TGenome }

  procedure GeneticRecombination (parentA, parentB       : TGenome;
                                  var numberOfCrossovers : Integer;
                                  out child              : TGenome);
    var
      childChromosome : TChromosome;
      chromosomeA     : TChromosome;
      chromosomeB     : TChromosome;
      i               : Integer;
      nCrossovers     : Integer;
      probMutate      : Real;
      probXover       : Real;
  begin
    if (parentA = nil) or (parentB = nil) then raise
      Exception.Create('X over requires two parents, one or both were nil');
    if (parentA.Chromosomes <> parentB.Chromosomes) then raise
      Exception.Create('parentA and parentB have incompatible genomes');
    child := parentA.Copy;
    probMutate := 0.0;
    probXover  := 0.0;
    child.GetProbabilities(probMutate, probXover);
    for i := 1 to parentA.Chromosomes do begin
      chromosomeA := parentA.Pop(i);
      chromosomeB := parentB.Pop(i);
      nCrossovers := numberOfCrossovers;
      childChromosome := nil;
      gaChromosomes.GeneticRecombination(probXover,
                                         chromosomeA, chromosomeB,
                                         numberOfCrossovers,
                                         childChromosome);
      if numberOfCrossovers > nCrossovers then
        child.Put(i, childChromosome)
      else if IsHeads(evenOdds) then
        child.Put(i, chromosomeA)
      else
        child.Put(i, chromosomeB)
    end;
    // clean up
    childChromosome.DeleteGenes;
    chromosomeA.DeleteGenes;
    chromosomeB.DeleteGenes
  end;

end.
