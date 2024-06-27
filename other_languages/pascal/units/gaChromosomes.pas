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
   Port began on - February 27, 2014
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
   This unit introduces the chromosome - a container of genes.
***************************************************************************** *)

// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

unit gaChromosomes;

interface

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    gaCore,    // genetic algorithm's core unit where types are defined
    gaGenes;   // genetic algorithm's gene class

  type
    TChromosome = class(TObject)
    private
      binary       : TBVector;
      chromosome   : array of TGene;
      gray         : TBVector;
      maxPhenotype : Real;
      minPhenotype : Real;
      rangeGenes   : Real;
      sigFigures   : Integer;
      procedure BinaryToGray;
      procedure GrayToBinary;
      function  BinaryToInteger : Int64;
      procedure IntegerToBinary (i : Int64);
      function  IntegerToPhenotype (i : Int64) : Real;
      function  PhenotypeToInteger (p : Real) : Int64;
    public
      function IsEqualTo (c : TChromosome) : Boolean;
      { compares two chromosomes to see if their genes are equivalent }

      function Print : String;
      { writes a chromosome to string as a concatination of gene values }

      function Copy : TChromosome;
      { returns a deep copy of the chromosome }

      procedure GetPhenotypeInfo (out minParameter, maxParameter : Real;
                                  out numberOfSignificantFigures : Integer);
      { returns the phenotype information stored by the chromosome }

      function Genes : Integer;
      { returns the number of individual genes in the chromosome }
      
      procedure DeleteGenes;
      { removes all genes from the chromosome }

      { chromosomes POP and PUSH genes index from 1 }

      function Pop (gene : Integer) : TGene;
      { extracts a gene from the chromosome at the specified location }

      procedure Put (gene : Integer; g : TGene);
      { inserts a gene into the chromosome at the specified location }

      procedure Mutate (probabilityOfMutation : Real;
                        var numberOfMutations : Integer);
      { assess whether or not to mutate each gene at specified probability,
        and if there are mutatations, increment the counter accordingly }

      function Decode : Real; 
      { extracts the phenotype held in the genetic code of the chromosome }

      procedure Encode (phenotype : Real);
      { inserts a pheontype into the genetic code of the chromosome }

      constructor Create (minParameter, maxParameter : Real;
                          numberOfSignificantFigures : Integer);
    end;

  procedure GeneticRecombination (probabilityOfCrossover : Real;
                                  parentA, parentB       : TChromosome;
                                  var numberOfCrossovers : Integer;
                                  out child              : TChromosome);
  { This takes two parent chromosomes and gives birth to a child chromosome
    where splitting and crossover occur at a probabilityOfCrossover, and
    if crossover occurs, then the associated counter is incremented.
    The location of splitting occurs at a random location along it length.}

implementation

  uses
    Classes,         // allows creation of objects
    SysUtils,        // the system's exception handler
    Math,            // basic math package
    gaProbabilities; // genetic algorithm's probability functions

  { private routines }

  { The encode/decode maps between a real and its Haploid representation. }

  { The algorithm for converting between binary and gray codes assumes the
    most significant bit (MSB) is at the left of the code, and associates
    with the [1] position in the binary and gray encodings.  The least
    significant bit associates with position [High], in other words, e.g.,
       code = [1|0|1|1|0|0|1|0] has a MSB of 1 and a LSB of 0. }

  procedure TChromosome.BinaryToGray;
    var
      i : Integer;
  begin
    gray[1] := binary[1];
    for i := 2 to Genes do
      gray[i] := binary[i - 1] xor binary[i]
  end;

  procedure TChromosome.GrayToBinary;
    var
      i : Integer;
  begin
    binary[1] := gray[1];
    for i := 2 to Genes do
      binary[i] := binary[i - 1] xor gray[i]
  end;

  function TChromosome.BinaryToInteger : Int64;
    var
      i     : Integer;
      int   : Int64;
      power : Int64;
  begin
    int   := 0;
    power := 1;
    for i := Genes downto 1 do begin
      if binary[i] then
        int := int + power;
      power := 2 * power
    end;
    BinaryToInteger := int
  end;

  procedure TChromosome.IntegerToBinary (i : Int64);
    var
      int  : Int64;
      j, k : Integer;
  begin
    j   := Genes;
    int := i;
    while int > 0 do begin
      if (int mod 2) = 0 then
        binary[j] := False
      else
        binary[j] := True;
      int := int div 2;
      Dec(j)
    end;
    // remaining higher-order binary bits are zeros
    for k := j downto 1 do
      binary[k] := False
  end;

  function TChromosome.IntegerToPhenotype (i : Int64) : Real;
    var
      phenotype : Real;
  begin
    phenotype := minPhenotype 
               + (i / rangeGenes) * (maxPhenotype - minPhenotype);
    if phenotype < minPhenotype then
      phenotype := minPhenotype;
    if phenotype > maxPhenotype then
      phenotype := maxPhenotype;
    IntegerToPhenotype := phenotype
  end;

  function TChromosome.PhenotypeToInteger (p : Real) : Int64;
    var
      fraction : Real;
      intValue : Int64;
      maxGene  : Int64;
  begin
    fraction := (p - minPhenotype) / (maxPhenotype - minPhenotype);
    intValue := Round(fraction * rangeGenes);
    if intValue < 0 then
      intValue := 0;
    maxGene := Round(rangeGenes);
    if intValue > maxGene then
      intValue := maxGene;
    PhenotypeToInteger := intValue
  end;

  { public routines }

  function TChromosome.IsEqualTo (c : TChromosome) : Boolean;
    var
      geneL : TGene;
      geneR : TGene;
      i     : Integer;
      truth : Boolean;
  begin
    if Genes = 0 then
      if (c = nil) or (c.Genes = 0) then
        truth := True
      else
        truth := False
    else begin
      i := 1;
      repeat
        geneL := chromosome[i].Copy;
        geneR := c.chromosome[i].Copy;
        if not geneL.IsEqualTo(geneR) then
          Exit(False);
        Inc(i)
      until i > Genes;
      truth := True
    end;
    // clean up
    geneL := nil;
    geneR := nil;
    // return the result
    IsEqualTo := truth
  end;

  function TChromosome.Print : String;
    var
      i : Integer;
      s : String;
  begin
    s := '';
    if Genes > 0 then
      for i := 1 to Genes do
        s := s + chromosome[i].Print;
    Print := s
  end;

  function TChromosome.Copy : TChromosome;
    var
      c : TChromosome;
      i : Integer;
  begin
    c := TChromosome.Create(minPhenotype, maxPhenotype, sigFigures);
    for i := 1 to Genes do 
      c.chromosome[i] := chromosome[i].Copy;
    c.maxPhenotype := maxPhenotype;
    c.minPhenotype := minPhenotype;
    c.rangeGenes   := rangeGenes;
    c.sigFigures   := sigFigures;
    // return the result
    Copy := c
  end;

  procedure TChromosome.GetPhenotypeInfo
                         (out minParameter, maxParameter : Real;
                          out numberOfSignificantFigures : Integer);
  begin
     minParameter := minPhenotype;
     maxParameter := maxPhenotype;
     numberOfSignificantFigures := sigFigures
  end;

  function TChromosome.Genes : Integer;
  begin
    if (chromosome = nil) or (Length(chromosome) < 2) then
      Genes := 0
    else
      Genes := High(chromosome)
  end;
  
  procedure TChromosome.DeleteGenes;
    var
      i : Integer;
  begin
    DelBVector(binary);
    DelBVector(gray);
    for i := Genes downto 0 do
      chromosome[i] := nil;
    SetLength(chromosome, 0)
  end;

  function TChromosome.Pop (gene : Integer) : TGene;
  begin
    if Genes = 0 then raise 
      Exception.Create('cannot pop a gene from an empty chromosome');
    if (gene < 1) or (gene > Genes) then raise
      Exception.Create('pop-from location is out of range');
    Pop := chromosome[gene].Copy
  end;

  procedure TChromosome.Put (gene : Integer; g : TGene);
  begin
    if Genes = 0 then raise 
      Exception.Create('cannot push a gene into an empty chromosome');
    if g = nil then raise
      Exception.Create('cannot push an empty gene into a chromosome');
    if (gene < 1) or (gene > Genes) then raise
      Exception.Create('push-to location is out of range');
    chromosome[gene] := g.Copy
  end;

  procedure TChromosome.Mutate (probabilityOfMutation : Real;
                                var numberOfMutations : Integer);
    var
      g : TGene;
      i : Integer;
      n : Integer;
  begin
    if Genes > 0 then
      for i := 1 to Genes do begin
        g := Pop(i);
        n := numberOfMutations;
        g.Mutate(probabilityOfMutation, numberOfMutations);
        if numberOfMutations > n then
          Put(i, g)
      end;
    // clean up
    g := nil
  end;

  function TChromosome.Decode : Real;
    var
      allele    : THaploid;
      i         : Integer;
      int       : Int64;
      phenotype : Real;
  begin
    if Genes = 0 then 
      phenotype := 0.0
    else begin
      for i := 1 to Genes do begin
        allele := chromosome[i].Pop;
        if allele = dominant then
          gray[i] := True
        else // recessive
          gray[i] := False
      end;
      GrayToBinary;
      int := BinaryToInteger;
      phenotype := IntegerToPhenotype(int)
    end;
    // return the result
    Decode := phenotype
  end;

  procedure TChromosome.Encode (phenotype : Real);
    var
      i   : Integer;
      int : Int64;
  begin
    if Genes > 0 then begin
      int := PhenotypeToInteger(phenotype);
      IntegerToBinary(int);
      BinaryToGray;
      for i := 1 to Genes do
        if gray[i] = True then 
          chromosome[i] := TGene.CreateWithGene(dominant)
        else
          chromosome[i] := TGene.CreateWithGene(recessive);
    end
  end;

  constructor TChromosome.Create (minParameter, maxParameter : Real;
                                  numberOfSignificantFigures : Integer);
    var
      bits          : Integer;
      bitsPerDecade : Integer;
      decades, i    : Integer;
      logDecades    : Real;
  begin
    if minParameter >= maxParameter then raise Exception.Create
      ('cannot create a chromosome unless maxParameter > minParameter');
    if minParameter > 0.0 then
      logDecades := Log10(maxParameter / minParameter)
    else if maxParameter < 0.0 then
      logDecades := Log10(minParameter / maxParameter)
    else if minParameter = 0.0 then
      logDecades := Log10(maxParameter)
    else { maxParameter = 0.0 }
      logDecades := Log10(-minParameter);
    if logDecades > 0.0 then
      decades := Ceil(logDecades)
    else
      decades := Abs(Floor(logDecades));
    if decades < 1 then
      decades := 1;
    if decades > 9 then
      decades := 9;
    case numberOfSignificantFigures of
    0..1 : bitsPerDecade := 7;
    2    : bitsPerDecade := 10;
    3    : bitsPerDecade := 14;
    4    : bitsPerDecade := 17;
    5    : bitsPerDecade := 20;
    6    : bitsPerDecade := 24;
    else   bitsPerDecade := 27
    end;
    bits := bitsPerDecade * decades;
    if bits > 62 then
      bits := 62;             // 63 will allow overflows to occur when rounding
    try
      SetLength(chromosome, Succ(bits));
      chromosome[0] := nil    // element [0] is not used
    except on E : EOutOfMemory do
      WriteLn('Memory error. Details: ' + E.ClassName + '/' + E.Message)
    end;
    binary := NewBVector(bits);     
    gray   := NewBVector(bits);    
    for i := 1 to bits do 
      Put(i, TGene.Create);
    minPhenotype := minParameter;
    maxPhenotype := maxParameter;
    rangeGenes   := IntPower(2.0, bits);
    sigFigures   := numberOfSignificantFigures
  end;

  { end of method declarations for TChromosome }

  procedure GeneticRecombination (probabilityOfCrossover : Real;
                                  parentA, parentB       : TChromosome;
                                  var numberOfCrossovers : Integer;
                                  out child              : TChromosome);
    var
      i     : Integer;
      xover : Integer;
  begin
    if (parentA = nil) or (parentB = nil) then raise
      Exception.Create('X over requires two parents, one or both were nil');
    if parentA.Genes <> parentB.Genes then raise
      Exception.Create('chromosome parents do not belong to the same species');
    if IsHeads(evenOdds) then begin
      // left side of chromosome belongs to parent A, the right to parent B
      child := parentA.Copy;
      if IsHeads(probabilityOfCrossover) then begin
        xover := RandomInteger(2, parentA.Genes - 1);
        for i := xover to parentA.Genes do
          child.Put(i, parentB.Pop(i));
        Inc(numberOfCrossovers)
      end
    end
    else begin
      // left side of chromosome belongs to parent B, the right to parent A
      child := parentB.Copy;
      if IsHeads(probabilityOfCrossover) then begin
        xover := RandomInteger(2, parentB.Genes - 1);
        for i := xover to parentB.Genes do
          child.Put(i, parentA.Pop(i));
        Inc(numberOfCrossovers)
      end
    end
  end;

end.
