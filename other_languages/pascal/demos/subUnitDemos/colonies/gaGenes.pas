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
   Port began on - February 26, 2014
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
   This unit introduces the notion of a gene, in particular, a haploid gene.
***************************************************************************** *)

// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

unit gaGenes;

interface

  const
    alphabet = 2.0;  { a haploid gene has two gene expressions }

  type
    THaploid = (recessive, dominant);

  type
    TGene = class(TObject)
    private
      allele : THaploid;
    public
      function IsEqualTo (g : TGene) : Boolean;
      { compares two genes to see if their allele have equivalent values }

      function Print : String;
      { writes a gene in a string format: 0 if recessive, 1 if dominant }

      function Copy : TGene;
      { returns a deep copy of the gene }

      procedure Put (haploid : THaploid);
      { assigns a haploid value to the gene }

      function Pop : THaploid;
      { returns the haploid value held by the gene }

      procedure Mutate (probabilityOfMutation : Real;
                        var numberOfMutations : Integer);
      { assesses whether or not to mutate at specified probability
        and if it mutates it increments the returned counter }

      constructor Create;
      constructor CreateWithGene (haploid : THaploid);
    end;

implementation

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    Classes,         // allows creation of objects
    SysUtils,        // the system's exception handler
    gaCore,          // genetic algorithm's core routines
    gaProbabilities; // genetic algorithm's probability functions

  function TGene.IsEqualTo (g : TGene) : Boolean;
  begin
    if ((allele = recessive) and (g.allele = recessive)) or
       ((allele = dominant)  and (g.allele = dominant))  then
      IsEqualTo := True
    else
      IsEqualTo := False
  end;

  function TGene.Print : String;
  begin
    if allele = recessive then
      Print := '0'
    else { allele = dominant }
      Print := '1'
  end;

  function TGene.Copy : TGene;
  begin
    Copy := TGene.CreateWithGene(allele)
  end;

  procedure TGene.Put (haploid : THaploid);
  begin
    if haploid = recessive then
      allele := recessive
    else
      allele := dominant
  end;

  function TGene.Pop : THaploid;
  begin
    if allele = recessive then
      Pop := recessive
    else
      Pop := dominant
  end;

  procedure TGene.Mutate (probabilityOfMutation : Real;
                          var numberOfMutations : Integer);
  begin
    if IsHeads(probabilityOfMutation) then
      begin
        if allele = recessive then
          allele := dominant
        else
          allele := recessive;
        Inc(numberOfMutations)
      end
  end;

  constructor TGene.Create;
  begin
    if IsHeads(evenOdds) then
      allele := dominant
    else
      allele := recessive
  end;

  constructor TGene.CreateWithGene (haploid : THaploid);
  begin
    allele := haploid
  end;

end.
