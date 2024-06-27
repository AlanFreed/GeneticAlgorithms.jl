// global directives
{$MODE OBJFPC}    // allows objects, classes, interfaces and exception handling
{$H+}             // use ANSI strings
{$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
{$PIC ON}         // create Position Independent Code - for Unix libraries

program testGenome;

uses
   gaCore, gaChromosomes, gaGenome;

  function NewRVector (len : Integer) : TRVector; external 'Arrays';
  // Creates a real vector indexing from 1 to len whose elements are all 0.0.

   procedure WriteGenome (var g : TGenome);
   var
      grayGenes  : TSVector;
      i          : Integer;
      phenotypes : TRVector;
   begin
      grayGenes  := g.Print;
      phenotypes := g.Decode;
      for i := 1 to g.Chromosomes do begin
         Write('Chromosome ', i, ' has values ');
         WriteLn(grayGenes[i], ' <=> ', phenotypes[i]:6:4)
      end;
   end;

   procedure Test;
   var
      genome1 : TGenome;
      genome2 : TGenome;
      genome3 : TGenome;
      i, n    : Integer;
      minX    : TRVector;
      maxX    : TRVector;
      probMut : Real;
      probXOv : Real;
      sigFig  : Integer;
   begin
      probMut := 0.0;
      probXOv := 0.0;
      WriteLn('Test the genome class.');
      WriteLn;
      Write('How many chromosomes are there to be? ');
      ReadLn(n);
      minX := NewRVector(n);
      maxX := NewRVector(n);
      for i := 1 to n do
         begin
            WriteLn;
            Write('Minimum real/phenotype value for chromosome ', i, ': ');
            ReadLn(minX[i]);
            Write('Maximum real/phenotype value for chromosome ', i, ': ');
            ReadLn(maxX[i])
         end;
      WriteLn;
      Write('Enter in the number of significant figures: ');
      ReadLn(sigFig);
      WriteLn;
      WriteLn('The first random genome has chromosomes');
      genome1 := TGenome.Create(minX, maxX, sigFig);
      WriteGenome(genome1);
      WriteLn('The second random genomoe has chromosomes');
      genome2 := TGenome.Create(minX, maxX, sigFig);
      writeGenome(genome2);
      WriteLn('A copy of the second genome has chromosomes');
      genome3 := genome2.Copy;
      WriteGenome(genome3);
      WriteLn;
      WriteLn('Test for equality:');
      Write('   this test should yield TRUE:  ');
      if genome2.IsEqualTo(genome3) then
         WriteLn('True')
      else
         WriteLn('False');
      Write('   this test should yield FALSE: ');
      if genome1.IsEqualTo(genome2) then
         WriteLn('True')
      else
         WriteLn('False');
      WriteLn;
      WriteLn('The first genome has chromosomes:');
      WriteGenome(genome1);
      WriteLn('which is comprised of ', genome1.Genes, ' genes.');
      WriteLn;
      WriteLn('The second genome has chromosomes:');
      WriteGenome(genome2);
      WriteLn('which is comprised of ', genome2.Genes, ' genes.');
      WriteLn;
      WriteLn('A mutation of this genome has genes ');
      n := 0;
      genome2.Mutate(n);
      WriteGenome(genome2);
      Write('incurring ', n);
      if n <> 1 then
         Write(' mutations ')
      else
         Write(' mutation ');
      genome2.GetProbabilities(probMut, probXOv);
      WriteLn('at a probability of ', probMut:6:4, '.');
      WriteLn;
      n := 0;
      gaGenome.GeneticRecombination(genome1, genome2, n, genome3);
      WriteLn('The first and second genomes mated to produce the child');
      WriteGenome(genome3);
      Write('incurring ', n);
      if n <> 1 then
         Write(' crossovers ')
      else
         Write(' crossover ');
      WriteLn('at a probability of ', probXOv:6:4, '.')
   end;

begin
  Test
end.

