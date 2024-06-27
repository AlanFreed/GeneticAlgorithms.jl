// global directives
{$MODE OBJFPC}    // allows objects, classes, interfaces and exception handling
{$H+}             // use ANSI strings
{$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
{$PIC ON}         // create Position Independent Code - for Unix libraries

program testChromosomes;

uses
   gaCore, gaChromosomes;

   procedure Test;
   var
      child  : TChromosome;
      chrom1 : TChromosome;
      chrom2 : TChromosome;
      chrom3 : TChromosome;
      minX   : Real;
      maxX   : Real;
      n      : Integer;
      sigFig : Integer;
      x      : Real;
   begin
      WriteLn('Test the chromosome class.');
      WriteLn;
      Write('Enter in the minimum real/phenotype value:  ');
      ReadLn(minX);
      Write('Enter in the maximum real/phenotype value:  ');
      ReadLn(maxX);
      Write('Enter in the number of significant figures: ');
      ReadLn(sigFig);
      WriteLn('Chromosome    Phenotype    Gene Expression');
      chrom1 := TChromosome.Create(minX, maxX, sigFig);
      x := chrom1.Decode;
      WriteLn('    1         ', x:8:4, '      ', chrom1.Print);
      chrom2 := TChromosome.Create(minX, maxX, sigFig);
      x := chrom2.Decode;
      WriteLn('    2         ', x:8:4, '      ', chrom2.Print);
      WriteLn;
      Write('Chromosome 1');
      if chrom1.IsEqualTo(chrom2) then 
        Write(' = ')
      else
        Write(' # ');
      WriteLn('chromosome 2');
      WriteLn;
      Write('Enter a real number to encode/decode: ');
      chrom3 := TChromosome.Create(minX, maxX, sigFig);
      ReadLn(x);
      if (x < minX) or (x > maxX) then
         begin
           WriteLn('   The value must lie within [', minX, ', ', maxX, '].');
           Write('   Try again.  Enter a number to encode/decode: ');
           ReadLn(x)
         end;
      chrom3.Encode(x);
      WriteLn('Gray encoding of this real is:    ' + chrom3.Print);
      x := chrom3.Decode;
      WriteLn('and it decodes back into:         ', x:8:4);
      n := 0;
      child := chrom2.Copy;
      x := child.Decode;
      WriteLn;
      WriteLn('A copy of chromosome 2 has data:');
      WriteLn('    copy      ', x:8:4, '      ', child.Print);
      Write('The copied chromosome');
      if child.IsEqualTo(chrom2) then 
        Write(' = ')
      else
        Write(' # ');
      WriteLn('chromosome 2');
      WriteLn;
      GeneticRecombination(0.8, chrom1, chrom2, n, child);
      x := child.Decode;
      WriteLn('A Xover with chromosomes 1 & 2 produces:');
      WriteLn('    3         ', x:8:4, '      ', child.Print);
      WriteLn('in ', n, ' crossovers');
      // clean up
      if child <> nil then
        child.Free;
      if chrom1 <> nil then
        chrom1.Free;
      if chrom2 <> nil then
        chrom2.Free;
      if chrom3 <> nil then
        chrom3.Free
   end;

begin
  Test
end.

