// global directives
{$MODE OBJFPC}    // allows objects, classes, interfaces and exception handling
{$H+}             // use ANSI strings
{$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
{$PIC ON}         // create Position Independent Code - for Unix libraries

program testGenes;

uses
    gaCore, gaGenes;

   procedure Test;
   var
      count       : Integer;
      gene1       : TGene;
      gene2       : TGene;
      haploid     : THaploid;
      nbrMutation : Integer;
      probability : Real;
   begin
      WriteLn('This demo tests the gene class.');
      WriteLn;
      gene1 := TGene.Create;
      if gene1 = nil then
         WriteLn('Gene 1 is nil :-(');
      WriteLn('Gene 1 has a value of ', gene1.Print, '.');
      gene2 := TGene.Create;
      WriteLn('Gene 2 has a value of ', gene2.Print, '.');
      Write('They are ');
      if gene1.IsEqualTo(gene2) then
         WriteLn('the same.')
      else
         WriteLn('different.');
      haploid := gene1.Pop;
      Write('gene 1 is a ');
      if haploid = dominant then
         WriteLn('dominant gene.')
      else
         WriteLn('recessive gene.');
      haploid := gene2.Pop;
      Write('gene 2 is a ');
      if haploid = dominant then
         WriteLn('dominant gene.')
      else
         WriteLn('recessive gene.');
      Write('Gene 1 mutated from ', gene1.Print, ' to ');
      count       := 0;
      nbrMutation := 0;
      probability := 0.001;
      repeat
         gene1.Mutate(probability, nbrMutation);
         Inc(count)
      until nbrMutation > 0;
      if count = 1 then
         WriteLn(gene1.Print, ' in ', count, ' attempt.')
      else
         WriteLn(gene1.Print, ' in ', count, ' attempts.');
      Write('This gene was set to dominant ');
      gene1.Put(dominant);
      WriteLn(gene1.Print);
      Write('and this one set to recessive ');
      gene2 := TGene.CreateWithGene(recessive);
      WriteLn(gene2.Print, '.');
   end;

begin
  Test
end.

