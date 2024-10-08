// global directives
{$MODE OBJFPC}    // allows objects, classes, interfaces and exception handling
{$H+}             // use ANSI strings
{$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
{$PIC ON}         // create Position Independent Code - for Unix libraries

program testCreatures;

  uses
    gaCore, gaCreatures;
    
  var
    nP, nV : Integer;
  
  function NewRVector (len : Integer) : TRVector; external 'Arrays';
  // Creates a real vector indexing from 1 to len whose elements are all 0.0.

  function LenRVector (var v : TRVector) : Integer; external 'Arrays';
  // Returns the length of real vector v.
  
  function MyFixedParams (var variedParams : TRVector) : TRVector;
    var
      fixedPar : TRVector;
      n        : Integer;
  begin
    if nP > nV then
      fixedPar := NewRVector(nP - nV)
    else
      fixedPar := nil;
    if fixedPar <> nil then begin
      for n := 1 to nP - nV do begin
        Write('Fixed parameter ', n, ' has value: ');
        ReadLn(fixedPar[n])
      end;
      WriteLn
    end;   
    MyFixedParams := fixedPar
  end;

  procedure Test;
    var
      alienPar  : TRVector;
      creature1 : TCreature;
      creature2 : TCreature;
      creature3 : TCreature;
      fixedPar  : TVectorFn;
      genome    : TSVector;
      maxFitPar : TRVector;
      minFitPar : TRVector;
      myPar     : TRVector;
      nbrMut    : Integer;
      nbrXov    : Integer;
      n         : Integer;
      sigFig    : Integer;
      varyPar   : TBVector;
      yes       : Char;
  begin
    myPar := nil;
    yes   := 'N';
    WriteLn('Test the creature class');
    WriteLn;
    Write('How many total parameters are there? ');
    ReadLn(nP);
    Write('Of these, how many are to be varied? ');
    ReadLn(nV);
    if nV > nP then begin
      WriteLn('input error: number varied parameters cannot exceed total number of parameters');
      halt
    end;
    alienPar  := nil;
    varyPar   := nil;
    fixedPar  := @MyFixedParams;
    minFitPar := nil;
    maxFitPar := nil;
    alienPar  := NewRVector(nV);
    varyPar   := NewBVector(nP);
    minFitPar := NewRVector(nV);
    maxFitPar := NewRVector(nV);
    WriteLn;
    if nP <> nV then
      for n := 1 to nP do begin
        Write('Is parameter ', n, ' to be varied? (Y or N) ');
        ReadLn(yes);
        if (yes = 'Y') or (yes = 'y') then
          varyPar[n] := True
        else
          varyPar[n] := False
      end
    else
      for n := 1 to nP do
        varyPar[n] := True;
    WriteLn;
    WriteLn('Enter the lower and upper bounds for the varied parameters.');
    for n := 1 to nV do begin
      WriteLn;
      Write('Minimum value for varied parameter ', n, ': ');
      ReadLn(minFitPar[n]);
      Write('Maximum value for varied parameter ', n, ': ');
      ReadLn(maxFitPar[n])
    end;
    WriteLn;
    Write('Enter in the number of significant figures: ');
    ReadLn(sigFig);
    creature1 := TCreature.Create
                 (varyPar, fixedPar, minFitPar, maxFitPar, sigFig);
    creature1.Procreate;
    WriteLn;
    WriteLn('A creature was procreated possessing parameters');
    creature1.GetAllParameters(myPar);
    creature1.GetGenome(genome);
    for n := 1 to LenRVector(myPar) do
      WriteLn('      ', genome[n], ' <=> ', myPar[n]:8:4);
    WriteLn;
    creature2 := TCreature.Create
                 (varyPar, fixedPar, minFitPar, maxFitPar, sigFig);
    for n := 1 to nV do
      alienPar[n] := (maxFitPar[n] + minFitPar[n]) / 2.0;
    creature2.Alien(alienPar);
    WriteLn('Another was created with alien DNA possessing parameters');
    creature2.GetAllParameters(myPar);
    creature2.GetGenome(genome);
    for n := 1 to LenRVector(myPar) do
      WriteLn('      ', genome[n], ' <=> ', myPar[n]:8:4);
    if nV > 1 then
      WriteLn('These should look similar as they are assigned midpoint values.');
    WriteLn;
    Write('These creatures are ');
    if creature1.IsEqualTo(creature2) then
      WriteLn('equal.')
    else
      WriteLn('not equal.');
    WriteLn;
    nbrMut := 0;
    nbrXov := 0;
    creature3 := TCreature.Create
                 (varyPar, fixedPar, minFitPar, maxFitPar, sigFig);
    creature3.Conceive(creature1, creature2, nbrMut, nbrXov);
    WriteLn('Finally a creature was conceived from these two:');
    WriteLn('   culminating in ', nbrMut, ' mutations and ', 
            nbrXov, ' crossover events.');
    WriteLn('   producing parameters:');
    creature3.GetAllParameters(myPar);
    creature3.GetGenome(genome);
    for n := 1 to LenRVector(myPar) do
      WriteLn('      ', genome[n], ' <=> ',  myPar[n]:8:4)
  end;

begin
  Test
end.

