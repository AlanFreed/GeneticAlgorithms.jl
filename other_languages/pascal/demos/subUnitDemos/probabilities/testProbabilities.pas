// global directives
{$MODE OBJFPC}    // allows objects, classes, interfaces and exception handling
{$H+}             // use ANSI strings
{$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
{$PIC ON}         // create Position Independent Code - for Unix libraries

program testProbabilities;

uses
   gaCore, gaProbabilities;

   procedure Test;
   var
      i : Integer;
   begin
      WriteLn('This demo tests the probability routines.');
      WriteLn;
      WriteLn('Ten coin flips with even odds:');
      for i := 1 to 10 do
         begin
            if IsHeads(evenOdds) then
               WriteLn('   heads')
            else
               WriteLn('   tails')
         end;
      WriteLn('Twenty random integers over [0, 9]:');
      for i := 1 to 5 do
         WriteLn(RandomInteger(0, 9), '   ', RandomInteger(0, 9), '   ',
                 RandomInteger(0, 9), '   ', RandomInteger(0, 9));
      WriteLn('And twenty random numbers distributed Johnson SB');
      WriteLn('with a mean of 0.5 and a standard deviation of 0.25');
      for i := 1 to 5 do
         WriteLn(RandomProbability(0.5, 0.25):8:4, '   ',
                 RandomProbability(0.5, 0.25):8:4, '   ',
                 RandomProbability(0.5, 0.25):8:4, '   ',
                 RandomProbability(0.5, 0.25):8:4, '   ');
   end;

begin
  Test
end.
