unit demo1form;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, TAGraph, TASeries, Forms, Controls, Graphics,
  Dialogs, gaData, gaModel, gaSolver;

type

  { TForm1 }

  TForm1 = class(TForm)
    Chart1            : TChart;
    Chart1LineSeries1 : TLineSeries;
    Chart1LineSeries2 : TLineSeries;
    Chart1LineSeries3 : TLineSeries;
    SaveDialog1       : TSaveDialog;
    procedure FormCreate(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end;

var
  Form1 : TForm1;

implementation

  {$R *.lfm}

  function NewRVector (len : Integer) : TRVector; external 'Arrays';
  // Creates a real vector indexing from 1 to len whose elements are all 0.0.

  { TForm1 }

  procedure TForm1.FormCreate(Sender: TObject);
    var
      b, m       : Real;
      genome     : TSVector;
      i          : Integer;
      parameters : TRVector;
      x, y       : TRVector;
  begin
    x := NewRVector(2);
    y := NewRVector(2);
    x[1] := xMin;
    x[2] := xMax;
    // plot the experimental data
    for i := 1 to dimS do
      Chart1LineSeries1.AddXY(xData[i], yData[i]);
    // plot the line predicted by the genetic algorithm
    parameters := nil;
    genome     := nil;
    GetEliteData(parameters, genome);
    m := parameters[1];
    b := parameters[2];
    for i := 1 to 2 do begin
      y[i] := m * x[i] + b;
      Chart1LineSeries2.AddXY(x[i], y[i])
    end;
    // plot the line predicted by linear regression
    for i := 1 to 2 do begin
      y[i] := a1 * x[i] + a0;
      Chart1LineSeries3.AddXY(x[i], y[i])
    end;
    // save the graphic to file in BMP, JPEG and PNG formats
    Chart1.CopyToClipboardBitmap;
    Chart1.SaveToBitmapFile('demo1.bmp');        // this file is large
    Chart1.SaveToFile(TJPEGImage, 'demo1.jpg');
    Chart1.SaveToFile(TPortableNetworkGraphic, 'demo1.png')
  end;

end.

