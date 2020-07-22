VERSION 5.00


Begin VB.Form Primary_Form 
   Caption         =   "frmPrimaryForm"
   ClientHeight    =   3135
   ClientLeft      =   60
   ClientTop       =   405
   ClientWidth     =   4680
   LinkTopic       =   "Form1"
   ScaleHeight     =   3135
   ScaleWidth      =   4680
   StartUpPosition =   3  'Windows Default
   Begin VB.CommandButton Run 
      Caption         =   "cmdRun"
      Height          =   735
      Left            =   240
      TabIndex        =   0
      Top             =   600
      Width           =   1335
   End
End
Attribute VB_Name = "Primary_Form"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Option Explicit



Private Sub cmdRun_Click()
    Call LayThongTin
    Call Defining_Center
    Call Defining_Radius_Lines
    Dim ii, jj As Integer
    For jj = 1 To cr + 1
        For ii = 1 To cr
            Call Descrerizating_Slices(rl(ii), Matrix_C((jj - 1) * cr + ii, 1), Matrix_C((jj - 1) * cr + ii, 2))
            'Xac dinh thong so hinh hoc can de tinh FS
            'Luc nay se co entry point , exit point va Radius
            
            
            
            
            
            
            
        Next
    Next
End Sub
